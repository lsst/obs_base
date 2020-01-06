# This file is part of obs_base.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
from __future__ import annotations

__all__ = ["RepoConverter"]

import os
import fnmatch
from dataclasses import dataclass
from collections import defaultdict
from abc import ABC, abstractmethod
from typing import (
    Any,
    Callable,
    Dict,
    Generic,
    Iterator,
    List,
    MutableMapping,
    Optional,
    Set,
    Tuple,
    TYPE_CHECKING,
    TypeVar,
)

from lsst.daf.butler import DatasetRef, Butler as Butler3, DataCoordinate, FileDataset, DatasetType
from lsst.sphgeom import RangeSet, Region

from .filePathParser import FilePathParser

if TYPE_CHECKING:
    from ..mapping import Mapping as CameraMapperMapping  # disambiguate from collections.abc.Mapping
    from .dataIdExtractor import DataIdExtractor
    from .convertRepo import ConvertRepoTask
    from lsst.daf.butler import StorageClass, Registry, SkyPixDimension


REPO_ROOT_FILES = ("registry.sqlite3", "_mapper", "repositoryCfg.yaml", "calibRegistry.sqlite3", "_parent")


T = TypeVar("T")


class MostRecentlyUsedStack(Generic[T]):
    """A simple container that maintains a most-recently-used ordering.
    """

    def __init__(self):
        self._elements = []

    def __iter__(self):
        # Iterate in reverse order so we can keep the most recent element used
        # at the end of the list.  We want to use the end rather than the
        # beginning because appending to lists is much more efficient than
        # inserting at the beginning.
        yield from reversed(self._elements)

    def apply(self, func: Callable[[T], Any]) -> Any:
        """Apply a function to elements until it returns a value that coerces
        to `True`, and move the corresponding element to the front of the
        stack.

        Parameters
        ----------
        func : callable
            Callable object.

        Returns
        -------
        value : `object`
            The first value returned by ``func`` that coerces to `True`.
        """
        for n, element in enumerate(self):
            result = func(element)
            if result:
                break
        else:
            return None
        # Move the extractor that matched to the back of the list (note that
        # n indexes from the back of the internal list).
        if n != 0:
            # i indexes from the front of the internal list.
            i = len(self._elements) - 1 - n
            assert self._elements[i] is element
            del self._elements[i]
            self._elements.append(element)
        return result

    def push(self, element):
        """Add a new element to the front of the stack.
        """
        self._elements.append(element)


@dataclass
class ConversionSubset:
    """A helper class for `ConvertRepoTask` and `RepoConverter` that maintains
    lists of related data ID values that should be included in the conversion.

    Parameters
    ----------
    instrument : `str`
        Instrument name used in Gen3 data IDs.
    visits : `set` of `int`
        Visit IDs that define the filter.
    """

    def __init__(self, instrument: str, visits: Set[int]):
        self.instrument = instrument
        self.visits = visits
        self.regions = None
        self.tracts = {}
        self.skypix = {}

    def addSkyMap(self, registry: Registry, name: str):
        """Populate the included tract IDs for the given skymap from those that
        overlap the visits the `ConversionSubset` was initialized with.

        Parameters
        ----------
        registry : `lsst.daf.butler.Registry`
            Registry that can be queried for visit/tract overlaps.
        name : `str`
            SkyMap name used in Gen3 data IDs.
        """
        tracts = set()
        self.tracts[name] = tracts
        for visit in self.visits:
            for dataId in registry.queryDimensions(["tract"], expand=False,
                                                   dataId={"skymap": name,
                                                           "instrument": self.instrument,
                                                           "visit": visit}):
                tracts.add(dataId["tract"])

    def addSkyPix(self, registry: Registry, dimension: SkyPixDimension):
        """Populate the included skypix IDs for the given dimension from those
        that overlap the visits the `ConversionSubset` was initialized with.

        Parameters
        ----------
        registry : `lsst.daf.butler.Registry`
            Registry that can be queried for visit regions.
        name : `str`
            SkyMap name used in Gen3 data IDs.
        """
        if self.regions is None:
            self.regions = []
            for visit in self.visits:
                dataId = registry.expandDataId(instrument=self.instrument, visit=visit)
                self.regions.append(dataId.region)
        ranges = RangeSet()
        for region in self.regions:
            ranges = ranges.union(dimension.pixelization.envelope(region))
        self.skypix[dimension] = ranges

    def isRelated(self, dataId: DataCoordinate) -> bool:
        """Test whether the given data ID is related to this subset and hence
        should be included in a repository conversion.

        Parameters
        ----------
        dataId : `lsst.daf.butler.DataCoordinate`
            Data ID to test.

        Returns
        -------
        related : `bool`
            `True` if this data ID should be included in a repository
            conversion.

        Notes
        -----
        More formally, this tests that the given data ID is not unrelated;
        if a data ID does not involve tracts, visits, or skypix dimensions,
        we always include it.
        """
        if self.visits is None:
            # We're not filtering at all.
            return True
        if "visit" in dataId.graph and dataId["visit"] not in self.visits:
            return False
        if "tract" in dataId.graph and dataId["tract"] not in self.tracts[dataId["skymap"]]:
            return False
        for dimension, ranges in self.skypix.items():
            if dimension in dataId.graph and not ranges.intersects(dataId[dimension]):
                return False
        return True

    # Class attributes that will be shadowed by public instance attributes;
    # defined here only for documentation purposes.

    instrument: str
    """The name of the instrument, as used in Gen3 data IDs (`str`).
    """

    visits: Set[int]
    """The set of visit IDs that should be included in the conversion (`set`
    of `int`).
    """

    regions: Optional[List[Region]]
    """Regions for all visits (`list` of `lsst.sphgeom.Region`).

    Set to `None` before it has been initialized.  Any code that attempts to
    use it when it is `None` has a logic bug.
    """

    tracts: Dict[str, Set[int]]
    """Tracts that should be included in the conversion, grouped by skymap
    name (`dict` mapping `str` to `set` of `int`).
    """

    skypix: Dict[SkyPixDimension, RangeSet]
    """SkyPix ranges that should be included in the conversion, grouped by
    dimension (`dict` mapping `SkyPixDimension` to `lsst.sphgeom.RangeSet`).
    """


class RepoConverter(ABC):
    """An abstract base class for objects that help `ConvertRepoTask` convert
    datasets from a single Gen2 repository.

    Parameters
    ----------
    task : `ConvertRepoTask`
        Task instance that is using this helper object.
    root : `str`
        Root of the Gen2 repo being converted.
    collections : `list` of `str`
        Gen3 collections with which all converted datasets should be
        associated.
    subset : `ConversionSubset, optional
        Helper object that implements a filter that restricts the data IDs that
        are converted.

    Notes
    -----
    `RepoConverter` defines the only public API users of its subclasses should
    use (`prep`, `insertDimensionRecords`, and `ingest`).  These delegate to
    several abstract methods that subclasses must implement.  In some cases,
    subclasses may reimplement the public methods as well, but are expected to
    delegate to ``super()`` either at the beginning or end of their own
    implementation.
    """

    def __init__(self, *, task: ConvertRepoTask, root: str, collections: List[str],
                 subset: Optional[ConversionSubset] = None):
        self.task = task
        self.root = root
        self.subset = subset
        self._collections = list(collections)
        self._extractors: MostRecentlyUsedStack[DataIdExtractor] = MostRecentlyUsedStack()
        self._skipParsers: MostRecentlyUsedStack[Tuple[FilePathParser, str, str]] = MostRecentlyUsedStack()
        self._fileDatasets: MutableMapping[DatasetType, List[FileDataset]] = defaultdict(list)

    @abstractmethod
    def isDatasetTypeSpecial(self, datasetTypeName: str) -> bool:
        """Test whether the given dataset is handled specially by this
        converter and hence should be ignored by generic base-class logic that
        searches for dataset types to convert.

        Parameters
        ----------
        datasetTypeName : `str`
            Name of the dataset type to test.

        Returns
        -------
        special : `bool`
            `True` if the dataset type is special.
        """
        raise NotImplementedError()

    @abstractmethod
    def isDirectorySpecial(self, subdirectory: str) -> bool:
        """Test whether the given directory is handled specially by this
        converter and hence should be ignored by generic base-class logic that
        searches for datasets to convert.

        Parameters
        ----------
        subdirectory : `str`
            Subdirectory.  This is only ever a single subdirectory, and it
            could appear anywhere within a repo root.  (A full path relative
            to the repo root might be more useful, but it is harder to
            implement, and we don't currently need it to identify any special
            directories).

        Returns
        -------
        special : `bool`
            `True` if the direct is special.
        """
        raise NotImplementedError()

    @abstractmethod
    def iterMappings(self) -> Iterator[Tuple[str, CameraMapperMapping]]:
        """Iterate over all `CameraMapper` `Mapping` objects that should be
        considered for conversion by this repository.

        This this should include any datasets that may appear in the
        repository, including those that are special (see
        `isDatasetTypeSpecial`) and those that are being ignored (see
        `ConvertRepoTask.isDatasetTypeIncluded`); this allows the converter
        to identify and hence skip these datasets quietly instead of warning
        about them as unrecognized.

        Yields
        ------
        datasetTypeName: `str`
            Name of the dataset type.
        mapping : `lsst.obs.base.mapping.Mapping`
            Mapping object used by the Gen2 `CameraMapper` to describe the
            dataset type.
        """
        raise NotImplementedError()

    @abstractmethod
    def makeDataIdExtractor(self, datasetTypeName: str, parser: FilePathParser,
                            storageClass: StorageClass) -> DataIdExtractor:
        """Construct a `DataIdExtractor` instance appropriate for a particular
        dataset type.

        Parameters
        ----------
        datasetTypeName : `str`
            Name of the dataset type; typically forwarded directly to
            the `DataIdExtractor` constructor.
        parser : `FilePathParser`
            Object that parses filenames into Gen2 data IDs; typically
            forwarded directly to the `DataIdExtractor` constructor.
        storageClass : `lsst.daf.butler.StorageClass`
            Storage class for this dataset type in the Gen3 butler; typically
            forwarded directly to the `DataIdExtractor` constructor.

        Returns
        -------
        extractor : `DataIdExtractor`
            A new `DataIdExtractor` instance.
        """
        raise NotImplementedError()

    def prep(self):
        """Perform preparatory work associated with the dataset types to be
        converted from this repository (but not the datasets themselves).

        Notes
        -----
        This should be a relatively fast operation that should not depend on
        the size of the repository.

        Subclasses may override this method, but must delegate to the base
        class implementation at some point in their own logic.
        More often, subclasses will specialize the behavior of
        `prepDatasetTypes` by overriding other methods to which the base class
        implementation delegates.
        These include:
         - `iterMappings`
         - `isDatasetTypeSpecial`
         - `makeDataIdExtractor`

        This should not perform any write operations to the Gen3 repository.
        It is guaranteed to be called before `insertDimensionData`.
        """
        self.task.log.info(f"Preparing other dataset types from root {self.root}.")
        for datasetTypeName, mapping in self.iterMappings():
            try:
                parser = FilePathParser.fromMapping(mapping)
            except RuntimeError:
                # No template, so there should be no way we'd get one of these
                # in the Gen2 repo anyway (and if we do, we'll still produce a
                # warning - just a less informative one than we might be able
                # to produce if we had a template).
                continue
            if (not self.task.isDatasetTypeIncluded(datasetTypeName) or
                    self.isDatasetTypeSpecial(datasetTypeName)):
                # User indicated not to include this data, but we still want
                # to recognize files of that type to avoid warning about them.
                self._skipParsers.push((parser, datasetTypeName, None))
                continue
            storageClass = self._guessStorageClass(datasetTypeName, mapping)
            if storageClass is None:
                # This may be a problem, but only if we actually encounter any
                # files corresponding to this dataset.  Of course, we need
                # to be able to parse those files in order to recognize that
                # situation.
                self._skipParsers.push((parser, datasetTypeName, "no storage class found."))
                continue
            self._extractors.push(self.makeDataIdExtractor(datasetTypeName, parser, storageClass))

    def iterDatasets(self) -> Iterator[FileDataset]:
        """Iterate over all datasets in the repository that should be
        ingested into the Gen3 repository.

        Subclasses may override this method, but must delegate to the base
        class implementation at some point in their own logic.

        Yields
        ------
        dataset : `FileDataset`
            Structures representing datasets to be ingested.  Paths should be
            absolute.
        ref : `lsst.daf.butler.DatasetRef`
            Reference for the Gen3 datasets, including a complete `DatasetType`
            and data ID.
        """
        for dirPath, subdirNamesInDir, fileNamesInDir in os.walk(self.root, followlinks=True):
            # Remove subdirectories that appear to be repositories themselves
            # from the walking
            def isRepoRoot(dirName):
                return any(os.path.exists(os.path.join(dirPath, dirName, f))
                           for f in REPO_ROOT_FILES)
            subdirNamesInDir[:] = [d for d in subdirNamesInDir
                                   if not isRepoRoot(d) and not self.isDirectorySpecial(d)]
            # Loop over files in this directory, and ask per-DatasetType
            # extractors if they recognize them and can extract a data ID;
            # if so, ingest.
            dirPathInRoot = dirPath[len(self.root) + len(os.path.sep):]
            for fileNameInDir in fileNamesInDir:
                if any(fnmatch.fnmatchcase(fileNameInDir, pattern)
                       for pattern in self.task.config.fileIgnorePatterns):
                    continue
                fileNameInRoot = os.path.join(dirPathInRoot, fileNameInDir)
                if fileNameInRoot in REPO_ROOT_FILES:
                    continue
                ref = self._extractDatasetRef(fileNameInRoot)
                if ref is not None:
                    if self.subset is None or self.subset.isRelated(ref.dataId):
                        yield FileDataset(path=os.path.join(self.root, fileNameInRoot), refs=ref)
                else:
                    self._handleUnrecognizedFile(fileNameInRoot)

    def findDatasets(self):
        self.task.log.info("Finding datasets from files in repo %s.", self.root)
        for dataset in self.iterDatasets():
            assert len(dataset.refs) == 1
            self._fileDatasets[dataset.refs[0].datasetType].append(dataset)

    def insertDimensionData(self):
        """Insert any dimension records uniquely derived from this repository
        into the registry.

        Subclasses may override this method, but may not need to; the default
        implementation does nothing.

        SkyMap and SkyPix dimensions should instead be handled by calling
        `ConvertRepoTask.useSkyMap` or `ConvertRepoTask.useSkyPix`, because
        these dimensions are in general shared by multiple Gen2 repositories.

        This method is guaranteed to be called between `prep` and `ingest`.
        """
        pass

    def handleDataIdExpansionFailure(self, dataset: FileDataset, err: LookupError):
        self.task.log.warn("Skipping ingestion for '%s': %s", dataset.path, err)
        return False

    def expandDataIds(self):
        for datasetType, datasetsForType in self._fileDatasets.items():
            self.task.log.info("Expanding data IDs for %s %s datasets.", len(datasetsForType),
                               datasetType.name)
            expanded = []
            for dataset in datasetsForType:
                try:
                    dataId = self.task.registry.expandDataId(dataset.ref.dataId)
                    dataset.ref = dataset.ref.expanded(dataId)
                    expanded.append(dataset)
                except LookupError as err:
                    if self.handleDataIdExpansionFailure(dataset, err):
                        expanded.append(dataset)
            datasetsForType[:] = expanded

    def ingest(self):
        """Insert converted datasets into the Gen3 repository.

        Subclasses may override this method, but must delegate to the base
        class implementation at some point in their own logic.

        This method is guaranteed to be called after both `prep` and
        `insertDimensionData`.
        """
        for datasetType, datasetsForType in self._fileDatasets.items():
            self.task.registry.registerDatasetType(datasetType)
            self.task.log.info("Ingesting %s %s datasets.", len(datasetsForType), datasetType.name)
            try:
                butler3, collections = self.getButler(datasetType.name)
            except LookupError as err:
                self.task.log.warn(str(err))
                continue
            try:
                butler3.ingest(*datasetsForType, transfer=self.task.config.transfer)
            except LookupError as err:
                raise LookupError(f"Error expanding data ID for dataset type {datasetType.name}.") from err
            for collection in collections:
                self.task.registry.associate(collection,
                                             [ref for dataset in datasetsForType for ref in dataset.refs])

    def getButler(self, datasetTypeName: str) -> Tuple[Butler3, List[str]]:
        """Create a new Gen3 Butler appropriate for a particular dataset type.

        This should be used exclusively by subclasses when obtaining a butler
        to use for dataset ingest (`ConvertRepoTask.butler3` should never be
        used directly).

        Parameters
        ----------
        datasetTypeName : `str`
            Name of the dataset type.

        Returns
        -------
        butler : `lsst.daf.butler.Butler`
            Gen3 Butler instance appropriate for ingesting the given dataset
            type.
        collections : `list` of `str`
            Collections the dataset should be associated with, in addition to
            the one used to define the `lsst.daf.butler.Run` used in
            ``butler``.
        """
        if datasetTypeName in self.task.config.collections:
            return (
                Butler3(butler=self.task.butler3, run=self.task.config.collections[datasetTypeName]),
                self._collections,
            )
        elif self._collections:
            return (
                Butler3(butler=self.task.butler3, run=self._collections[0]),
                self._collections[1:],
            )
        else:
            raise LookupError("No collection configured for dataset type {datasetTypeName}.")

    def _extractDatasetRef(self, fileNameInRoot: str) -> Optional[DatasetRef]:
        """Extract a `DatasetRef` from a file name.

        This method is for internal use by `RepoConverter` itself (not its
        subclasses).

        Parameters
        ----------
        fileNameInRoot : `str`
            Name of the file to be ingested, relative to the repository root.

        Returns
        -------
        ref : `lsst.daf.butler.DatasetRef` or `None`
            Reference for the Gen3 datasets, including a complete `DatasetType`
            and data ID.  `None` if the converter does not recognize the
            file as one to be converted.
        """
        def closure(extractor):
            try:
                dataId = extractor.apply(fileNameInRoot)
            except LookupError as err:
                raise RuntimeError(f"Error extracting data ID for {extractor.datasetType.name} "
                                   f"on file {fileNameInRoot}.") from err
            if dataId is None:
                return None
            else:
                return DatasetRef(extractor.datasetType, dataId=dataId)
        return self._extractors.apply(closure)

    def _handleUnrecognizedFile(self, fileNameInRoot: str):
        """Generate appropriate warnings (or not) for files not matched by
        `_extractDatasetRef`.

        This method is for internal use by `RepoConverter` itself (not its
        subclasses).

        Parameters
        ----------
        fileNameInRoot : `str`
            Name of the file, relative to the repository root.
        """
        def closure(skipTuple):
            parser, datasetTypeName, message = skipTuple
            if parser(fileNameInRoot) is not None:
                if message is not None:
                    self.task.log.warn("Skipping dataset %s file %s: %s", datasetTypeName,
                                       fileNameInRoot, message)
                return True
            return False
        if not self._skipParsers.apply(closure):
            self.task.log.warn("Skipping unrecognized file %s.", fileNameInRoot)

    def _guessStorageClass(self, datasetTypeName: str, mapping: CameraMapperMapping
                           ) -> Optional[StorageClass]:
        """Infer the Gen3 `StorageClass` from a dataset from a combination of
        configuration and Gen2 dataset type information.

        datasetTypeName: `str`
            Name of the dataset type.
        mapping : `lsst.obs.base.mapping.Mapping`
            Mapping object used by the Gen2 `CameraMapper` to describe the
            dataset type.
        """
        storageClassName = self.task.config.storageClasses.get(datasetTypeName)
        if storageClassName is None and mapping.python is not None:
            storageClassName = self.task.config.storageClasses.get(mapping.python, None)
        if storageClassName is None and mapping.persistable is not None:
            storageClassName = self.task.config.storageClasses.get(mapping.persistable, None)
        if storageClassName is None and mapping.python is not None:
            unqualified = mapping.python.split(".")[-1]
            storageClassName = self.task.config.storageClasses.get(unqualified, None)
        if storageClassName is not None:
            storageClass = self.task.butler3.storageClasses.getStorageClass(storageClassName)
        else:
            try:
                storageClass = self.task.butler3.storageClasses.getStorageClass(mapping.persistable)
            except KeyError:
                storageClass = None
            if storageClass is None and mapping.python is not None:
                try:
                    storageClass = self.task.butler3.storageClasses.getStorageClass(unqualified)
                except KeyError:
                    pass
            if storageClass is None:
                self.task.log.debug("No StorageClass found for %s; skipping.", datasetTypeName)
            else:
                self.task.log.debug("Using StorageClass %s for %s.", storageClass.name, datasetTypeName)
        return storageClass

    # Class attributes that will be shadowed by public instance attributes;
    # defined here only for documentation purposes.

    task: ConvertRepoTask
    """The parent task that constructed and uses this converter
    (`ConvertRepoTask`).
    """

    root: str
    """Root path to the Gen2 repository this converter manages (`str`).

    This is a complete path, not relative to some other repository root.
    """

    subset: Optional[ConversionSubset]
    """An object that represents a filter to be applied to the datasets that
    are converted (`ConversionSubset` or `None`).
    """
