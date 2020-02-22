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

from dataclasses import dataclass
from collections import defaultdict
from abc import ABC, abstractmethod
import fnmatch
import re
from typing import (
    Dict,
    Iterator,
    List,
    MutableMapping,
    Optional,
    Set,
    Tuple,
    Union,
    TYPE_CHECKING,
)

from lsst.utils import doImport
from lsst.daf.butler import DataCoordinate, FileDataset, DatasetType
from lsst.sphgeom import RangeSet, Region
from .repoWalker import RepoWalker

if TYPE_CHECKING:
    from ..mapping import Mapping as CameraMapperMapping  # disambiguate from collections.abc.Mapping
    from .convertRepo import ConvertRepoTask
    from .scanner import PathElementHandler
    from lsst.daf.butler import StorageClass, Registry, SkyPixDimension, FormatterParameter


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

    def __init__(self, *, task: ConvertRepoTask, root: str, run: Optional[str],
                 subset: Optional[ConversionSubset] = None):
        self.task = task
        self.root = root
        self.subset = subset
        self._run = run
        self._repoWalker = None  # Created in prep
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
    def makeRepoWalkerTarget(self, datasetTypeName: str, template: str, keys: Dict[str, type],
                             storageClass: StorageClass,
                             formatter: FormatterParameter = None,
                             targetHandler: Optional[PathElementHandler] = None,
                             ) -> RepoWalker.Target:
        """Make a struct that identifies a dataset type to be extracted by
        walking the repo directory structure.

        Parameters
        ----------
        datasetTypeName : `str`
            Name of the dataset type (the same in both Gen2 and Gen3).
        template : `str`
            The full Gen2 filename template.
        keys : `dict` [`str`, `type`]
            A dictionary mapping Gen2 data ID key to the type of its value.
        storageClass : `lsst.daf.butler.StorageClass`
            Gen3 storage class for this dataset type.
        formatter : `lsst.daf.butler.Formatter` or `str`, optional
            A Gen 3 formatter class or fully-qualified name.
        targetHandler : `PathElementHandler`, optional
            Specialist target handler to use for this dataset type.

        Returns
        -------
        target : `RepoWalker.Target`
            A struct containing information about the target dataset (much of
            it simplify forwarded from the arguments).
        """
        raise NotImplementedError()

    def getSpecialDirectories(self) -> List[str]:
        """Return a list of directory paths that should not be searched for
        files.

        These may be directories that simply do not contain datasets (or
        contain datasets in another repository), or directories whose datasets
        are handled specially by a subclass.

        Returns
        -------
        directories : `list` [`str`]
            The full paths of directories to skip, relative to the repository
            root.
        """
        return []

    def prep(self):
        """Perform preparatory work associated with the dataset types to be
        converted from this repository (but not the datasets themselves).

        Notes
        -----
        This should be a relatively fast operation that should not depend on
        the size of the repository.

        Subclasses may override this method, but must delegate to the base
        class implementation at some point in their own logic.
        More often, subclasses will specialize the behavior of `prep` by
        overriding other methods to which the base class implementation
        delegates.  These include:
         - `iterMappings`
         - `isDatasetTypeSpecial`
         - `getSpecialDirectories`
         - `makeRepoWalkerTarget`

        This should not perform any write operations to the Gen3 repository.
        It is guaranteed to be called before `insertDimensionData`.
        """
        self.task.log.info(f"Preparing other dataset types from root {self.root}.")
        walkerInputs: List[Union[RepoWalker.Target, RepoWalker.Skip]] = []
        for datasetTypeName, mapping in self.iterMappings():
            try:
                template = mapping.template
            except RuntimeError:
                # No template for this dataset in this mapper, so there's no
                # way there should be instances of this dataset in this repo.
                continue
            extensions = [""]
            skip = False
            message = None
            storageClass = None
            if (not self.task.isDatasetTypeIncluded(datasetTypeName)
                    or self.isDatasetTypeSpecial(datasetTypeName)):
                # User indicated not to include this data, but we still want
                # to recognize files of that type to avoid warning about them.
                skip = True
            else:
                storageClass = self._guessStorageClass(datasetTypeName, mapping)
                if storageClass is None:
                    # This may be a problem, but only if we actually encounter any
                    # files corresponding to this dataset.  Of course, we need
                    # to be able to parse those files in order to recognize that
                    # situation.
                    message = f"no storage class found for {datasetTypeName}"
                    skip = True
            # Handle files that are compressed on disk, but the gen2 template is just `.fits`
            if template.endswith(".fits"):
                extensions.extend((".gz", ".fz"))
            for extension in extensions:
                if skip:
                    walkerInput = RepoWalker.Skip(
                        template=template+extension,
                        keys=mapping.keys(),
                        message=message,
                    )
                    self.task.log.debug("Skipping template in walker: %s", template)
                else:
                    assert message is None
                    targetHandler = self.task.config.targetHandlerClasses.get(datasetTypeName)
                    if targetHandler is not None:
                        targetHandler = doImport(targetHandler)
                    walkerInput = self.makeRepoWalkerTarget(
                        datasetTypeName=datasetTypeName,
                        template=template+extension,
                        keys=mapping.keys(),
                        storageClass=storageClass,
                        formatter=self.task.config.formatterClasses.get(datasetTypeName),
                        targetHandler=targetHandler,
                    )
                    self.task.log.debug("Adding template to walker: %s", template)
                walkerInputs.append(walkerInput)

        for dirPath in self.getSpecialDirectories():
            walkerInputs.append(
                RepoWalker.Skip(
                    template=dirPath,  # not really a template, but that's fine; it's relative to root.
                    keys={},
                    message=None,
                    isForFiles=True,
                )
            )
        fileIgnoreRegExTerms = []
        for pattern in self.task.config.fileIgnorePatterns:
            fileIgnoreRegExTerms.append(fnmatch.translate(pattern))
        if fileIgnoreRegExTerms:
            fileIgnoreRegEx = re.compile("|".join(fileIgnoreRegExTerms))
        else:
            fileIgnoreRegEx = None
        self._repoWalker = RepoWalker(walkerInputs, fileIgnoreRegEx=fileIgnoreRegEx)

    def iterDatasets(self) -> Iterator[FileDataset]:
        """Iterate over datasets in the repository that should be ingested into
        the Gen3 repository.

        The base class implementation yields nothing; the datasets handled by
        the `RepoConverter` base class itself are read directly in
        `findDatasets`.

        Subclasses should override this method if they support additional
        datasets that are handled some other way.

        Yields
        ------
        dataset : `FileDataset`
            Structures representing datasets to be ingested.  Paths should be
            absolute.
        """
        yield from ()

    def findDatasets(self):
        assert self._repoWalker, "prep() must be called before findDatasets."
        self.task.log.info("Adding special datasets in repo %s.", self.root)
        for dataset in self.iterDatasets():
            assert len(dataset.refs) == 1
            self._fileDatasets[dataset.refs[0].datasetType].append(dataset)
        self.task.log.info("Finding datasets from files in repo %s.", self.root)
        self._fileDatasets.update(
            self._repoWalker.walk(
                self.root,
                log=self.task.log,
                predicate=(self.subset.isRelated if self.subset is not None else None)
            )
        )

    def insertDimensionData(self):
        """Insert any dimension records uniquely derived from this repository
        into the registry.

        Subclasses may override this method, but may not need to; the default
        implementation does nothing.

        SkyMap and SkyPix dimensions should instead be handled by calling
        `ConvertRepoTask.useSkyMap` or `ConvertRepoTask.useSkyPix`, because
        these dimensions are in general shared by multiple Gen2 repositories.

        This method is guaranteed to be called between `prep` and
        `expandDataIds`.
        """
        pass

    def expandDataIds(self):
        """Expand the data IDs for all datasets to be inserted.

        Subclasses may override this method, but must delegate to the base
        class implementation if they do.

        This involves queries to the registry, but not writes.  It is
        guaranteed to be called between `insertDimensionData` and `ingest`.
        """
        import itertools
        for datasetType, datasetsForType in self._fileDatasets.items():
            self.task.log.info("Expanding data IDs for %s %s datasets.", len(datasetsForType),
                               datasetType.name)
            expanded = []
            for dataset in datasetsForType:
                for i, ref in enumerate(dataset.refs):
                    try:
                        dataId = self.task.registry.expandDataId(ref.dataId)
                        dataset.refs[i] = ref.expanded(dataId)
                    except LookupError as err:
                        self.task.log.warn("Skipping ingestion for '%s': %s", dataset.path, err)
                        # Remove skipped datasets from multi-extension FileDatasets
                        dataset.refs[i] = None  # We will strip off the `None`s after the loop.
                dataset.refs[:] = itertools.filterfalse(lambda x: x is None, dataset.refs)
                if dataset.refs:
                    expanded.append(dataset)

            datasetsForType[:] = expanded

    def ingest(self):
        """Insert converted datasets into the Gen3 repository.

        Subclasses may override this method, but must delegate to the base
        class implementation at some point in their own logic.

        This method is guaranteed to be called after `expandDataIds`.
        """
        for datasetType, datasetsForType in self._fileDatasets.items():
            self.task.registry.registerDatasetType(datasetType)
            self.task.log.info("Ingesting %s %s datasets.", len(datasetsForType), datasetType.name)
            try:
                run = self.getRun(datasetType.name)
            except LookupError:
                self.task.log.warn(f"No run configured for dataset type {datasetType.name}.")
                continue
            try:
                self.task.registry.registerRun(run)
                self.task.butler3.ingest(*datasetsForType, transfer=self.task.config.transfer, run=run)
            except LookupError as err:
                raise LookupError(f"Error expanding data ID for dataset type {datasetType.name}.") from err

    def getRun(self, datasetTypeName: str) -> str:
        """Return the name of the run to insert instances of the given dataset
        type into in this collection.

        Parameters
        ----------
        datasetTypeName : `str`
            Name of the dataset type.

        Returns
        -------
        run : `str`
            Name of the `~lsst.daf.butler.CollectionType.RUN` collection.
        """
        assert self._run is not None, "Method must be overridden if self._run is allowed to be None"
        return self._run

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
