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

__all__ = ["ConvertRepoConfig", "ConvertRepoTask", "ConvertRepoSkyMapConfig", "Rerun"]

import os
import fnmatch
from dataclasses import dataclass
from typing import Iterable, Optional, List, Dict

from lsst.utils import doImport
from lsst.daf.butler import (
    Butler as Butler3,
    CollectionType,
    SkyPixDimension
)
from lsst.pex.config import Config, ConfigurableField, ConfigDictField, DictField, ListField, Field
from lsst.pipe.base import Task
from lsst.skymap import skyMapRegistry, BaseSkyMap

from ..ingest import RawIngestTask
from ..defineVisits import DefineVisitsTask
from .repoConverter import ConversionSubset
from .rootRepoConverter import RootRepoConverter
from .calibRepoConverter import CalibRepoConverter
from .standardRepoConverter import StandardRepoConverter


@dataclass
class ConfiguredSkyMap:
    """Struct containing information about a skymap that may appear in a Gen2
    repository.
    """

    name: str
    """Name of the skymap used in Gen3 data IDs.
    """

    sha1: bytes
    """Hash computed by `BaseSkyMap.getSha1`.
    """

    instance: BaseSkyMap
    """Name of the skymap used in Gen3 data IDs.
    """

    used: bool = False
    """Whether this skymap has been found in at least one repository being
    converted.
    """


@dataclass
class Rerun:
    """Specification for a Gen2 processing-output repository to convert.
    """

    path: str
    """Absolute or relative (to the root repository) path to the Gen2
    repository (`str`).
    """

    runName: str
    """Name of the `~lsst.daf.butler.CollectionType.RUN` collection datasets
    will be inserted into (`str`).
    """

    chainName: Optional[str]
    """Name of a `~lsst.daf.butler.CollectionType.CHAINED` collection that will
    combine this repository's datasets with those of its parent repositories
    (`str`, optional).
    """

    parents: List[str]
    """Collection names associated with parent repositories, used to define the
    chained collection (`list` [ `str` ]).

    Ignored if `chainName` is `None`.  Runs used in the root repo are
    automatically included.
    """


class ConvertRepoSkyMapConfig(Config):
    """Sub-config used to hold the parameters of a SkyMap.

    Notes
    -----
    This config only needs to exist because we can't put a
    `~lsst.pex.config.RegistryField` directly inside a
    `~lsst.pex.config.ConfigDictField`.

    It needs to have its only field named "skyMap" for compatibility with the
    configuration of `lsst.pipe.tasks.MakeSkyMapTask`, which we want so we can
    use one config file in an obs package to configure both.

    This name leads to unfortunate repetition with the field named
    "skymap" that holds it - "skyMap[name].skyMap" - but that seems
    unavoidable.
    """
    skyMap = skyMapRegistry.makeField(
        doc="Type and parameters for the SkyMap itself.",
        default="dodeca",
    )


class ConvertRepoConfig(Config):
    raws = ConfigurableField(
        "Configuration for subtask responsible for ingesting raws and adding "
        "exposure dimension entries.",
        target=RawIngestTask,
    )
    defineVisits = ConfigurableField(
        "Configuration for the subtask responsible for defining visits from "
        "exposures.",
        target=DefineVisitsTask,
    )
    skyMaps = ConfigDictField(
        "Mapping from Gen3 skymap name to the parameters used to construct a "
        "BaseSkyMap instance.  This will be used to associate names with "
        "existing skymaps found in the Gen2 repo.",
        keytype=str,
        itemtype=ConvertRepoSkyMapConfig,
        default={}
    )
    rootSkyMapName = Field(
        "Name of a Gen3 skymap (an entry in ``self.skyMaps``) to assume for "
        "datasets in the root repository when no SkyMap is found there. ",
        dtype=str,
        optional=True,
        default=None,
    )
    runs = DictField(
        "A mapping from dataset type name to the RUN collection they should "
        "be inserted into.  This must include all datasets that can be found "
        "in the root repository; other repositories will use per-repository "
        "runs.",
        keytype=str,
        itemtype=str,
        default={
            "deepCoadd_skyMap": "skymaps",
            "brightObjectMask": "masks",
        }
    )
    storageClasses = DictField(
        "Mapping from dataset type name or Gen2 policy entry (e.g. 'python' "
        "or 'persistable') to the Gen3 StorageClass name.",
        keytype=str,
        itemtype=str,
        default={
            "bias": "ExposureF",
            "dark": "ExposureF",
            "flat": "ExposureF",
            "defects": "Defects",
            "BaseSkyMap": "SkyMap",
            "BaseCatalog": "Catalog",
            "BackgroundList": "Background",
            "raw": "Exposure",
            "MultilevelParquetTable": "DataFrame",
            "ParquetTable": "DataFrame",
            "SkyWcs": "Wcs",
        }
    )
    formatterClasses = DictField(
        "Mapping from dataset type name to formatter class. "
        "By default these are derived from the formatters listed in the"
        " Gen3 datastore configuration.",
        keytype=str,
        itemtype=str,
        default={}
    )
    targetHandlerClasses = DictField(
        "Mapping from dataset type name to target handler class.",
        keytype=str,
        itemtype=str,
        default={}
    )
    doRegisterInstrument = Field(
        "If True (default), add dimension records for the Instrument and its "
        "filters and detectors to the registry instead of assuming they are "
        "already present.",
        dtype=bool,
        default=True,
    )
    doWriteCuratedCalibrations = Field(
        "If True (default), ingest human-curated calibrations directly via "
        "the Instrument interface.  Note that these calibrations are never "
        "converted from Gen2 repositories.",
        dtype=bool,
        default=True,
    )
    refCats = ListField(
        "The names of reference catalogs (subdirectories under ref_cats) to "
        "be converted",
        dtype=str,
        default=[]
    )
    fileIgnorePatterns = ListField(
        "Filename globs that should be ignored instead of being treated as "
        "datasets.",
        dtype=str,
        default=["README.txt", "*~?", "butler.yaml", "gen3.sqlite3",
                 "registry.sqlite3", "calibRegistry.sqlite3", "_mapper",
                 "_parent", "repositoryCfg.yaml"]
    )
    rawDatasetType = Field(
        "Gen2 dataset type to use for raw data.",
        dtype=str,
        default="raw",
    )
    datasetIncludePatterns = ListField(
        "Glob-style patterns for dataset type names that should be converted.",
        dtype=str,
        default=["*"]
    )
    datasetIgnorePatterns = ListField(
        "Glob-style patterns for dataset type names that should not be "
        "converted despite matching a pattern in datasetIncludePatterns.",
        dtype=str,
        default=[]
    )
    ccdKey = Field(
        "Key used for the Gen2 equivalent of 'detector' in data IDs.",
        dtype=str,
        default="ccd",
    )
    relatedOnly = Field(
        "If True (default), only convert datasets that are related to the "
        "ingested visits.  Ignored unless a list of visits is passed to "
        "run().",
        dtype=bool,
        default=False,
    )
    curatedCalibrations = ListField(
        "Dataset types that are handled by `Instrument.writeCuratedCalibrations()` "
        "and thus should not be converted using the standard calibration "
        "conversion system.",
        dtype=str,
        default=["camera",
                 "transmission_sensor",
                 "transmission_filter",
                 "transmission_optics",
                 "transmission_atmosphere",
                 "bfKernel"]
    )
    instrument = Field(
        doc=("Fully-qualified Python name of the `Instrument` subclass for "
             "all converted datasets."),
        dtype=str,
        optional=False,
        default=None,
    )

    @property
    def transfer(self):
        return self.raws.transfer

    @transfer.setter
    def transfer(self, value):
        self.raws.transfer = value

    def setDefaults(self):
        self.transfer = None

    # TODO: check that there are no collection overrides for curated
    # calibrations, since we don't have a good way to utilize them.


class ConvertRepoTask(Task):
    """A task that converts one or more related Gen2 data repositories to a
    single Gen3 data repository (with multiple collections).

    Parameters
    ----------
    config: `ConvertRepoConfig`
        Configuration for this task.
    butler3: `lsst.daf.butler.Butler`
        A writeable Gen3 Butler instance that represents the data repository
        that datasets will be ingested into.  If the 'raw' dataset is
        configured to be included in the conversion, ``butler3.run`` should be
        set to the name of the collection raws should be ingested into, and
        ``butler3.collections`` should include a calibration collection from
        which the ``camera`` dataset can be loaded, unless a calibration repo
        is converted and ``doWriteCuratedCalibrations`` is `True`.
    **kwargs
        Other keyword arguments are forwarded to the `Task` constructor.

    Notes
    -----
    Most of the work of converting repositories is delegated to instances of
    the `RepoConverter` hierarchy.  The `ConvertRepoTask` instance itself holds
    only state that is relevant for all Gen2 repositories being ingested, while
    each `RepoConverter` instance holds only state relevant for the conversion
    of a single Gen2 repository.  Both the task and the `RepoConverter`
    instances are single use; `ConvertRepoTask.run` and most `RepoConverter`
    methods may only be called once on a particular instance.
    """

    ConfigClass = ConvertRepoConfig

    _DefaultName = "convertRepo"

    def __init__(self, config=None, *, butler3: Butler3, **kwargs):
        config.validate()  # Not a CmdlineTask nor PipelineTask, so have to validate the config here.
        super().__init__(config, **kwargs)
        self.butler3 = butler3
        self.registry = self.butler3.registry
        self.universe = self.registry.dimensions
        if self.isDatasetTypeIncluded("raw"):
            self.makeSubtask("raws", butler=butler3)
            self.makeSubtask("defineVisits", butler=butler3)
        else:
            self.raws = None
            self.defineVisits = None
        self.instrument = doImport(self.config.instrument)()
        self._configuredSkyMapsBySha1 = {}
        self._configuredSkyMapsByName = {}
        for name, config in self.config.skyMaps.items():
            instance = config.skyMap.apply()
            self._populateSkyMapDicts(name, instance)
        self._usedSkyPix = set()
        self.translatorFactory = self.instrument.makeDataIdTranslatorFactory()
        self.translatorFactory.log = self.log.getChild("translators")

    def _populateSkyMapDicts(self, name, instance):
        struct = ConfiguredSkyMap(name=name, sha1=instance.getSha1(), instance=instance)
        self._configuredSkyMapsBySha1[struct.sha1] = struct
        self._configuredSkyMapsByName[struct.name] = struct

    def isDatasetTypeIncluded(self, datasetTypeName: str):
        """Return `True` if configuration indicates that the given dataset type
        should be converted.

        This method is intended to be called primarily by the
        `RepoConverter` instances used interally by the task.

        Parameters
        ----------
        datasetTypeName: str
            Name of the dataset type.

        Returns
        -------
        included : `bool`
            Whether the dataset should be included in the conversion.
        """
        return (
            any(fnmatch.fnmatchcase(datasetTypeName, pattern)
                for pattern in self.config.datasetIncludePatterns)
            and not any(fnmatch.fnmatchcase(datasetTypeName, pattern)
                        for pattern in self.config.datasetIgnorePatterns)
        )

    def useSkyMap(self, skyMap: BaseSkyMap, skyMapName: str) -> str:
        """Indicate that a repository uses the given SkyMap.

        This method is intended to be called primarily by the
        `RepoConverter` instances used interally by the task.

        Parameters
        ----------
        skyMap : `lsst.skymap.BaseSkyMap`
            SkyMap instance being used, typically retrieved from a Gen2
            data repository.
        skyMapName : `str`
            The name of the gen2 skymap, for error reporting.

        Returns
        -------
        name : `str`
            The name of the skymap in Gen3 data IDs.

        Raises
        ------
            LookupError
                Raised if the specified skymap cannot be found.
        """
        sha1 = skyMap.getSha1()
        if sha1 not in self._configuredSkyMapsBySha1:
            self._populateSkyMapDicts(skyMapName, skyMap)
        try:
            struct = self._configuredSkyMapsBySha1[sha1]
        except KeyError as err:
            msg = f"SkyMap '{skyMapName}' with sha1={sha1} not included in configuration."
            raise LookupError(msg) from err
        struct.used = True
        return struct.name

    def registerUsedSkyMaps(self, subset: Optional[ConversionSubset]):
        """Register all skymaps that have been marked as used.

        This method is intended to be called primarily by the
        `RepoConverter` instances used interally by the task.

        Parameters
        ----------
        subset : `ConversionSubset`, optional
            Object that will be used to filter converted datasets by data ID.
            If given, it will be updated with the tracts of this skymap that
            overlap the visits in the subset.
        """
        for struct in self._configuredSkyMapsBySha1.values():
            if struct.used:
                struct.instance.register(struct.name, self.registry)
                if subset is not None and self.config.relatedOnly:
                    subset.addSkyMap(self.registry, struct.name)

    def useSkyPix(self, dimension: SkyPixDimension):
        """Indicate that a repository uses the given SkyPix dimension.

        This method is intended to be called primarily by the
        `RepoConverter` instances used interally by the task.

        Parameters
        ----------
        dimension : `lsst.daf.butler.SkyPixDimension`
            Dimension represening a pixelization of the sky.
        """
        self._usedSkyPix.add(dimension)

    def registerUsedSkyPix(self, subset: Optional[ConversionSubset]):
        """Register all skymaps that have been marked as used.

        This method is intended to be called primarily by the
        `RepoConverter` instances used interally by the task.

        Parameters
        ----------
        subset : `ConversionSubset`, optional
            Object that will be used to filter converted datasets by data ID.
            If given, it will be updated with the pixelization IDs that
            overlap the visits in the subset.
        """
        if subset is not None and self.config.relatedOnly:
            for dimension in self._usedSkyPix:
                subset.addSkyPix(self.registry, dimension)

    def run(self, root: str, *,
            calibs: Dict[str, str] = None,
            reruns: List[Rerun],
            visits: Optional[Iterable[int]] = None):
        """Convert a group of related data repositories.

        Parameters
        ----------
        root : `str`
            Complete path to the root Gen2 data repository.  This should be
            a data repository that includes a Gen2 registry and any raw files
            and/or reference catalogs.
        calibs : `dict`
            Dictionary mapping calibration repository path to the
            `~lsst.daf.butler.CollectionType.RUN` collection that converted
            datasets within it should be inserted into.
        reruns : `list` of `Rerun`
            Specifications for rerun (processing output) collections to
            convert.
        visits : iterable of `int`, optional
            The integer IDs of visits to convert.  If not provided, all visits
            in the Gen2 root repository will be converted.
        """
        if calibs is None:
            calibs = {}
        if visits is not None:
            subset = ConversionSubset(instrument=self.instrument.getName(), visits=frozenset(visits))
        else:
            if self.config.relatedOnly:
                self.log.warn("config.relatedOnly is True but all visits are being ingested; "
                              "no filtering will be done.")
            subset = None

        # Make converters for all Gen2 repos.
        converters = []
        rootConverter = RootRepoConverter(task=self, root=root, subset=subset)
        converters.append(rootConverter)
        for calibRoot, run in calibs.items():
            if not os.path.isabs(calibRoot):
                calibRoot = os.path.join(rootConverter.root, calibRoot)
            converter = CalibRepoConverter(task=self, root=calibRoot, run=run,
                                           mapper=rootConverter.mapper,
                                           subset=rootConverter.subset)
            converters.append(converter)
        for spec in reruns:
            runRoot = spec.path
            if not os.path.isabs(runRoot):
                runRoot = os.path.join(rootConverter.root, runRoot)
            converter = StandardRepoConverter(task=self, root=runRoot, run=spec.runName,
                                              subset=rootConverter.subset)
            converters.append(converter)

        # Register the instrument if we're configured to do so.
        if self.config.doRegisterInstrument:
            # Allow registration to fail on the assumption that this means
            # we are reusing a butler
            try:
                self.instrument.register(self.registry)
            except Exception:
                pass

        # Run raw ingest (does nothing if we weren't configured to convert the
        # 'raw' dataset type).
        rootConverter.runRawIngest()

        # Write curated calibrations to all calibration repositories.
        # Add new collections to the list of collections the butler was
        # initialized to pass to DefineVisitsTask, to deal with the (likely)
        # case the only 'camera' dataset in the repo will be one we're adding
        # here.
        if self.config.doWriteCuratedCalibrations:
            for run in calibs.values():
                butler3 = Butler3(butler=self.butler3, run=run)
                self.instrument.writeCuratedCalibrations(butler3)

        # Define visits (also does nothing if we weren't configurd to convert
        # the 'raw' dataset type).
        rootConverter.runDefineVisits()

        # Walk Gen2 repos to find datasets convert.
        for converter in converters:
            converter.prep()

        # Insert dimensions needed by any converters.  In practice this is just
        # calibration_labels right now, because exposures and visits (and
        # things related to them) are handled by RawIngestTask and
        # DefineVisitsTask earlier and skymaps are handled later.
        #
        # Note that we do not try to filter dimensions down to just those
        # related to the given visits, even if config.relatedOnly is True; we
        # need them in the Gen3 repo in order to be able to know which datasets
        # to convert, because Gen2 alone doesn't know enough about the
        # relationships between data IDs.
        for converter in converters:
            converter.insertDimensionData()

        # Insert dimensions that are potentially shared by all Gen2
        # repositories (and are hence managed directly by the Task, rather
        # than a converter instance).
        # This also finishes setting up the (shared) converter.subsets object
        # that is used to filter data IDs for config.relatedOnly.
        self.registerUsedSkyMaps(rootConverter.subset)
        self.registerUsedSkyPix(rootConverter.subset)

        # Look for datasets, generally by scanning the filesystem.
        # This requires dimensions to have already been inserted so we can use
        # dimension information to identify related datasets.
        for converter in converters:
            converter.findDatasets()

        # Expand data IDs.
        for converter in converters:
            converter.expandDataIds()

        # Actually ingest datasets.
        for converter in converters:
            converter.ingest()

        # Add chained collections for reruns.
        for spec in reruns:
            if spec.chainName is not None:
                self.butler3.registry.registerCollection(spec.chainName, type=CollectionType.CHAINED)
                chain = [spec.runName]
                chain.extend(spec.parents)
                chain.extend(rootConverter.getCollectionChain())
                self.log.info("Defining %s from chain %s.", spec.chainName, chain)
                self.butler3.registry.setCollectionChain(spec.chainName, chain)
