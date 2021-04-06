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

__all__ = ["CalibRepo", "ConvertRepoConfig", "ConvertRepoTask", "ConvertRepoSkyMapConfig", "Rerun"]

import os
import fnmatch
from dataclasses import dataclass
from multiprocessing import Pool
from typing import Iterable, Optional, List, Tuple

from lsst.daf.butler import (
    Butler as Butler3,
    ButlerURI,
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
from .._instrument import Instrument


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


def _dropPrefix(s: str, prefix: str) -> Tuple[str, bool]:
    """If ``s`` starts with ``prefix``, return the rest of ``s`` and `True`.
    Otherwise return ``s`` and `False`.
    """
    if s.startswith(prefix):
        return s[len(prefix):], True
    return s, False


@dataclass
class Rerun:
    """Specification for a Gen2 processing-output repository to convert.
    """

    path: str
    """Absolute or relative (to the root repository) path to the Gen2
    repository (`str`).
    """

    runName: Optional[str]
    """Name of the `~lsst.daf.butler.CollectionType.RUN` collection datasets
    will be inserted into (`str` or `None`).

    If `None`, a name will be guessed by calling `guessCollectionNames`.
    """

    chainName: Optional[str]
    """Name of a `~lsst.daf.butler.CollectionType.CHAINED` collection that will
    combine this repository's datasets with those of its parent repositories
    (`str` or `None`).

    If `None`, a name will be guessed by calling `guessCollectionNames`.
    """

    parents: List[str]
    """Collection names associated with parent repositories, used to define the
    chained collection (`list` [ `str` ]).

    Ignored if `chainName` is `None`.  Runs used in the root repo are
    automatically included.
    """

    def guessCollectionNames(self, instrument: Instrument, root: str) -> None:
        """Update `runName` and `chainName` with guesses that match Gen3 naming
        conventions.

        If `chainName` is not `None`, and `runName` is, `runName` will be set
        from it.  If `runName` is already set, nothing will be changed, and
        if `chainName` is `None`, no chained collection will be created.

        Parameters
        ----------
        instrument : `Instrument`
            Instrument object for the repository being converted.
        root : `str`
            Path to the root repository.  If this is present at the start of
            ``self.path``, it will be stripped as part of generating the run
            name.

        Raises
        ------
        ValueError
            Raised if the appropriate collection names cannot be inferred.
        """
        if self.runName is not None:
            return
        if self.chainName is None:
            if os.path.isabs(self.path):
                rerunURI = ButlerURI(self.path)
                rootURI = ButlerURI(root)
                chainName = rerunURI.relative_to(rootURI)
                if chainName is None:
                    raise ValueError(
                        f"Cannot guess run name collection for rerun at '{self.path}': "
                        f"no clear relationship to root '{root}'."
                    )
            else:
                chainName = self.path
            chainName, _ = _dropPrefix(chainName, "rerun/")
            chainName, isPersonal = _dropPrefix(chainName, "private/")
            if isPersonal:
                chainName = f"u/{chainName}"
            else:
                chainName, _ = _dropPrefix(chainName, "shared/")
                chainName = instrument.makeCollectionName("runs", chainName)
            self.chainName = chainName
        self.runName = f"{self.chainName}/direct"


@dataclass
class CalibRepo:
    """Specification for a Gen2 calibration repository to convert.
    """

    path: Optional[str]
    """Absolute or relative (to the root repository) path to the Gen2
    repository (`str` or `None`).

    If `None`, no calibration datasets will be converted from Gen2, but
    curated calibrations may still be written.
    """

    curated: bool = True
    """If `True`, write curated calibrations into the associated
    ``CALIBRATION`` collection (`bool`).
    """

    labels: Tuple[str, ...] = ()
    """Extra strings to insert into collection names, including both the
    ``RUN`` collections that datasets are ingested directly into and the
    ``CALIBRATION`` collection that associates them with validity ranges.

    An empty tuple will directly populate the default calibration collection
    for this instrument with the converted datasets, and is incompatible with
    ``default=False``.  This is a good choice for test data repositories where
    only one ``CALIBRATION`` collection will ever exist.  In other cases, this
    should be a non-empty tuple, so the default calibration collection can
    actually be a ``CHAINED`` collection pointer that points to the current
    recommended ``CALIBRATION`` collection.
    """

    default: bool = True
    """If `True`, the created ``CALIBRATION`` collection should be the default
    for this instrument.

    This field may only be `True` for one converted calibration collection if
    more than one is passed to `ConvertRepoTask.run`.  It defaults to `True`
    because the vast majority of the time only one calibration collection is
    being converted.  If ``labels`` is not empty, ``default=True`` will cause
    a ``CHAINED`` collection that points to the converted ``CALIBRATION``
    collection to be defined.  If ``labels`` is empty, ``default`` *must* be
    `True` and no ``CHAINED`` collection pointer is necessary.
    """

    def __post_init__(self) -> None:
        if not self.labels and not self.default:
            raise ValueError("labels=() requires default=True")


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
        default={},
    )
    runsForced = DictField(
        "Like ``runs``, but is used even when the dataset is present in a "
        "non-root repository (i.e. rerun), overriding the non-root "
        "repository's main collection.",
        keytype=str,
        itemtype=str,
        default={
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
            "crosstalk": "CrosstalkCalib",
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
    datasetTemplateOverrides = DictField(
        "Overrides for Gen2 filename templates, keyed by dataset type. "
        "This can be used to support conversions of Gen2 repos whose mapper "
        "templates were modified in obs_* packages since the datasets were "
        "written.",
        keytype=str,
        itemtype=str,
        default={},
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
    doExpandDataIds = Field(
        "If True (default), expand data IDs to include extra metadata before "
        "ingesting them. "
        "This may be required in order to associate calibration datasets with "
        "validity ranges or populate file templates, so setting this to False "
        "is considered advanced usage (and it may not always work).  When it "
        "does, it can provide a considerable speedup.",
        dtype=bool,
        default=True,
    )
    doMakeUmbrellaCollection = Field(
        "If True (default), define an '<instrument>/defaults' CHAINED "
        "collection that includes everything found in the root repo as well "
        "as the default calibration collection.",
        dtype=bool,
        default=True,
    )
    extraUmbrellaChildren = ListField(
        "Additional child collections to include in the umbrella collection. "
        "Ignored if doMakeUmbrellaCollection=False.",
        dtype=str,
        default=[]
    )

    @property
    def transfer(self):
        return self.raws.transfer

    @transfer.setter
    def transfer(self, value):
        self.raws.transfer = value

    def setDefaults(self):
        self.transfer = None

    def validate(self):
        super().validate()
        if self.relatedOnly and not self.doExpandDataIds():
            raise ValueError("relatedOnly requires doExpandDataIds.")


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
    instrument : `lsst.obs.base.Instrument`
        The Gen3 instrument that should be used for this conversion.
    dry_run : `bool`, optional
        If `True` (`False` is default), make no changes to the Gen3 data
        repository while running as many steps as possible.  This option is
        best used with a read-only ``butler3`` argument to ensure unexpected
        edge cases respect this argument (and fail rather than write if they
        do not).
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

    def __init__(self, config=None, *, butler3: Butler3, instrument: Instrument, dry_run: bool = False,
                 **kwargs):
        config.validate()  # Not a CmdlineTask nor PipelineTask, so have to validate the config here.
        super().__init__(config, **kwargs)
        # Make self.butler3 one that doesn't have any collections associated
        # with it - those are needed by RawIngestTask and DefineVisitsTask, but
        # we don't want them messing with converted datasets, because those
        # have their own logic for figuring out which collections to write to.
        self.butler3 = Butler3(butler=butler3)
        self.registry = self.butler3.registry
        self.universe = self.registry.dimensions
        if self.isDatasetTypeIncluded("raw"):
            self.makeSubtask("raws", butler=butler3)
            self.makeSubtask("defineVisits", butler=butler3)
        else:
            self.raws = None
            self.defineVisits = None
        self.instrument = instrument
        self._configuredSkyMapsBySha1 = {}
        self._configuredSkyMapsByName = {}
        for name, config in self.config.skyMaps.items():
            instance = config.skyMap.apply()
            self._populateSkyMapDicts(name, instance)
        self._usedSkyPix = set()
        self.translatorFactory = self.instrument.makeDataIdTranslatorFactory()
        self.translatorFactory.log = self.log.getChild("translators")
        self.dry_run = dry_run

    def _reduce_kwargs(self):
        # Add extra parameters to pickle
        return dict(**super()._reduce_kwargs(), butler3=self.butler3, instrument=self.instrument)

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
                if not self.dry_run:
                    try:
                        # If the skymap isn't registerd, this will raise.
                        self.butler3.registry.expandDataId(skymap=struct.name)
                    except LookupError:
                        self.log.info("Registering skymap %s.", struct.name)
                        struct.instance.register(struct.name, self.butler3)
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
            calibs: Optional[List[CalibRepo]] = None,
            reruns: Optional[List[Rerun]] = None,
            visits: Optional[Iterable[int]] = None,
            pool: Optional[Pool] = None,
            processes: int = 1):
        """Convert a group of related data repositories.

        Parameters
        ----------
        root : `str`
            Complete path to the root Gen2 data repository.  This should be
            a data repository that includes a Gen2 registry and any raw files
            and/or reference catalogs.
        calibs : `list` of `CalibRepo`
            Specifications for Gen2 calibration repos to convert.  If `None`
            (default), curated calibrations only will be written to the default
            calibration collection for this instrument; set to ``()`` explictly
            to disable this.
        reruns : `list` of `Rerun`
            Specifications for rerun (processing output) repos to convert.  If
            `None` (default), no reruns are converted.
        visits : iterable of `int`, optional
            The integer IDs of visits to convert.  If not provided, all visits
            in the Gen2 root repository will be converted.
        pool : `multiprocessing.Pool`, optional
            If not `None`, a process pool with which to parallelize some
            operations.
        processes : `int`, optional
            The number of processes to use for conversion.
        """
        if pool is None and processes > 1:
            pool = Pool(processes)
        if calibs is None:
            calibs = [CalibRepo(path=None)]
        elif calibs and not self.config.doExpandDataIds:
            raise ValueError("Cannot convert calib repos with config.doExpandDataIds=False.")
        if visits is not None:
            subset = ConversionSubset(instrument=self.instrument.getName(), visits=frozenset(visits))
        else:
            if self.config.relatedOnly:
                self.log.warn("config.relatedOnly is True but all visits are being ingested; "
                              "no filtering will be done.")
            subset = None
        if (not self.config.doExpandDataIds
                and self.butler.datastore.needs_expanded_data_ids(self.config.transfer)):
            self.log.warn("config.doExpandDataIds=False but datastore reports that expanded data "
                          "IDs may be needed.",
                          self.config.transfer)

        # Check that at most one CalibRepo is marked as default, to fail before
        # we actually write anything.
        defaultCalibRepos = [c.path for c in calibs if c.default]
        if len(defaultCalibRepos) > 1:
            raise ValueError(f"Multiple calib repos marked as default: {defaultCalibRepos}.")

        # Make converters for all Gen2 repos.
        converters = []
        # Start with the root repo, which must always be given even if we are
        # not configured to convert anything from it.
        rootConverter = RootRepoConverter(task=self, root=root, subset=subset, instrument=self.instrument)
        converters.append(rootConverter)
        # Calibration repos are next.
        for spec in calibs:
            calibRoot = spec.path
            if calibRoot is not None:
                if not os.path.isabs(calibRoot):
                    calibRoot = os.path.join(rootConverter.root, calibRoot)
                converter = CalibRepoConverter(task=self, root=calibRoot,
                                               labels=spec.labels,
                                               instrument=self.instrument,
                                               mapper=rootConverter.mapper,
                                               subset=rootConverter.subset)
                converters.append(converter)
            # CalibRepo entries that don't have a path are just there for
            # curated calibs and maybe to set up a collection pointer; that's
            # handled further down (after we've done everything we can that
            # doesn't involve actually writing to the output Gen3 repo).
        # And now reruns.
        rerunConverters = {}
        for spec in reruns:
            runRoot = spec.path
            if not os.path.isabs(runRoot):
                runRoot = os.path.join(rootConverter.root, runRoot)
            spec.guessCollectionNames(self.instrument, rootConverter.root)
            converter = StandardRepoConverter(task=self, root=runRoot, run=spec.runName,
                                              instrument=self.instrument, subset=rootConverter.subset)
            converters.append(converter)
            rerunConverters[spec.runName] = converter

        # Walk Gen2 repos to find datasets to convert.
        for converter in converters:
            converter.prep()

        # Register the instrument if we're configured to do so.
        if self.config.doRegisterInstrument and not self.dry_run:
            self.instrument.register(self.registry)

        # Run raw ingest (does nothing if we weren't configured to convert the
        # 'raw' dataset type).
        rootConverter.runRawIngest(pool=pool)

        # Write curated calibrations to all calibration collections where they
        # were requested (which may be implicit, by passing calibs=None).  Also
        # set up a CHAINED collection that points to the default CALIBRATION
        # collection if one is needed.
        if not self.dry_run:
            for spec in calibs:
                if spec.curated:
                    self.instrument.writeCuratedCalibrations(self.butler3, labels=spec.labels)
                if spec.default and spec.labels:
                    # This is guaranteed to be True at most once in the loop by
                    # logic at the top of this method.
                    defaultCalibName = self.instrument.makeCalibrationCollectionName()
                    self.butler3.registry.registerCollection(defaultCalibName, CollectionType.CHAINED)
                    recommendedCalibName = self.instrument.makeCalibrationCollectionName(*spec.labels)
                    self.butler3.registry.registerCollection(recommendedCalibName, CollectionType.CALIBRATION)
                    self.butler3.registry.setCollectionChain(defaultCalibName, [recommendedCalibName])

        # Define visits (also does nothing if we weren't configurd to convert
        # the 'raw' dataset type).
        rootConverter.runDefineVisits(pool=pool)

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
        if self.config.doExpandDataIds:
            for converter in converters:
                converter.expandDataIds()

        if self.dry_run:
            return

        # Actually ingest datasets.
        for converter in converters:
            converter.ingest()

        # Perform any post-ingest processing.
        for converter in converters:
            converter.finish()

        # Make the umbrella collection, if desired.
        if self.config.doMakeUmbrellaCollection:
            umbrella = self.instrument.makeUmbrellaCollectionName()
            self.registry.registerCollection(umbrella, CollectionType.CHAINED)
            children = list(self.registry.getCollectionChain(umbrella))
            children.extend(rootConverter.getCollectionChain())
            children.append(self.instrument.makeCalibrationCollectionName())
            if BaseSkyMap.SKYMAP_RUN_COLLECTION_NAME not in children:
                # Ensure the umbrella collection includes the global skymap
                # collection, even if it's currently empty.
                self.registry.registerRun(BaseSkyMap.SKYMAP_RUN_COLLECTION_NAME)
                children.append(BaseSkyMap.SKYMAP_RUN_COLLECTION_NAME)
            children.extend(self.config.extraUmbrellaChildren)
            self.log.info("Defining %s from chain %s.", umbrella, children)
            self.registry.setCollectionChain(umbrella, children)

        # Add chained collections for reruns.
        for spec in reruns:
            if spec.chainName is not None:
                self.butler3.registry.registerCollection(spec.chainName, type=CollectionType.CHAINED)
                chain = [spec.runName]
                chain.extend(rerunConverters[spec.runName].getCollectionChain())
                for parent in spec.parents:
                    chain.append(parent)
                    parentConverter = rerunConverters.get(parent)
                    if parentConverter is not None:
                        chain.extend(parentConverter.getCollectionChain())
                chain.extend(rootConverter.getCollectionChain())
                if len(calibs) == 1:
                    # Exactly one calibration repo being converted, so it's
                    # safe-ish to assume that's the one the rerun used.
                    chain.append(self.instrument.makeCalibrationCollectionName(*calibs[0].labels))
                self.log.info("Defining %s from chain %s.", spec.chainName, chain)
                self.butler3.registry.setCollectionChain(spec.chainName, chain, flatten=True)
