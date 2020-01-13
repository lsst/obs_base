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

__all__ = ["ConvertRepoConfig", "ConvertRepoTask", "ConvertRepoSkyMapConfig"]

import os
import fnmatch
from dataclasses import dataclass
from typing import Iterable, Optional, List, Dict

from lsst.utils import doImport
from lsst.daf.butler import (
    Butler as Butler3,
    SkyPixDimension
)
from lsst.pex.config import Config, ConfigurableField, ConfigDictField, DictField, ListField, Field
from lsst.pipe.base import Task
from lsst.skymap import skyMapRegistry, BaseSkyMap

from ..ingest import RawIngestTask
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
        "visit and exposure dimension entries.",
        target=RawIngestTask,
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
    collections = DictField(
        "Special collections (values) for certain dataset types (keys).  "
        "These are used in addition to rerun collections for datasets in "
        "reruns.  The 'raw' dataset must have an entry here if it is to be "
        "converted.",
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
            "BaseSkyMap": "SkyMap",
            "BaseCatalog": "Catalog",
            "BackgroundList": "Background",
            "raw": "Exposure",
            "MultilevelParquetTable": "DataFrame",
        }
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

    @property
    def transfer(self):
        return self.raws.transfer

    @transfer.setter
    def transfer(self, value):
        self.raws.transfer = value

    @property
    def instrument(self):
        return self.raws.instrument

    @instrument.setter
    def instrument(self, value):
        self.raws.instrument = value

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
        Gen3 Butler instance that represents the data repository datasets will
        be ingested into.  The collection and/or run associated with this
        Butler will be ignored in favor of collections/runs passed via config
        or to `run`.
    kwds
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

    def __init__(self, config=None, *, butler3: Butler3, **kwds):
        super().__init__(config, **kwds)
        self.butler3 = butler3
        self.registry = self.butler3.registry
        self.universe = self.registry.dimensions
        if self.isDatasetTypeIncluded("raw"):
            self.makeSubtask("raws", butler=butler3)
            self.instrument = self.raws.instrument
        else:
            self.raws = None
            self.instrument = doImport(self.config.instrument)()
        self._configuredSkyMapsBySha1 = {}
        self._configuredSkyMapsByName = {}
        for name, config in self.config.skyMaps.items():
            instance = config.skyMap.apply()
            struct = ConfiguredSkyMap(name=name, sha1=instance.getSha1(), instance=instance)
            self._configuredSkyMapsBySha1[struct.sha1] = struct
            self._configuredSkyMapsByName[struct.name] = struct
        self._usedSkyPix = set()

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
                for pattern in self.config.datasetIncludePatterns) and
            not any(fnmatch.fnmatchcase(datasetTypeName, pattern)
                    for pattern in self.config.datasetIgnorePatterns)
        )

    def useSkyMap(self, skyMap: BaseSkyMap) -> str:
        """Indicate that a repository uses the given SkyMap.

        This method is intended to be called primarily by the
        `RepoConverter` instances used interally by the task.

        Parameters
        ----------
        skyMap : `lsst.skymap.BaseSkyMap`
            SkyMap instance being used, typically retrieved from a Gen2
            data repository.

        Returns
        -------
        name : `str`
            The name of the skymap in Gen3 data IDs.
        """
        sha1 = skyMap.getSha1()
        try:
            struct = self._configuredSkyMapsBySha1[sha1]
        except KeyError as err:
            raise LookupError(f"SkyMap with sha1={sha1} not included in configuration.") from err
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

    def run(self, root: str, collections: List[str], *,
            calibs: Dict[str, List[str]] = None,
            reruns: Dict[str, List[str]] = None,
            visits: Optional[Iterable[int]] = None):
        """Convert a group of related data repositories.

        Parameters
        ----------
        root : `str`
            Complete path to the root Gen2 data repository.  This should be
            a data repository that includes a Gen2 registry and any raw files
            and/or reference catalogs.
        collections : `list` of `str`
            Gen3 collections that datasets from the root repository should be
            associated with.  This should include any rerun collection that
            these datasets should also be considered to be part of; because of
            structural difference between Gen2 parent/child relationships and
            Gen3 collections, these cannot be reliably inferred.
        calibs : `dict`
            Dictionary mapping calibration repository path to the collections
            that the repository's datasets should be associated with.  The path
            may be relative to ``root`` or absolute.  Collections should
            include child repository collections as appropriate (see
            documentation for ``collections``).
        reruns : `dict`
            Dictionary mapping rerun repository path to the collections that
            the repository's datasets should be associated with.  The path may
            be relative to ``root`` or absolute.  Collections should include
            child repository collections as appropriate (see documentation for
            ``collections``).
        visits : iterable of `int`, optional
            The integer IDs of visits to convert.  If not provided, all visits
            in the Gen2 root repository will be converted.
        """

        if calibs is None:
            calibs = {}
        if reruns is None:
            reruns = {}
        if visits is not None:
            subset = ConversionSubset(instrument=self.instrument.getName(), visits=frozenset(visits))
        else:
            if self.config.relatedOnly:
                self.log.warn("config.relatedOnly is True but all visits are being ingested; "
                              "no filtering will be done.")
            subset = None

        # We can't wrap database writes sanely in transactions (yet) because we
        # keep initializing new Butler instances just so we can write into new
        # runs/collections, and transactions are managed at the Butler level.
        # DM-21246 should let us fix this, assuming we actually want to keep
        # the transaction open that long.
        if self.config.doRegisterInstrument:
            self.instrument.register(self.registry)

        # Make and prep converters for all Gen2 repos.  This should not modify
        # the Registry database or filesystem at all, though it may query it.
        # The prep() calls here will be some of the slowest ones, because
        # that's when we walk the filesystem.
        converters = []
        rootConverter = RootRepoConverter(task=self, root=root, collections=collections, subset=subset)
        rootConverter.prep()
        converters.append(rootConverter)

        for root, collections in calibs.items():
            if not os.path.isabs(root):
                root = os.path.join(rootConverter.root, root)
            converter = CalibRepoConverter(task=self, root=root, collections=collections,
                                           mapper=rootConverter.mapper,
                                           subset=rootConverter.subset)
            converter.prep()
            converters.append(converter)

        for root, collections in reruns.items():
            if not os.path.isabs(root):
                root = os.path.join(rootConverter.root, root)
            converter = StandardRepoConverter(task=self, root=root, collections=collections,
                                              subset=rootConverter.subset)
            converter.prep()
            converters.append(converter)

        # Actual database writes start here.  We can't wrap these sanely in
        # transactions (yet) because we keep initializing new Butler instances
        # just so we can write into new runs/collections, and transactions
        # are managed at the Butler level (DM-21246 should let us fix this).

        # Insert dimensions needed by any converters.  These are only the
        # dimensions that a converter expects to be uniquely derived from the
        # Gen2 repository it is reponsible for - e.g. visits, exposures, and
        # calibration_labels.
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
