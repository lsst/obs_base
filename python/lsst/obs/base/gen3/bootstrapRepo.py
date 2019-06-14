# This file is part of obs_base.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
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


__all__ = ("BootstrapRepoConfig", "BootstrapRepoTask", "BootstrapRepoInputs",
           "BootstrapRepoSkyMapConfig", "BootstrapRepoRefCatConfig")

import os.path
from dataclasses import dataclass
from typing import List
import glob

from lsst import sphgeom
from lsst.daf.butler import Butler, DatasetType
from lsst.daf.butler.instrument import Instrument
from lsst.pex.config import Config, Field, ConfigurableField, ConfigDictField, ConfigField
from lsst.pipe.base import Task
from lsst.obs.base.gen3 import RawIngestTask, makeTransferChoiceField
from lsst.skymap import skyMapRegistry
from lsst.meas.algorithms import DatasetConfig

from .repoConverter import RepoConverter
from .calibRepoConverter import CalibRepoConverter


class BootstrapRepoSkyMapConfig(Config):
    datasetTypeName = Field(("DatasetType used to write the SkyMap instance.  If None, the instance will "
                             "not be written, and only the Registry will be modified."),
                            dtype=str, default="deepCoadd_skyMap", optional=True)
    collection = Field(("Butler collection the SkyMap instance should be written to.  If None, the "
                        "collection used to initialize the butler will be used."),
                       dtype=str, default="skymaps", optional=True)
    skyMap = skyMapRegistry.makeField(
        doc="Type and parameters for the SkyMap itself.",
        default="dodeca",
    )


class BootstrapRepoRefCatConfig(Config):
    datasetTypeName = Field(("DatasetType used to write the catalog shards.."),
                            dtype=str, default="ref_cat")
    filterByRawRegions = Field(("If True, do not ingest shards that do not overlap visits.  "
                                "Does not guarantee that all ingested shards will overlap a visit."),
                               dtype=bool, default=True)
    collection = Field(("Butler collection the reference catalog should be written to.  If None, the "
                        "collection used to initialize the butler will be used.  May also be a string with "
                        "the format placeholder '{name}', which will be replaced with the reference "
                        "catalog name (i.e. the key of the configuration dictionary,"),
                       dtype=str, default="refcats/{name}", optional=True)
    transfer = makeTransferChoiceField(default="symlink")


class BootstrapRepoGenericIngestConfig(Config):
    collection = Field(("Butler collection that datasets should be ingested into.  "
                        "If None, the collection used to initialize the butler will be used."),
                       dtype=str, default=None, optional=True)
    transfer = makeTransferChoiceField(default="symlink")


class BootstrapRepoBrightObjectMasksConfig(BootstrapRepoGenericIngestConfig):
    skymap = Field("SkyMap dimension name used to define the tracts and patches for bright object masks.",
                   dtype=str, default=None, optional=False)
    filterByRawRegions = Field(("If True, do not ingest files that do not overlap visits.  "
                                "Does not guarantee that all ingested files will overlap a visit."),
                               dtype=bool, default=True)


class BootstrapRepoConfig(Config):
    raws = ConfigurableField(target=RawIngestTask,
                             doc=("Configuration for subtask responsible for ingesting raws and adding "
                                  "visit and exposure dimension entries."))
    skymaps = ConfigDictField(doc=("SkyMap definitions to register and ingest into the repo, keyed by "
                                   "skymap dimension name."),
                              keytype=str,
                              itemtype=BootstrapRepoSkyMapConfig,
                              default={})
    refCats = ConfigDictField(doc=("Reference catalogs to ingest into the repo, keyed by their subdirectory "
                                   "within the overall reference catalog root."),
                              keytype=str,
                              itemtype=BootstrapRepoRefCatConfig,
                              default={})
    brightObjectMasks = ConfigField(doc="Configuration for ingesting brightObjectMask files.",
                                    dtype=BootstrapRepoBrightObjectMasksConfig)
    calibrations = ConfigField(doc="Configuration for ingesting and creating master calibration products.",
                               dtype=BootstrapRepoGenericIngestConfig)

    def setDefaults(self):
        self.raws.transfer = "symlink"


@dataclass
class BootstrapRepoInputs:
    """Simple struct that aggregates all non-config inputs to
    `BootstrapRepoTask`.

    Generally, this stuct contains inputs that depend on the organization
    of the input files on a particular system, while the config includes
    everything else.  The exception is the ``instrument`` attribute, which
    cannot be included in the config because it's expected that driver code
    will actually use it (via
    `~lsst.daf.butler.instrument.Instrument.applyConfigOverrides`) to define
    the config.
    """

    instrument: Instrument
    """Instrument subclass instance for the raws and calibrations to be
    included in the initial repo.
    """

    raws: List[str]
    """List of filenames for raw files to ingest (complete paths).
    """

    refCatRoot: str
    """Root of the directory containing the reference catalogs, with immediate
    subdirectories that correspond to different reference catalogs.
    """

    brightObjectMaskRoot: str
    """Root of the Gen2 repository containing bright object masks.
    """

    calibRoot: str
    """Root of the Gen2 calibraion repository containing flats, biases,
    darks, and fringes.
    """


class BootstrapRepoTask(Task):
    """A Task that populates a Gen3 repo with the minimum content needed to
    run the DRP pipelines.

    BootstrapRepoTask currently relies on Gen2 data repository information
    for both bright object masks and master calibrations, but nothing else;
    unlike dedicated Gen2->Gen3 conversion code, it will be updated in the
    future as more pure-Gen3 approaches become available.

    Like other Gen3 Tasks that are not PipelineTasks, BootstrapRepoTask does
    not yet have a dedicated, general-purpose command-line driver.  At least
    for now, it is instead expected that custom driver scripts will be written
    for different contexts and predefined datasets.

    Parameters
    ----------
    config : `BootstrapRepoConfig`
        Configuration for the task.
    butler : `lsst.daf.butler.Butler`
        Gen3 Butler defining the repository to populate.  New butlers with
        different output collections will be created as necessary from this
        butler to match the output collections defined in the configuration.
    kwds
        Additional keyword arguments are forwarded to the
        `lsst.pipe.base.Task` constructor.
    """

    ConfigClass = BootstrapRepoConfig

    _DefaultName = "bootstrapRepo"

    def __init__(self, config=None, *, butler, **kwds):
        super().__init__(config, **kwds)
        self.butler = butler
        self.makeSubtask("raws", butler=self.butler)
        self.skyMaps = {}

    def getButler(self, collection=None):
        """Create a new butler that writes into the given collection.

        Parameters
        ----------
        collection : `str`, optional
            The new output collection.  If `None`, ``self.butler`` is returned
            directly.

        Returns
        -------
        butler : `lsst.daf.butler.Butler`
            Butler instance pointing at the same repository as
            ``self.butler``, but possibly a different collection.
        """
        if collection is not None:
            return Butler(butler=self.butler, run=collection)
        return self.butler

    def run(self, inputs):
        """Run all steps involved in populating the new repository.

        Parameters
        ----------
        inputs : `BootstrapRepoInputs`
            Filenames and paths for the data to be ingested.
        """
        self.bootstrapInstrument(inputs.instrument)
        self.bootstrapCalibrations(inputs.instrument, inputs.calibRoot)
        self.bootstrapRaws(inputs.raws)
        self.bootstrapRefCats(inputs.refCatRoot)
        self.bootstrapSkyMaps()
        self.bootstrapBrightObjectMasks(inputs.instrument, inputs.brightObjectMaskRoot)

    def bootstrapInstrument(self, instrument):
        """Add an instrument, associated metadata, and human-curated
        calibrations to the repository.

        Parameters
        ----------
        instrument : `lsst.daf.butler.instrument.Instrument`
            Instrument class that defines detectors, physical filters, and
            curated calibrations to ingest.
        """
        self.log.info("Registering instrument '%s' and adding curated calibrations.", instrument.getName())
        with self.butler.transaction():
            instrument.register(self.butler.registry)
            instrument.writeCuratedCalibrations(self.getButler(self.config.calibrations.collection))

    def bootstrapSkyMaps(self):
        """Add configured SkyMaps to the repository.

        This both registers skymap dimension entries (the skymap, tract, and
        patch tables, and their associated join tables) and adds a
        ``<something>Coadd_skyMap`` dataset.
        """
        for name, config in self.config.skymaps.items():
            self.log.info("Registering skymap '%s'.", name)
            with self.butler.transaction():
                skyMap = config.skyMap.apply()
                skyMap.register(name, self.butler.registry)
                if config.datasetTypeName is not None:
                    datasetType = DatasetType(config.datasetTypeName, dimensions=["skymap"],
                                              storageClass="SkyMap",
                                              universe=self.butler.registry.dimensions)
                    self.butler.registry.registerDatasetType(datasetType)
                    self.getButler(config.collection).put(skyMap, datasetType, skymap=name)
            self.skyMaps[name] = skyMap

    def bootstrapRaws(self, files):
        """Ingest raw images.

        This step must be run after `bootstrapInstrument`, but may be run
        multiple times with different arguments (which may be overlapping if
        the nested `RawIngestTask` is configured to ignore duplicates).

        Parameters
        ----------
        files : sequence of `str`
            The complete path names of the files to be ingested.
        """
        self.log.info("Ingesting raw images.")
        return self.raws.run(files)  # transaction handled internally, according to config.

    def computeRawSkyPixels(self):
        """Compute and return the skypix dimension entries that overlap
        already-ingested visits.
        """
        # TODO: provide a non-SQL way to efficiently perform this query?
        return list(
            row["skypix"] for row in self.butler.registry.query(
                "SELECT DISTINCT skypix FROM visit_skypix_join"
            )
        )

    def bootstrapRefCats(self, root):
        """Ingest reference catalogs.

        This step must be run after `bootstrapRaws` if the
        ``filterByRawRegions`` config option is `True` for any reference
        catalog.

        Parameters
        ----------
        root : `str`
            Root of the directory containing the reference catalogs, with
            immediate subdirectories that correspond to different reference
            catalogs.
        """
        if not self.config.refCats:
            return
        if any(config.filterByRawRegions for config in self.config.refCats.values()):
            rawSkyPixels = self.computeRawSkyPixels()
        datasetType = DatasetType("ref_cat", dimensions=["skypix"], storageClass="SimpleCatalog",
                                  universe=self.butler.registry.dimensions)
        self.butler.registry.registerDatasetType(datasetType)
        for name, config in self.config.refCats.items():
            self.log.info("Ingesting reference catalog '%s'.", name)
            with self.butler.transaction():
                onDiskConfig = DatasetConfig()
                onDiskConfig.load(os.path.join(root, name, "config.py"))
                if onDiskConfig.indexer.name != "HTM":
                    raise ValueError(f"Reference catalog '{name}' uses unsupported "
                                     f"pixelization '{onDiskConfig.indexer.name}'.")
                if not isinstance(self.butler.registry.pixelization, sphgeom.HtmPixelization):
                    raise ValueError(f"Registry uses unsupported pixelization class "
                                     f"{self.butler.registry.pixelization.__class__}.")
                if onDiskConfig.indexer["HTM"].depth != self.butler.registry.pixelization.getLevel():
                    raise ValueError(f"Registry HTM level {self.butler.registry.pixelization.getLevel()} "
                                     f"does not match reference catalog level {onDiskConfig.indexer.depth}.")
                butler = self.getButler(config.collection.format(name))
                if config.filterByRawRegions:
                    missing = []
                    for index in rawSkyPixels:
                        path = os.path.join(root, name, f"{index}.fits")
                        if os.path.exists(path):
                            butler.ingest(path, datasetType, transfer=config.transfer, skypix=index)
                        else:
                            missing.append(index)
                    if missing:
                        self.log.warn("Some overlapping reference catalog shards missing: %s", missing)
                else:
                    for path in glob.glob(os.path.join(root, name, "*.fits")):
                        if path.endswith("master_schema.fits"):
                            continue
                        _, filename = os.path.split(path)
                        basename, _ = os.path.splitext(filename)
                        try:
                            index = int(basename)
                        except ValueError:
                            self.log.warn("Unrecognized file in reference catalog root: '%s'.", path)
                            continue
                        butler.ingest(path, datasetType, transfer=config.transfer, skypix=index)

    def computeRawTracts(self, skymap):
        """Compute and return the tract dimension entries that overlap
        already-ingested visits.
        """
        # TODO: provide a non-SQL way to efficiently perform this query?
        return list(
            row["tract"] for row in self.butler.registry.query(
                "SELECT DISTINCT tract FROM visit_tract_join WHERE skymap=:skymap",
                skymap=skymap
            )
        )

    def bootstrapBrightObjectMasks(self, instrument, root):
        """Ingest bright object masks from a Gen2 data repository.

        This step must be run after `bootstrapRaws` if the
        ``filterByRawRegions`` config option is `True` for any reference
        catalog, and must always be run after `bootstrapSkyMaps`.

        Parameters
        ----------
        root : `str`
            Root of the Gen2 repository containing bright object masks.
        instrument : `lsst.daf.butler.instrument.Instrument`
            Instrument subclass instance; used to relate Gen2 filter
            strings to Gen3 physical_filters and abstract_filters.
        """
        self.log.info("Ingesting bright object masks.")
        butler = self.getButler(self.config.brightObjectMasks.collection)
        baseDataId = {
            "skymap": self.config.brightObjectMasks.skymap,
            "instrument": instrument.getName()
        }
        converter = RepoConverter(root, universe=butler.registry.dimensions, baseDataId=baseDataId,
                                  skyMap=self.skyMaps[self.config.brightObjectMasks.skymap])
        converter.addDatasetType("brightObjectMask", "ObjectMaskCatalog")
        if self.config.brightObjectMasks.filterByRawRegions:
            for tract in self.computeRawTracts(self.config.brightObjectMasks.skymap):
                with self.butler.transaction():
                    converter.convertRepo(butler, directory=f"{root}/deepCoadd/BrightObjectMasks/{tract:d}",
                                          transfer=self.config.brightObjectMasks.transfer)
        else:
            with self.butler.transaction():
                converter.convertRepo(butler, transfer=self.config.brightObjectMasks.transfer)

    def bootstrapCalibrations(self, instrument, root):
        """Ingest master calibrations from a Gen2 calibration data repository.

        At present, all master calibrations in the Gen2 repostory are
        transferred, even those unrelated to the ingested raws.

        This step must be run after `bootstrapInstrument`.

        Parameters
        ----------
        instrument : `lsst.daf.butler.instrument.Instrument`
            Instrument subclass instance for the raws and calibrations to be
            included in the initial repo.
        root : `str`
            Root of the Gen2 calibration data repository.
        """
        self.log.info("Ingesting calibrations.")
        baseDataId = {"instrument": instrument.getName()}
        butler = self.getButler(self.config.calibrations.collection)
        converter = CalibRepoConverter(root, universe=butler.registry.dimensions, baseDataId=baseDataId)
        converter.addDatasetType("flat", "MaskedImageF")
        converter.addDatasetType("bias", "ImageF")
        converter.addDatasetType("dark", "ImageF")
        converter.addDatasetType("sky", "ExposureF")
        converter.addDatasetType("fringe", "ExposureF")
        # TODO, DM-16805: No StorageClass/Formatter for yBackground in Gen3.
        with self.butler.transaction():
            converter.convertRepo(butler, transfer=self.config.brightObjectMasks.transfer)
