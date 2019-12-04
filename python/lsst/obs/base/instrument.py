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

__all__ = ("Instrument", "makeExposureRecordFromObsInfo", "makeVisitRecordFromObsInfo",
           "addUnboundedCalibrationLabel")

import os.path
from datetime import datetime
from abc import ABCMeta, abstractmethod

import lsst.log

from lsst.daf.butler import DatasetType, DataCoordinate
from lsst.pipe.tasks import read_curated_calibs


class Instrument(metaclass=ABCMeta):
    """Base class for instrument-specific logic for the Gen3 Butler.

    Concrete instrument subclasses should be directly constructable with no
    arguments.
    """

    filterDefinitions = None
    """`lsst.obs.base.filters.FilterDefinitions` defining the filters used by
    this instrument.
    """

    configPaths = []
    """Paths to config files to read for specific Tasks.

    The paths in this list should contain files of the form `task.py`, for
    each of the Tasks that requires special configuration.
    """

    dataPath = None
    """Path to the ``obs_*_data/instrument`` directory containing human-curated
    standardized text calibration products, e.g. ``defects/``.

    For example, for HyperSuprimeCam, this would be ``$OBS_SUBARU_DATA/hsc``.
    """

    @property
    @abstractmethod
    def filterDefinitions(self):
        """`~lsst.obs.base.FilterDefinitionCollection`, defining the filters
        for this instrument.
        """
        return None

    def __init__(self, *args, **kwargs):
        self.log = lsst.log.Log()
        self.filterDefinitions.defineFilters()

    @classmethod
    @abstractmethod
    def getName(cls):
        raise NotImplementedError()

    @abstractmethod
    def getCamera(self):
        """Retrieve the cameraGeom representation of this instrument.

        This is a temporary API that should go away once obs_ packages have
        a standardized approach to writing versioned cameras to a Gen3 repo.
        It is more future proof to use ``butler.get('camera')`` to get the
        camera geometry.
        """
        raise NotImplementedError()

    @abstractmethod
    def register(self, registry):
        """Insert instrument, physical_filter, and detector entries into a
        `Registry`.
        """
        raise NotImplementedError()

    def _registerFilters(self, registry):
        """Register the physical and abstract filter Dimension relationships.
        This should be called in the ``register`` implementation.

        Parameters
        ----------
        registry : `lsst.daf.butler.core.Registry`
            The registry to add dimensions to.
        """
        self.log.info("Registering physical_filter dimensions...")
        registry.insertDimensionData(
            "physical_filter",
            *[
                {
                    "instrument": self.getName(),
                    "name": filter.physical_filter,
                    "abstract_filter": filter.abstract_filter,
                }
                for filter in self.filterDefinitions
            ]
        )

    @abstractmethod
    def getRawFormatter(self, dataId):
        """Return the Formatter class that should be used to read a particular
        raw file.

        Parameters
        ----------
        dataId : `DataCoordinate`
            Dimension-based ID for the raw file or files being ingested.

        Returns
        -------
        formatter : `Formatter` class
            Class to be used that reads the file into an
            `lsst.afw.image.Exposure` instance.
        """
        raise NotImplementedError()

    def _getBrighterFatterKernel(self):
        """Return the brighter-fatter kernel as a `numpy.ndarray`, or `None`
        if your instrument does not have brighter-fatter data.
        """
        return None

    def _writeBrighterFatterKernel(self, butler, calibrationLabelDataId):
        """Write a brighter-fatter kernel to the butler with the specified
        calibrationLabelDataId.

        Specialize `_getBrighterFatterKernel` to get the necessary data.
        """
        bfKernel = self._getBrighterFatterKernel()
        if bfKernel is None:
            self.log.info("No brighter-fatter kernels to write.")
            return

        self.log.info("Writing brighter-fatter kernels...")
        datasetType = DatasetType("bfKernel", ("instrument", "calibration_label"), "NumpyArray",
                                  universe=butler.registry.dimensions)
        butler.registry.registerDatasetType(datasetType)
        butler.put(bfKernel, datasetType, calibrationLabelDataId)

    def _writeCamera(self, butler, calibrationLabelDataId):
        """Write an `lsst.afw.cameraGeom.Camera` to the butler with the
        specified calibrationLabelDataId.
        """
        self.log.info("Writing Camera geometry...")
        datasetType = DatasetType("camera", ("instrument", "calibration_label"), "Camera",
                                  universe=butler.registry.dimensions)
        butler.registry.registerDatasetType(datasetType)
        camera = self.getCamera()
        butler.put(camera, datasetType, calibrationLabelDataId)

    def _getDefects(self):
        """Return defects loaded from the obs data directory.
        """
        defectsPath = os.path.join(self.dataPath, "defects")
        if not os.path.exists(defectsPath):
            self.log.info("No defects to read.")
            return None
        camera = self.getCamera()
        self.log.info("Reading defects...")
        data_by_chip, calib_type = read_curated_calibs.read_all(defectsPath, camera)
        if calib_type != "defects":
            raise TypeError(f"Loaded {calib_type} instead of `defects` from: {defectsPath}")
        return data_by_chip

    def _writeDefects(self, butler):
        """Write a collection of `lsst.meas.algorithms.Defects` to the butler.
        """
        defectsDict = self._getDefects()
        if defectsDict is None:
            self.log.info("No defects to write...")
            return

        self.log.info("Writing defects...")
        datasetType = DatasetType("defects", ("instrument", "detector", "calibration_label"), "DefectsList",
                                  universe=butler.registry.dimensions)
        butler.registry.registerDatasetType(datasetType)

        dimensionRecords = []
        datasetRecords = []
        # First loop just gathers up the things we want to insert, so we
        # can do some bulk inserts and minimize the time spent in transaction.
        for det in defectsDict:
            times = sorted([k for k in defectsDict[det]])
            defects = [defectsDict[det][time] for time in times]
            # Use an "infinite" end time for the "oldest" defect
            times.append(datetime.max)
            for defect, beginTime, endTime in zip(defects, times[:-1], times[1:]):
                md = defect.getMetadata()
                detectorId = md['DETECTOR']
                calibrationLabel = f"defect/{md['CALIBDATE']}/{detectorId}"
                dataId = DataCoordinate.standardize(
                    universe=butler.registry.dimensions,
                    instrument=self.getName(),
                    calibration_label=calibrationLabel,
                    detector=detectorId,
                )
                datasetRecords.append((defect, dataId))
                dimensionRecords.append({
                    "instrument": self.getName(),
                    "name": calibrationLabel,
                    "datetime_begin": beginTime,
                    "datetime_end": endTime,
                })
        # Second loop actually does the inserts and filesystem writes.
        with butler.transaction():
            butler.registry.insertDimensionData("calibration_label", *dimensionRecords)
            # TODO: vectorize these puts, once butler APIs for that become
            # available.
            for defect, dataId in datasetRecords:
                butler.put(defect, datasetType, dataId)

    def _getDetectorTransmission(self, unboundedDataId):
        """Return a dictionary of detector tramission curves, keyed on their
        DataCoordinates.
        """
        return None

    def _getOpticsTransmission(self, unboundedDataId):
        """Return a dictionary of optics tramission curves, keyed on their
        DataCoordinates.
        """
        return None

    def _getFilterTransmission(self, unboundedDataId):
        """Return a dictionary of filter tramission curves, keyed on their
        DataCoordinates.
        """
        return None

    def _getAtmosphereTransmission(self, unboundedDataId):
        """Return a dictionary of atmosphere tramission curves, keyed on their
        DataCoordinates.
        """
        return None

    def _writeTransmission(self, butler, unboundedDataId):
        """Write a collection of `lsst.afw.image.TransmissionCurves` to the
        butler, loaded via the varoius ``_get*Tranmission`` methods, if they
        have been implemented.
        """
        detectorTransmission = self._getDetectorTransmission(unboundedDataId)

        if detectorTransmission is not None:
            self.log.info("Writing detector transmission curves...")
            # NOTE: we want to rename "transmission_sensor" to "transmission_detector" once we retire gen2.
            datasetType = DatasetType("transmission_sensor",
                                      ("instrument", "detector", "calibration_label"),
                                      "TransmissionCurve",
                                      universe=butler.registry.dimensions)
            self._putTransmissionCurves(butler, detectorTransmission, datasetType)
        else:
            self.log.info("No detector transmission curves to write.")

        opticsTransmission = self._getOpticsTransmission(unboundedDataId)
        if opticsTransmission is not None:
            self.log.info("Writing optics transmission curves...")
            datasetType = DatasetType("transmission_optics",
                                      ("instrument", "calibration_label"),
                                      "TransmissionCurve",
                                      universe=butler.registry.dimensions)
            print(opticsTransmission)
            self._putTransmissionCurves(butler, opticsTransmission, datasetType)
        else:
            self.log.info("No optics transmission curves to write.")

        filterTransmission = self._getFilterTransmission(unboundedDataId)
        if filterTransmission is not None:
            self.log.info("Writing filter transmission curves...")
            datasetType = DatasetType("transmission_filter",
                                      ("instrument", "physical_filter", "calibration_label"),
                                      "TransmissionCurve",
                                      universe=butler.registry.dimensions)
            self._putTransmissionCurves(butler, filterTransmission, datasetType)
        else:
            self.log.info("No filter transmission curves to write.")

        atmosphereTransmission = self._getAtmosphereTransmission(unboundedDataId)
        if atmosphereTransmission is not None:
            self.log.info("Writing atmosphere transmission curves...")
            datasetType = DatasetType("transmission_atmosphere", ("instrument",),
                                      "TransmissionCurve",
                                      universe=butler.registry.dimensions)
            self._putTransmissionCurves(butler, atmosphereTransmission, datasetType)
        else:
            self.log.info("No atmosphere transmission curves to write.")

    def _putTransmissionCurves(self, butler, transmissionCurves, datasetType):
        """Put a dictionary of transmissionCurves into the Butler.

        Parameters
        ----------
        butler : `lsst.daf.butler.Butler`
            The butler to write the data to.
        transmissionCurves : `dict` [`lsst.daf.butler.DataCoordinate`, `lsst.afw.image.TransmissionCurve`]
            The transmission curves to be written, with their corresponding
            dataIds.
        datasetType : `lsst.daf.butler.DatasetType`
            The category of dataset represented by ``transmissionCurves`` to
            be written to the butler.
        """
        butler.registry.registerDatasetType(datasetType)
        for dataId, transmissionCurve in transmissionCurves.items():
            butler.put(transmissionCurve, datasetType, dataId)

    def writeInstrumentSignatureData(self, butler):
        """Write Instrument Signature Datasets to the given Butler with
        the appropriate validity ranges.
        """
        unboundedDataId = addUnboundedCalibrationLabel(butler.registry, self.getName())

        self._writeCamera(butler, unboundedDataId)
        self._writeBrighterFatterKernel(butler, unboundedDataId)
        self._writeDefects(butler)
        self._writeTransmission(butler, unboundedDataId)

    def applyConfigOverrides(self, name, config):
        """Apply instrument-specific overrides for a task config.

        Parameters
        ----------
        name : `str`
            Name of the object being configured; typically the _DefaultName
            of a Task.
        config : `lsst.pex.config.Config`
            Config instance to which overrides should be applied.
        """
        for root in self.configPaths:
            path = os.path.join(root, f"{name}.py")
            if os.path.exists(path):
                config.load(path)


def makeExposureRecordFromObsInfo(obsInfo, universe):
    """Construct an exposure DimensionRecord from
    `astro_metadata_translator.ObservationInfo`.

    Parameters
    ----------
    obsInfo : `astro_metadata_translator.ObservationInfo`
        A `~astro_metadata_translator.ObservationInfo` object corresponding to
        the exposure.
    universe : `DimensionUniverse`
        Set of all known dimensions.

    Returns
    -------
    record : `DimensionRecord`
        A record containing exposure metadata, suitable for insertion into
        a `Registry`.
    """
    dimension = universe["exposure"]
    return dimension.RecordClass.fromDict({
        "instrument": obsInfo.instrument,
        "id": obsInfo.exposure_id,
        "name": obsInfo.observation_id,
        "datetime_begin": obsInfo.datetime_begin.to_datetime(),
        "datetime_end": obsInfo.datetime_end.to_datetime(),
        "exposure_time": obsInfo.exposure_time.to_value("s"),
        "dark_time": obsInfo.dark_time.to_value("s"),
        "observation_type": obsInfo.observation_type,
        "physical_filter": obsInfo.physical_filter,
        "visit": obsInfo.visit_id,
    })


def makeVisitRecordFromObsInfo(obsInfo, universe, *, region=None):
    """Construct a visit `DimensionRecord` from
    `astro_metadata_translator.ObservationInfo`.

    Parameters
    ----------
    obsInfo : `astro_metadata_translator.ObservationInfo`
        A `~astro_metadata_translator.ObservationInfo` object corresponding to
        the exposure.
    universe : `DimensionUniverse`
        Set of all known dimensions.
    region : `lsst.sphgeom.Region`, optional
        Spatial region for the visit.

    Returns
    -------
    record : `DimensionRecord`
        A record containing visit metadata, suitable for insertion into a
        `Registry`.
    """
    dimension = universe["visit"]
    return dimension.RecordClass.fromDict({
        "instrument": obsInfo.instrument,
        "id": obsInfo.visit_id,
        "name": obsInfo.observation_id,
        "datetime_begin": obsInfo.datetime_begin.to_datetime(),
        "datetime_end": obsInfo.datetime_end.to_datetime(),
        "exposure_time": obsInfo.exposure_time.to_value("s"),
        "physical_filter": obsInfo.physical_filter,
        "region": region,
    })


def addUnboundedCalibrationLabel(registry, instrumentName):
    """Add a special 'unbounded' calibration_label dimension entry for the
    given camera that is valid for any exposure.

    If such an entry already exists, this function just returns a `DataId`
    for the existing entry.

    Parameters
    ----------
    registry : `Registry`
        Registry object in which to insert the dimension entry.
    instrumentName : `str`
        Name of the instrument this calibration label is associated with.

    Returns
    -------
    dataId : `DataId`
        New or existing data ID for the unbounded calibration.
    """
    d = dict(instrument=instrumentName, calibration_label="unbounded")
    try:
        return registry.expandDataId(d)
    except LookupError:
        pass
    entry = d.copy()
    entry["datetime_begin"] = datetime.min
    entry["datetime_end"] = datetime.max
    registry.insertDimensionData("calibration_label", entry)
    return registry.expandDataId(d)
