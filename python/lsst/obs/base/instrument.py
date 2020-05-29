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

__all__ = ("Instrument", "makeExposureRecordFromObsInfo", "addUnboundedCalibrationLabel", "loadCamera")

import os.path
from abc import ABCMeta, abstractmethod
from typing import Any, Tuple, TYPE_CHECKING
import astropy.time

from lsst.afw.cameraGeom import Camera
from lsst.daf.butler import Butler, DataId, TIMESPAN_MIN, TIMESPAN_MAX, DatasetType, DataCoordinate
from lsst.utils import getPackageDir, doImport

if TYPE_CHECKING:
    from .gen2to3 import TranslatorFactory
    from lsst.daf.butler import Registry

# To be a standard text curated calibration means that we use a
# standard definition for the corresponding DatasetType.
StandardCuratedCalibrationDatasetTypes = {
    "defects": {"dimensions": ("instrument", "detector", "calibration_label"),
                "storageClass": "Defects"},
    "qe_curve": {"dimensions": ("instrument", "detector", "calibration_label"),
                 "storageClass": "QECurve"},
}


class Instrument(metaclass=ABCMeta):
    """Base class for instrument-specific logic for the Gen3 Butler.

    Concrete instrument subclasses should be directly constructable with no
    arguments.
    """

    configPaths = ()
    """Paths to config files to read for specific Tasks.

    The paths in this list should contain files of the form `task.py`, for
    each of the Tasks that requires special configuration.
    """

    policyName = None
    """Instrument specific name to use when locating a policy or configuration
    file in the file system."""

    obsDataPackage = None
    """Name of the package containing the text curated calibration files.
    Usually a obs _data package.  If `None` no curated calibration files
    will be read. (`str`)"""

    standardCuratedDatasetTypes = tuple(StandardCuratedCalibrationDatasetTypes)
    """The dataset types expected to be obtained from the obsDataPackage.
    These dataset types are all required to have standard definitions and
    must be known to the base class.  Clearing this list will prevent
    any of these calibrations from being stored. If a dataset type is not
    known to a specific instrument it can still be included in this list
    since the data package is the source of truth.
    """

    @property
    @abstractmethod
    def filterDefinitions(self):
        """`~lsst.obs.base.FilterDefinitionCollection`, defining the filters
        for this instrument.
        """
        return None

    def __init__(self, *args, **kwargs):
        self.filterDefinitions.reset()
        self.filterDefinitions.defineFilters()
        self._obsDataPackageDir = None

    @classmethod
    @abstractmethod
    def getName(cls):
        """Return the short (dimension) name for this instrument.

        This is not (in general) the same as the class name - it's what is used
        as the value of the "instrument" field in data IDs, and is usually an
        abbreviation of the full name.
        """
        raise NotImplementedError()

    @abstractmethod
    def getCamera(self):
        """Retrieve the cameraGeom representation of this instrument.

        This is a temporary API that should go away once obs_ packages have
        a standardized approach to writing versioned cameras to a Gen3 repo.
        """
        raise NotImplementedError()

    @abstractmethod
    def register(self, registry):
        """Insert instrument, physical_filter, and detector entries into a
        `Registry`.
        """
        raise NotImplementedError()

    @property
    def obsDataPackageDir(self):
        """The root of the obs package that provides specializations for
        this instrument (`str`).
        """
        if self.obsDataPackage is None:
            return None
        if self._obsDataPackageDir is None:
            # Defer any problems with locating the package until
            # we need to find it.
            self._obsDataPackageDir = getPackageDir(self.obsDataPackage)
        return self._obsDataPackageDir

    @staticmethod
    def fromName(name: str, registry: Registry) -> Instrument:
        """Given an instrument name and a butler, retrieve a corresponding
        instantiated instrument object.

        Parameters
        ----------
        name : `str`
            Name of the instrument (must match the return value of `getName`).
        registry : `lsst.daf.butler.Registry`
            Butler registry to query to find the information.

        Returns
        -------
        instrument : `Instrument`
            An instance of the relevant `Instrument`.

        Notes
        -----
        The instrument must be registered in the corresponding butler.

        Raises
        ------
        LookupError
            Raised if the instrument is not known to the supplied registry.
        ModuleNotFoundError
            Raised if the class could not be imported.  This could mean
            that the relevant obs package has not been setup.
        TypeError
            Raised if the class name retrieved is not a string.
        """
        dimensions = list(registry.queryDimensions("instrument", dataId={"instrument": name}))
        cls = dimensions[0].records["instrument"].class_name
        if not isinstance(cls, str):
            raise TypeError(f"Unexpected class name retrieved from {name} instrument dimension (got {cls})")
        instrument = doImport(cls)
        return instrument()

    def _registerFilters(self, registry):
        """Register the physical and abstract filter Dimension relationships.
        This should be called in the ``register`` implementation.

        Parameters
        ----------
        registry : `lsst.daf.butler.core.Registry`
            The registry to add dimensions to.
        """
        for filter in self.filterDefinitions:
            # fix for undefined abstract filters causing trouble in the registry:
            if filter.abstract_filter is None:
                abstract_filter = filter.physical_filter
            else:
                abstract_filter = filter.abstract_filter

            registry.insertDimensionData("physical_filter",
                                         {"instrument": self.getName(),
                                          "name": filter.physical_filter,
                                          "abstract_filter": abstract_filter
                                          })

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

    def writeCuratedCalibrations(self, butler):
        """Write human-curated calibration Datasets to the given Butler with
        the appropriate validity ranges.

        Parameters
        ----------
        butler : `lsst.daf.butler.Butler`
            Butler to use to store these calibrations.

        Notes
        -----
        Expected to be called from subclasses.  The base method calls
        ``writeCameraGeom`` and ``writeStandardTextCuratedCalibrations``.
        """
        self.writeCameraGeom(butler)
        self.writeStandardTextCuratedCalibrations(butler)

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

    def writeCameraGeom(self, butler):
        """Write the default camera geometry to the butler repository
        with an infinite validity range.

        Parameters
        ----------
        butler : `lsst.daf.butler.Butler`
            Butler to receive these calibration datasets.
        """

        datasetType = DatasetType("camera", ("instrument", "calibration_label"), "Camera",
                                  universe=butler.registry.dimensions)
        butler.registry.registerDatasetType(datasetType)
        unboundedDataId = addUnboundedCalibrationLabel(butler.registry, self.getName())
        camera = self.getCamera()
        butler.put(camera, datasetType, unboundedDataId)

    def writeStandardTextCuratedCalibrations(self, butler):
        """Write the set of standardized curated text calibrations to
        the repository.

        Parameters
        ----------
        butler : `lsst.daf.butler.Butler`
            Butler to receive these calibration datasets.
        """

        for datasetTypeName in self.standardCuratedDatasetTypes:
            # We need to define the dataset types.
            if datasetTypeName not in StandardCuratedCalibrationDatasetTypes:
                raise ValueError(f"DatasetType {datasetTypeName} not in understood list"
                                 f" [{'.'.join(StandardCuratedCalibrationDatasetTypes)}]")
            definition = StandardCuratedCalibrationDatasetTypes[datasetTypeName]
            datasetType = DatasetType(datasetTypeName,
                                      universe=butler.registry.dimensions,
                                      **definition)
            self._writeSpecificCuratedCalibrationDatasets(butler, datasetType)

    def _writeSpecificCuratedCalibrationDatasets(self, butler, datasetType):
        """Write standardized curated calibration datasets for this specific
        dataset type from an obs data package.

        Parameters
        ----------
        butler : `lsst.daf.butler.Butler`
            Gen3 butler in which to put the calibrations.
        datasetType : `lsst.daf.butler.DatasetType`
            Dataset type to be put.

        Notes
        -----
        This method scans the location defined in the ``obsDataPackageDir``
        class attribute for curated calibrations corresponding to the
        supplied dataset type.  The directory name in the data package must
        match the name of the dataset type. They are assumed to use the
        standard layout and can be read by
        `~lsst.pipe.tasks.read_curated_calibs.read_all` and provide standard
        metadata.
        """
        if self.obsDataPackageDir is None:
            # if there is no data package then there can't be datasets
            return

        calibPath = os.path.join(self.obsDataPackageDir, self.policyName,
                                 datasetType.name)

        if not os.path.exists(calibPath):
            return

        # Register the dataset type
        butler.registry.registerDatasetType(datasetType)

        # obs_base can't depend on pipe_tasks but concrete obs packages
        # can -- we therefore have to defer import
        from lsst.pipe.tasks.read_curated_calibs import read_all

        camera = self.getCamera()
        calibsDict = read_all(calibPath, camera)[0]  # second return is calib type
        endOfTime = TIMESPAN_MAX
        dimensionRecords = []
        datasetRecords = []
        for det in calibsDict:
            times = sorted([k for k in calibsDict[det]])
            calibs = [calibsDict[det][time] for time in times]
            times = [astropy.time.Time(t, format="datetime", scale="utc") for t in times]
            times += [endOfTime]
            for calib, beginTime, endTime in zip(calibs, times[:-1], times[1:]):
                md = calib.getMetadata()
                calibrationLabel = f"{datasetType.name}/{md['CALIBDATE']}/{md['DETECTOR']}"
                dataId = DataCoordinate.standardize(
                    universe=butler.registry.dimensions,
                    instrument=self.getName(),
                    calibration_label=calibrationLabel,
                    detector=md["DETECTOR"],
                )
                datasetRecords.append((calib, dataId))
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
            for calib, dataId in datasetRecords:
                butler.put(calib, datasetType, dataId)

    @abstractmethod
    def makeDataIdTranslatorFactory(self) -> TranslatorFactory:
        """Return a factory for creating Gen2->Gen3 data ID translators,
        specialized for this instrument.

        Derived class implementations should generally call
        `TranslatorFactory.addGenericInstrumentRules` with appropriate
        arguments, but are not required to (and may not be able to if their
        Gen2 raw data IDs are sufficiently different from the HSC/DECam/CFHT
        norm).

        Returns
        -------
        factory : `TranslatorFactory`.
            Factory for `Translator` objects.
        """
        raise NotImplementedError("Must be implemented by derived classes.")


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
        "group_name": obsInfo.exposure_group,
        "group_id": obsInfo.visit_id,
        "datetime_begin": obsInfo.datetime_begin,
        "datetime_end": obsInfo.datetime_end,
        "exposure_time": obsInfo.exposure_time.to_value("s"),
        "dark_time": obsInfo.dark_time.to_value("s"),
        "observation_type": obsInfo.observation_type,
        "physical_filter": obsInfo.physical_filter,
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
    entry["datetime_begin"] = TIMESPAN_MIN
    entry["datetime_end"] = TIMESPAN_MAX
    registry.insertDimensionData("calibration_label", entry)
    return registry.expandDataId(d)


def loadCamera(butler: Butler, dataId: DataId, *, collections: Any = None) -> Tuple[Camera, bool]:
    """Attempt to load versioned camera geometry from a butler, but fall back
    to obtaining a nominal camera from the `Instrument` class if that fails.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
        Butler instance to attempt to query for and load a ``camera`` dataset
        from.
    dataId : `dict` or `DataCoordinate`
        Data ID that identifies at least the ``instrument`` and ``exposure``
        dimensions.
    collections : Any, optional
        Collections to be searched, overriding ``self.butler.collections``.
        Can be any of the types supported by the ``collections`` argument
        to butler construction.

    Returns
    -------
    camera : `lsst.afw.cameraGeom.Camera`
        Camera object.
    versioned : `bool`
        If `True`, the camera was obtained from the butler and should represent
        a versioned camera from a calibration repository.  If `False`, no
        camera datasets were found, and the returned camera was produced by
        instantiating the appropriate `Instrument` class and calling
        `Instrument.getCamera`.
    """
    if collections is None:
        collections = butler.collections
    # Registry would do data ID expansion internally if we didn't do it first,
    # but we might want an expanded data ID ourselves later, so we do it here
    # to ensure it only happens once.
    # This will also catch problems with the data ID not having keys we need.
    dataId = butler.registry.expandDataId(dataId, graph=butler.registry.dimensions["exposure"].graph)
    cameraRefs = list(butler.registry.queryDatasets("camera", dataId=dataId, collections=collections,
                                                    deduplicate=True))
    if cameraRefs:
        assert len(cameraRefs) == 1, "Should be guaranteed by deduplicate=True above."
        return butler.getDirect(cameraRefs[0]), True
    instrument = Instrument.fromName(dataId["instrument"], butler.registry)
    return instrument.getCamera(), False
