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
from typing import Any, Optional, Set, Sequence, Tuple, TYPE_CHECKING
import astropy.time
from functools import lru_cache

from lsst.afw.cameraGeom import Camera
from lsst.daf.butler import (
    Butler,
    CollectionType,
    DataCoordinate,
    DataId,
    DatasetType,
    Timespan,
)
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
    "crosstalk": {"dimensions": ("instrument", "detector", "calibration_label"),
                  "storageClass": "CrosstalkCalib"},
}


class Instrument(metaclass=ABCMeta):
    """Base class for instrument-specific logic for the Gen3 Butler.

    Concrete instrument subclasses should be directly constructable with no
    arguments.
    """

    configPaths: Sequence[str] = ()
    """Paths to config files to read for specific Tasks.

    The paths in this list should contain files of the form `task.py`, for
    each of the Tasks that requires special configuration.
    """

    policyName: Optional[str] = None
    """Instrument specific name to use when locating a policy or configuration
    file in the file system."""

    obsDataPackage: Optional[str] = None
    """Name of the package containing the text curated calibration files.
    Usually a obs _data package.  If `None` no curated calibration files
    will be read. (`str`)"""

    standardCuratedDatasetTypes: Set[str] = frozenset(StandardCuratedCalibrationDatasetTypes)
    """The dataset types expected to be obtained from the obsDataPackage.

    These dataset types are all required to have standard definitions and
    must be known to the base class.  Clearing this list will prevent
    any of these calibrations from being stored. If a dataset type is not
    known to a specific instrument it can still be included in this list
    since the data package is the source of truth. (`set` of `str`)
    """

    additionalCuratedDatasetTypes: Set[str] = frozenset()
    """Curated dataset types specific to this particular instrument that do
    not follow the standard organization found in obs data packages.

    These are the instrument-specific dataset types written by
    `writeAdditionalCuratedCalibrations` in addition to the calibrations
    found in obs data packages that follow the standard scheme.
    (`set` of `str`)"""

    @property
    @abstractmethod
    def filterDefinitions(self):
        """`~lsst.obs.base.FilterDefinitionCollection`, defining the filters
        for this instrument.
        """
        return None

    def __init__(self):
        self.filterDefinitions.reset()
        self.filterDefinitions.defineFilters()

    @classmethod
    @abstractmethod
    def getName(cls):
        """Return the short (dimension) name for this instrument.

        This is not (in general) the same as the class name - it's what is used
        as the value of the "instrument" field in data IDs, and is usually an
        abbreviation of the full name.
        """
        raise NotImplementedError()

    @classmethod
    @lru_cache()
    def getCuratedCalibrationNames(cls) -> Set[str]:
        """Return the names of all the curated calibration dataset types.

        Returns
        -------
        names : `set` of `str`
            The dataset type names of all curated calibrations. This will
            include the standard curated calibrations even if the particular
            instrument does not support them.

        Notes
        -----
        The returned list does not indicate whether a particular dataset
        is present in the Butler repository, simply that these are the
        dataset types that are handled by ``writeCuratedCalibrations``.
        """

        # Camera is a special dataset type that is also handled as a
        # curated calibration.
        curated = {"camera"}

        # Make a cursory attempt to filter out curated dataset types
        # that are not present for this instrument
        for datasetTypeName in cls.standardCuratedDatasetTypes:
            calibPath = cls._getSpecificCuratedCalibrationPath(datasetTypeName)
            if calibPath is not None:
                curated.add(datasetTypeName)

        curated.update(cls.additionalCuratedDatasetTypes)
        return frozenset(curated)

    @abstractmethod
    def getCamera(self):
        """Retrieve the cameraGeom representation of this instrument.

        This is a temporary API that should go away once ``obs_`` packages have
        a standardized approach to writing versioned cameras to a Gen3 repo.
        """
        raise NotImplementedError()

    @abstractmethod
    def register(self, registry):
        """Insert instrument, physical_filter, and detector entries into a
        `Registry`.
        """
        raise NotImplementedError()

    @classmethod
    @lru_cache()
    def getObsDataPackageDir(cls):
        """The root of the obs data package that provides specializations for
        this instrument.

        returns
        -------
        dir : `str`
            The root of the relevat obs data package.
        """
        if cls.obsDataPackage is None:
            return None
        return getPackageDir(cls.obsDataPackage)

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
        records = list(registry.queryDimensionRecords("instrument", instrument=name))
        if not records:
            raise LookupError(f"No registered instrument with name '{name}'.")
        cls = records[0].class_name
        if not isinstance(cls, str):
            raise TypeError(f"Unexpected class name retrieved from {name} instrument dimension (got {cls})")
        instrument = doImport(cls)
        return instrument()

    @staticmethod
    def importAll(registry: Registry) -> None:
        """Import all the instruments known to this registry.

        This will ensure that all metadata translators have been registered.

        Parameters
        ----------
        registry : `lsst.daf.butler.Registry`
            Butler registry to query to find the information.

        Notes
        -----
        It is allowed for a particular instrument class to fail on import.
        This might simply indicate that a particular obs package has
        not been setup.
        """
        records = list(registry.queryDimensionRecords("instrument"))
        for record in records:
            cls = record.class_name
            try:
                doImport(cls)
            except Exception:
                pass

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

    def writeCuratedCalibrations(self, butler, run=None):
        """Write human-curated calibration Datasets to the given Butler with
        the appropriate validity ranges.

        Parameters
        ----------
        butler : `lsst.daf.butler.Butler`
            Butler to use to store these calibrations.
        run : `str`
            Run to use for this collection of calibrations. If `None` the
            collection name is worked out automatically from the instrument
            name and other metadata.

        Notes
        -----
        Expected to be called from subclasses.  The base method calls
        ``writeCameraGeom`` and ``writeStandardTextCuratedCalibrations``.
        """
        # Need to determine the run for ingestion based on the instrument
        # name and eventually the data package version. The camera geom
        # is currently special in that it is not in the _data package.
        if run is None:
            run = self.makeCollectionName("calib")
        butler.registry.registerCollection(run, type=CollectionType.RUN)
        self.writeCameraGeom(butler, run=run)
        self.writeStandardTextCuratedCalibrations(butler, run=run)
        self.writeAdditionalCuratedCalibrations(butler, run=run)

    def writeAdditionalCuratedCalibrations(self, butler, run=None):
        """Write additional curated calibrations that might be instrument
        specific and are not part of the standard set.

        Default implementation does nothing.

        Parameters
        ----------
        butler : `lsst.daf.butler.Butler`
            Butler to use to store these calibrations.
        run : `str`, optional
            Name of the run to use to override the default run associated
            with this Butler.
        """
        return

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

    def writeCameraGeom(self, butler, run=None):
        """Write the default camera geometry to the butler repository
        with an infinite validity range.

        Parameters
        ----------
        butler : `lsst.daf.butler.Butler`
            Butler to receive these calibration datasets.
        run : `str`, optional
            Name of the run to use to override the default run associated
            with this Butler.
        """

        datasetType = DatasetType("camera", ("instrument", "calibration_label"), "Camera",
                                  universe=butler.registry.dimensions)
        butler.registry.registerDatasetType(datasetType)
        unboundedDataId = addUnboundedCalibrationLabel(butler.registry, self.getName())
        camera = self.getCamera()
        butler.put(camera, datasetType, unboundedDataId, run=run)

    def writeStandardTextCuratedCalibrations(self, butler, run=None):
        """Write the set of standardized curated text calibrations to
        the repository.

        Parameters
        ----------
        butler : `lsst.daf.butler.Butler`
            Butler to receive these calibration datasets.
        run : `str`, optional
            Name of the run to use to override the default run associated
            with this Butler.
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
            self._writeSpecificCuratedCalibrationDatasets(butler, datasetType, run=run)

    @classmethod
    def _getSpecificCuratedCalibrationPath(cls, datasetTypeName):
        """Return the path of the curated calibration directory.

        Parameters
        ----------
        datasetTypeName : `str`
            The name of the standard dataset type to find.

        Returns
        -------
        path : `str`
            The path to the standard curated data directory.  `None` if the
            dataset type is not found or the obs data package is not
            available.
        """
        if cls.getObsDataPackageDir() is None:
            # if there is no data package then there can't be datasets
            return None

        calibPath = os.path.join(cls.getObsDataPackageDir(), cls.policyName,
                                 datasetTypeName)

        if os.path.exists(calibPath):
            return calibPath

        return None

    def _writeSpecificCuratedCalibrationDatasets(self, butler, datasetType, run=None):
        """Write standardized curated calibration datasets for this specific
        dataset type from an obs data package.

        Parameters
        ----------
        butler : `lsst.daf.butler.Butler`
            Gen3 butler in which to put the calibrations.
        datasetType : `lsst.daf.butler.DatasetType`
            Dataset type to be put.
        run : `str`, optional
            Name of the run to use to override the default run associated
            with this Butler.

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
        calibPath = self._getSpecificCuratedCalibrationPath(datasetType.name)
        if calibPath is None:
            return

        # Register the dataset type
        butler.registry.registerDatasetType(datasetType)

        # obs_base can't depend on pipe_tasks but concrete obs packages
        # can -- we therefore have to defer import
        from lsst.pipe.tasks.read_curated_calibs import read_all

        camera = self.getCamera()
        calibsDict = read_all(calibPath, camera)[0]  # second return is calib type
        dimensionRecords = []
        datasetRecords = []
        for det in calibsDict:
            times = sorted([k for k in calibsDict[det]])
            calibs = [calibsDict[det][time] for time in times]
            times = [astropy.time.Time(t, format="datetime", scale="utc") for t in times]
            times += [None]
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
                    "timespan": Timespan(beginTime, endTime),
                })

        # Second loop actually does the inserts and filesystem writes.
        with butler.transaction():
            butler.registry.insertDimensionData("calibration_label", *dimensionRecords)
            # TODO: vectorize these puts, once butler APIs for that become
            # available.
            for calib, dataId in datasetRecords:
                butler.put(calib, datasetType, dataId, run=run)

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

    @classmethod
    def makeDefaultRawIngestRunName(cls) -> str:
        """Make the default instrument-specific run collection string for raw
        data ingest.

        Returns
        -------
        coll : `str`
            Run collection name to be used as the default for ingestion of
            raws.
        """
        return cls.makeCollectionName("raw/all")

    @classmethod
    def makeCollectionName(cls, label: str) -> str:
        """Get the instrument-specific collection string to use as derived
        from the supplied label.

        Parameters
        ----------
        label : `str`
            String to be combined with the instrument name to form a
            collection name.

        Returns
        -------
        name : `str`
            Collection name to use that includes the instrument name.
        """
        return f"{cls.getName()}/{label}"


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

    ra, dec, sky_angle, zenith_angle = (None, None, None, None)
    if obsInfo.tracking_radec is not None:
        icrs = obsInfo.tracking_radec.icrs
        ra = icrs.ra.degree
        dec = icrs.dec.degree
        if obsInfo.boresight_rotation_coord == "sky":
            sky_angle = obsInfo.boresight_rotation_angle.degree
    if obsInfo.altaz_begin is not None:
        zenith_angle = obsInfo.altaz_begin.zen.degree

    return dimension.RecordClass(
        instrument=obsInfo.instrument,
        id=obsInfo.exposure_id,
        name=obsInfo.observation_id,
        group_name=obsInfo.exposure_group,
        group_id=obsInfo.visit_id,
        datetime_begin=obsInfo.datetime_begin,
        datetime_end=obsInfo.datetime_end,
        exposure_time=obsInfo.exposure_time.to_value("s"),
        dark_time=obsInfo.dark_time.to_value("s"),
        observation_type=obsInfo.observation_type,
        physical_filter=obsInfo.physical_filter,
        science_program=obsInfo.science_program,
        target_name=obsInfo.object,
        tracking_ra=ra,
        tracking_dec=dec,
        sky_angle=sky_angle,
        zenith_angle=zenith_angle,
    )


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
    entry["timespan"] = Timespan(None, None)
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
