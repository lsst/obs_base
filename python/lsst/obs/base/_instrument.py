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

__all__ = ("Instrument", "makeExposureRecordFromObsInfo", "loadCamera")

import os.path
from abc import ABCMeta, abstractmethod
from collections import defaultdict
import datetime
from typing import Any, Optional, Set, Sequence, Tuple, TYPE_CHECKING, Union
from functools import lru_cache

import astropy.time

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
    "defects": {"dimensions": ("instrument", "detector"), "storageClass": "Defects"},
    "qe_curve": {"dimensions": ("instrument", "detector"), "storageClass": "QECurve"},
    "crosstalk": {"dimensions": ("instrument", "detector"), "storageClass": "CrosstalkCalib"},
    "linearizer": {"dimensions": ("instrument", "detector"), "storageClass": "Linearizer"},
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

        This is a temporary API that should go away once ``obs`` packages have
        a standardized approach to writing versioned cameras to a Gen3 repo.
        """
        raise NotImplementedError()

    @abstractmethod
    def register(self, registry):
        """Insert instrument, physical_filter, and detector entries into a
        `Registry`.

        Implementations should guarantee that registration is atomic (the
        registry should not be modified if any error occurs) and idempotent at
        the level of individual dimension entries; new detectors and filters
        should be added, but changes to any existing record should not be.
        This can generally be achieved via a block like::

            with registry.transaction():
                registry.syncDimensionData("instrument", ...)
                registry.syncDimensionData("detector", ...)
                self.registerFilters(registry)

        Raises
        ------
        lsst.daf.butler.registry.ConflictingDefinitionError
            Raised if any existing record has the same key but a different
            definition as one being registered.
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
        This should be called in the `register` implementation, within
        a transaction context manager block.

        Parameters
        ----------
        registry : `lsst.daf.butler.core.Registry`
            The registry to add dimensions to.
        """
        for filter in self.filterDefinitions:
            # fix for undefined abstract filters causing trouble in the
            # registry:
            if filter.band is None:
                band = filter.physical_filter
            else:
                band = filter.band

            registry.syncDimensionData("physical_filter",
                                       {"instrument": self.getName(),
                                        "name": filter.physical_filter,
                                        "band": band
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

    def writeCuratedCalibrations(self, butler: Butler, collection: Optional[str] = None,
                                 labels: Sequence[str] = ()) -> None:
        """Write human-curated calibration Datasets to the given Butler with
        the appropriate validity ranges.

        Parameters
        ----------
        butler : `lsst.daf.butler.Butler`
            Butler to use to store these calibrations.
        collection : `str`, optional
            Name to use for the calibration collection that associates all
            datasets with a validity range.  If this collection already exists,
            it must be a `~CollectionType.CALIBRATION` collection, and it must
            not have any datasets that would conflict with those inserted by
            this method.  If `None`, a collection name is worked out
            automatically from the instrument name and other metadata by
            calling ``makeCalibrationCollectionName``, but this
            default name may not work well for long-lived repositories unless
            ``labels`` is also provided (and changed every time curated
            calibrations are ingested).
        labels : `Sequence` [ `str` ], optional
            Extra strings to include in collection names, after concatenating
            them with the standard collection name delimeter.  If provided,
            these are inserted into the names of the `~CollectionType.RUN`
            collections that datasets are inserted directly into, as well the
            `~CollectionType.CALIBRATION` collection if it is generated
            automatically (i.e. if ``collection is None``).  Usually this is
            just the name of the ticket on which the calibration collection is
            being created.

        Notes
        -----
        Expected to be called from subclasses.  The base method calls
        ``writeCameraGeom``, ``writeStandardTextCuratedCalibrations``,
        and ``writeAdditionalCuratdCalibrations``.
        """
        # Delegate registration of collections (and creating names for them)
        # to other methods so they can be called independently with the same
        # preconditions.  Collection registration is idempotent, so this is
        # safe, and while it adds a bit of overhead, as long as it's one
        # registration attempt per method (not per dataset or dataset type),
        # that's negligible.
        self.writeCameraGeom(butler, collection, labels=labels)
        self.writeStandardTextCuratedCalibrations(butler, collection, labels=labels)
        self.writeAdditionalCuratedCalibrations(butler, collection, labels=labels)

    def writeAdditionalCuratedCalibrations(self, butler: Butler, collection: Optional[str] = None,
                                           labels: Sequence[str] = ()) -> None:
        """Write additional curated calibrations that might be instrument
        specific and are not part of the standard set.

        Default implementation does nothing.

        Parameters
        ----------
        butler : `lsst.daf.butler.Butler`
            Butler to use to store these calibrations.
        collection : `str`, optional
            Name to use for the calibration collection that associates all
            datasets with a validity range.  If this collection already exists,
            it must be a `~CollectionType.CALIBRATION` collection, and it must
            not have any datasets that would conflict with those inserted by
            this method.  If `None`, a collection name is worked out
            automatically from the instrument name and other metadata by
            calling ``makeCalibrationCollectionName``, but this
            default name may not work well for long-lived repositories unless
            ``labels`` is also provided (and changed every time curated
            calibrations are ingested).
        labels : `Sequence` [ `str` ], optional
            Extra strings to include in collection names, after concatenating
            them with the standard collection name delimeter.  If provided,
            these are inserted into the names of the `~CollectionType.RUN`
            collections that datasets are inserted directly into, as well the
            `~CollectionType.CALIBRATION` collection if it is generated
            automatically (i.e. if ``collection is None``).  Usually this is
            just the name of the ticket on which the calibration collection is
            being created.
        """
        return

    def writeCameraGeom(self, butler: Butler, collection: Optional[str] = None,
                        labels: Sequence[str] = ()) -> None:
        """Write the default camera geometry to the butler repository and
        associate it with the appropriate validity range in a calibration
        collection.

        Parameters
        ----------
        butler : `lsst.daf.butler.Butler`
            Butler to use to store these calibrations.
        collection : `str`, optional
            Name to use for the calibration collection that associates all
            datasets with a validity range.  If this collection already exists,
            it must be a `~CollectionType.CALIBRATION` collection, and it must
            not have any datasets that would conflict with those inserted by
            this method.  If `None`, a collection name is worked out
            automatically from the instrument name and other metadata by
            calling ``makeCalibrationCollectionName``, but this
            default name may not work well for long-lived repositories unless
            ``labels`` is also provided (and changed every time curated
            calibrations are ingested).
        labels : `Sequence` [ `str` ], optional
            Extra strings to include in collection names, after concatenating
            them with the standard collection name delimeter.  If provided,
            these are inserted into the names of the `~CollectionType.RUN`
            collections that datasets are inserted directly into, as well the
            `~CollectionType.CALIBRATION` collection if it is generated
            automatically (i.e. if ``collection is None``).  Usually this is
            just the name of the ticket on which the calibration collection is
            being created.
        """
        if collection is None:
            collection = self.makeCalibrationCollectionName(*labels)
        butler.registry.registerCollection(collection, type=CollectionType.CALIBRATION)
        run = self.makeUnboundedCalibrationRunName(*labels)
        butler.registry.registerRun(run)
        datasetType = DatasetType("camera", ("instrument",), "Camera", isCalibration=True,
                                  universe=butler.registry.dimensions)
        butler.registry.registerDatasetType(datasetType)
        camera = self.getCamera()
        ref = butler.put(camera, datasetType, {"instrument": self.getName()}, run=run)
        butler.registry.certify(collection, [ref], Timespan(begin=None, end=None))

    def writeStandardTextCuratedCalibrations(self, butler: Butler, collection: Optional[str] = None,
                                             labels: Sequence[str] = ()) -> None:
        """Write the set of standardized curated text calibrations to
        the repository.

        Parameters
        ----------
        butler : `lsst.daf.butler.Butler`
            Butler to receive these calibration datasets.
        collection : `str`, optional
            Name to use for the calibration collection that associates all
            datasets with a validity range.  If this collection already exists,
            it must be a `~CollectionType.CALIBRATION` collection, and it must
            not have any datasets that would conflict with those inserted by
            this method.  If `None`, a collection name is worked out
            automatically from the instrument name and other metadata by
            calling ``makeCalibrationCollectionName``, but this
            default name may not work well for long-lived repositories unless
            ``labels`` is also provided (and changed every time curated
            calibrations are ingested).
        labels : `Sequence` [ `str` ], optional
            Extra strings to include in collection names, after concatenating
            them with the standard collection name delimeter.  If provided,
            these are inserted into the names of the `~CollectionType.RUN`
            collections that datasets are inserted directly into, as well the
            `~CollectionType.CALIBRATION` collection if it is generated
            automatically (i.e. if ``collection is None``).  Usually this is
            just the name of the ticket on which the calibration collection is
            being created.
        """
        if collection is None:
            collection = self.makeCalibrationCollectionName(*labels)
        butler.registry.registerCollection(collection, type=CollectionType.CALIBRATION)
        runs = set()
        for datasetTypeName in self.standardCuratedDatasetTypes:
            # We need to define the dataset types.
            if datasetTypeName not in StandardCuratedCalibrationDatasetTypes:
                raise ValueError(f"DatasetType {datasetTypeName} not in understood list"
                                 f" [{'.'.join(StandardCuratedCalibrationDatasetTypes)}]")
            definition = StandardCuratedCalibrationDatasetTypes[datasetTypeName]
            datasetType = DatasetType(datasetTypeName,
                                      universe=butler.registry.dimensions,
                                      isCalibration=True,
                                      **definition)
            self._writeSpecificCuratedCalibrationDatasets(butler, datasetType, collection, runs=runs,
                                                          labels=labels)

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

    def _writeSpecificCuratedCalibrationDatasets(self, butler: Butler, datasetType: DatasetType,
                                                 collection: str, runs: Set[str], labels: Sequence[str]):
        """Write standardized curated calibration datasets for this specific
        dataset type from an obs data package.

        Parameters
        ----------
        butler : `lsst.daf.butler.Butler`
            Gen3 butler in which to put the calibrations.
        datasetType : `lsst.daf.butler.DatasetType`
            Dataset type to be put.
        collection : `str`
            Name of the `~CollectionType.CALIBRATION` collection that
            associates all datasets with validity ranges.  Must have been
            registered prior to this call.
        runs : `set` [ `str` ]
            Names of runs that have already been registered by previous calls
            and need not be registered again.  Should be updated by this
            method as new runs are registered.
        labels : `Sequence` [ `str` ]
            Extra strings to include in run names when creating them from
            ``CALIBDATE`` metadata, via calls to `makeCuratedCalibrationName`.
            Usually this is the name of the ticket on which the calibration
            collection is being created.

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

        # Read calibs, registering a new run for each CALIBDATE as needed.
        # We try to avoid registering runs multiple times as an optimization
        # by putting them in the ``runs`` set that was passed in.
        camera = self.getCamera()
        calibsDict = read_all(calibPath, camera)[0]  # second return is calib type
        datasetRecords = []
        for det in calibsDict:
            times = sorted([k for k in calibsDict[det]])
            calibs = [calibsDict[det][time] for time in times]
            times = [astropy.time.Time(t, format="datetime", scale="utc") for t in times]
            times += [None]
            for calib, beginTime, endTime in zip(calibs, times[:-1], times[1:]):
                md = calib.getMetadata()
                run = self.makeCuratedCalibrationRunName(md['CALIBDATE'], *labels)
                if run not in runs:
                    butler.registry.registerRun(run)
                    runs.add(run)
                dataId = DataCoordinate.standardize(
                    universe=butler.registry.dimensions,
                    instrument=self.getName(),
                    detector=md["DETECTOR"],
                )
                datasetRecords.append((calib, dataId, run, Timespan(beginTime, endTime)))

        # Second loop actually does the inserts and filesystem writes.  We
        # first do a butler.put on each dataset, inserting it into the run for
        # its calibDate.  We remember those refs and group them by timespan, so
        # we can vectorize the certify calls as much as possible.
        refsByTimespan = defaultdict(list)
        with butler.transaction():
            for calib, dataId, run, timespan in datasetRecords:
                refsByTimespan[timespan].append(butler.put(calib, datasetType, dataId, run=run))
            for timespan, refs in refsByTimespan.items():
                butler.registry.certify(collection, refs, timespan)

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

    @staticmethod
    def formatCollectionTimestamp(timestamp: Union[str, datetime.datetime]) -> str:
        """Format a timestamp for use in a collection name.

        Parameters
        ----------
        timestamp : `str` or `datetime.datetime`
            Timestamp to format.  May be a date or datetime string in extended
            ISO format (assumed UTC), with or without a timezone specifier, a
            datetime string in basic ISO format with a timezone specifier, a
            naive `datetime.datetime` instance (assumed UTC) or a
            timezone-aware `datetime.datetime` instance (converted to UTC).
            This is intended to cover all forms that string ``CALIBDATE``
            metadata values have taken in the past, as well as the format this
            method itself writes out (to enable round-tripping).

        Returns
        -------
        formatted : `str`
            Standardized string form for the timestamp.
        """
        if isinstance(timestamp, str):
            if "-" in timestamp:
                # extended ISO format, with - and : delimiters
                timestamp = datetime.datetime.fromisoformat(timestamp)
            else:
                # basic ISO format, with no delimiters (what this method
                # returns)
                timestamp = datetime.datetime.strptime(timestamp, "%Y%m%dT%H%M%S%z")
        if not isinstance(timestamp, datetime.datetime):
            raise TypeError(f"Unexpected date/time object: {timestamp!r}.")
        if timestamp.tzinfo is not None:
            timestamp = timestamp.astimezone(datetime.timezone.utc)
        return f"{timestamp:%Y%m%dT%H%M%S}Z"

    @staticmethod
    def makeCollectionTimestamp() -> str:
        """Create a timestamp string for use in a collection name from the
        current time.

        Returns
        -------
        formatted : `str`
            Standardized string form of the current time.
        """
        return Instrument.formatCollectionTimestamp(datetime.datetime.now(tz=datetime.timezone.utc))

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
        return cls.makeCollectionName("raw", "all")

    @classmethod
    def makeUnboundedCalibrationRunName(cls, *labels: str) -> str:
        """Make a RUN collection name appropriate for inserting calibration
        datasets whose validity ranges are unbounded.

        Parameters
        ----------
        *labels : `str`
            Extra strings to be included in the base name, using the default
            delimiter for collection names.  Usually this is the name of the
            ticket on which the calibration collection is being created.

        Returns
        -------
        name : `str`
            Run collection name.
        """
        return cls.makeCollectionName("calib", *labels, "unbounded")

    @classmethod
    def makeCuratedCalibrationRunName(cls, calibDate: str, *labels: str) -> str:
        """Make a RUN collection name appropriate for inserting curated
        calibration datasets with the given ``CALIBDATE`` metadata value.

        Parameters
        ----------
        calibDate : `str`
            The ``CALIBDATE`` metadata value.
        *labels : `str`
            Strings to be included in the collection name (before
            ``calibDate``, but after all other terms), using the default
            delimiter for collection names.  Usually this is the name of the
            ticket on which the calibration collection is being created.

        Returns
        -------
        name : `str`
            Run collection name.
        """
        return cls.makeCollectionName("calib", *labels, "curated", cls.formatCollectionTimestamp(calibDate))

    @classmethod
    def makeCalibrationCollectionName(cls, *labels: str) -> str:
        """Make a CALIBRATION collection name appropriate for associating
        calibration datasets with validity ranges.

        Parameters
        ----------
        *labels : `str`
            Strings to be appended to the base name, using the default
            delimiter for collection names.  Usually this is the name of the
            ticket on which the calibration collection is being created.

        Returns
        -------
        name : `str`
            Calibration collection name.
        """
        return cls.makeCollectionName("calib", *labels)

    @staticmethod
    def makeRefCatCollectionName(*labels: str) -> str:
        """Return a global (not instrument-specific) name for a collection that
        holds reference catalogs.

        With no arguments, this returns the name of the collection that holds
        all reference catalogs (usually a ``CHAINED`` collection, at least in
        long-lived repos that may contain more than one reference catalog).

        Parameters
        ----------
        *labels : `str`
            Strings to be added to the global collection name, in order to
            define a collection name for one or more reference catalogs being
            ingested at the same time.

        Returns
        -------
        name : `str`
            Collection name.

        Notes
        -----
        This is a ``staticmethod``, not a ``classmethod``, because it should
        be the same for all instruments.
        """
        return "/".join(("refcats",) + labels)

    @classmethod
    def makeUmbrellaCollectionName(cls) -> str:
        """Return the name of the umbrella ``CHAINED`` collection for this
        instrument that combines all standard recommended input collections.

        This method should almost never be overridden by derived classes.

        Returns
        -------
        name : `str`
            Name for the umbrella collection.
        """
        return cls.makeCollectionName("defaults")

    @classmethod
    def makeCollectionName(cls, *labels: str) -> str:
        """Get the instrument-specific collection string to use as derived
        from the supplied labels.

        Parameters
        ----------
        *labels : `str`
            Strings to be combined with the instrument name to form a
            collection name.

        Returns
        -------
        name : `str`
            Collection name to use that includes the instrument name.
        """
        return "/".join((cls.getName(),) + labels)


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
        obs_id=obsInfo.observation_id,
        group_name=obsInfo.exposure_group,
        group_id=obsInfo.visit_id,
        datetime_begin=obsInfo.datetime_begin,
        datetime_end=obsInfo.datetime_end,
        exposure_time=obsInfo.exposure_time.to_value("s"),
        # we are not mandating that dark_time be calculable
        dark_time=obsInfo.dark_time.to_value("s") if obsInfo.dark_time is not None else None,
        observation_type=obsInfo.observation_type,
        observation_reason=obsInfo.observation_reason,
        day_obs=obsInfo.observing_day,
        seq_num=obsInfo.observation_counter,
        physical_filter=obsInfo.physical_filter,
        science_program=obsInfo.science_program,
        target_name=obsInfo.object,
        tracking_ra=ra,
        tracking_dec=dec,
        sky_angle=sky_angle,
        zenith_angle=zenith_angle,
    )


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
    try:
        cameraRef = butler.get("camera", dataId=dataId, collections=collections)
        return cameraRef, True
    except LookupError:
        pass
    instrument = Instrument.fromName(dataId["instrument"], butler.registry)
    return instrument.getCamera(), False
