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
from abc import abstractmethod
from collections import defaultdict
from functools import lru_cache
from typing import (
    TYPE_CHECKING,
    AbstractSet,
    Any,
    FrozenSet,
    Iterable,
    Mapping,
    Optional,
    Sequence,
    Set,
    Tuple,
    Union,
)

import astropy.time
from lsst.afw.cameraGeom import Camera, CameraBuilder, Detector, DetectorBuilder
from lsst.daf.butler import (
    Butler,
    CollectionType,
    DataCoordinate,
    DataId,
    DatasetRef,
    DatasetType,
    DimensionRecord,
    DimensionUniverse,
    Timespan,
)
from lsst.daf.butler.registry import DataIdError
from lsst.pipe.base import Instrument as InstrumentBase
from lsst.utils import getPackageDir

from ._camera_loaders import CameraLoader, DetectorLoader

if TYPE_CHECKING:
    from astro_metadata_translator import ObservationInfo
    from lsst.daf.butler import Registry

    from .filters import FilterDefinitionCollection
    from .gen2to3 import TranslatorFactory

# To be a standard text curated calibration means that we use a
# standard definition for the corresponding DatasetType.
StandardCuratedCalibrationDatasetTypes = {
    "defects": {"dimensions": ("instrument", "detector"), "storageClass": "Defects"},
    "qe_curve": {"dimensions": ("instrument", "detector"), "storageClass": "QECurve"},
    "crosstalk": {"dimensions": ("instrument", "detector"), "storageClass": "CrosstalkCalib"},
    "linearizer": {"dimensions": ("instrument", "detector"), "storageClass": "Linearizer"},
    "bfk": {"dimensions": ("instrument", "detector"), "storageClass": "BrighterFatterKernel"},
}


class Instrument(InstrumentBase):
    """Rubin-specified base for instrument-specific logic for the Gen3 Butler.

    Parameters
    ----------
    collection_prefix : `str`, optional
        Prefix for collection names to use instead of the intrument's own name.
        This is primarily for use in simulated-data repositories, where the
        instrument name may not be necessary and/or sufficient to distinguish
        between collections.

    Notes
    -----
    Concrete instrument subclasses must have the same construction signature as
    the base class.
    """

    policyName: Optional[str] = None
    """Instrument specific name to use when locating a policy or configuration
    file in the file system."""

    obsDataPackage: Optional[str] = None
    """Name of the package containing the text curated calibration files.
    Usually a obs _data package.  If `None` no curated calibration files
    will be read. (`str`)"""

    standardCuratedDatasetTypes: AbstractSet[str] = frozenset(StandardCuratedCalibrationDatasetTypes)
    """The dataset types expected to be obtained from the obsDataPackage.

    These dataset types are all required to have standard definitions and
    must be known to the base class.  Clearing this list will prevent
    any of these calibrations from being stored. If a dataset type is not
    known to a specific instrument it can still be included in this list
    since the data package is the source of truth. (`set` of `str`)
    """

    additionalCuratedDatasetTypes: AbstractSet[str] = frozenset()
    """Curated dataset types specific to this particular instrument that do
    not follow the standard organization found in obs data packages.

    These are the instrument-specific dataset types written by
    `writeAdditionalCuratedCalibrations` in addition to the calibrations
    found in obs data packages that follow the standard scheme.
    (`set` of `str`)"""

    @property
    @abstractmethod
    def filterDefinitions(self) -> FilterDefinitionCollection:
        """`~lsst.obs.base.FilterDefinitionCollection`, defining the filters
        for this instrument.
        """
        raise NotImplementedError()

    def __init__(self, collection_prefix: Optional[str] = None):
        super().__init__(collection_prefix=collection_prefix)
        self.filterDefinitions.reset()
        self.filterDefinitions.defineFilters()

    @classmethod
    @lru_cache()
    def getCuratedCalibrationNames(cls) -> FrozenSet[str]:
        """Return the names of all the curated calibration dataset types.

        Returns
        -------
        names : `frozenset` of `str`
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
    def getCamera(self) -> Camera:
        """Retrieve the cameraGeom representation of this instrument.

        This is a temporary API that should go away once ``obs`` packages have
        a standardized approach to writing versioned cameras to a Gen3 repo.
        """
        raise NotImplementedError()

    @classmethod
    @lru_cache()
    def getObsDataPackageDir(cls) -> Optional[str]:
        """The root of the obs data package that provides specializations for
        this instrument.

        returns
        -------
        dir : `str` or `None`
            The root of the relevant obs data package, or `None` if this
            instrument does not have one.
        """
        if cls.obsDataPackage is None:
            return None
        return getPackageDir(cls.obsDataPackage)

    def _registerFilters(self, registry: Registry, update: bool = False) -> None:
        """Register the physical and abstract filter Dimension relationships.
        This should be called in the `register` implementation, within
        a transaction context manager block.

        Parameters
        ----------
        registry : `lsst.daf.butler.core.Registry`
            The registry to add dimensions to.
        update : `bool`, optional
            If `True` (`False` is default), update existing records if they
            differ from the new ones.
        """
        for filter in self.filterDefinitions:
            # fix for undefined abstract filters causing trouble in the
            # registry:
            if filter.band is None:
                band = filter.physical_filter
            else:
                band = filter.band

            registry.syncDimensionData(
                "physical_filter",
                {"instrument": self.getName(), "name": filter.physical_filter, "band": band},
                update=update,
            )

    def writeCuratedCalibrations(
        self, butler: Butler, collection: Optional[str] = None, labels: Sequence[str] = ()
    ) -> None:
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

    def writeAdditionalCuratedCalibrations(
        self, butler: Butler, collection: Optional[str] = None, labels: Sequence[str] = ()
    ) -> None:
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

    def writeCameraGeom(
        self, butler: Butler, collection: Optional[str] = None, labels: Sequence[str] = ()
    ) -> None:
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
        datasetType = DatasetType(
            "camera", ("instrument",), "Camera", isCalibration=True, universe=butler.registry.dimensions
        )
        butler.registry.registerDatasetType(datasetType)
        camera = self.getCamera()
        ref = butler.put(camera, datasetType, {"instrument": self.getName()}, run=run)
        butler.registry.certify(collection, [ref], Timespan(begin=None, end=None))

    def writeStandardTextCuratedCalibrations(
        self, butler: Butler, collection: Optional[str] = None, labels: Sequence[str] = ()
    ) -> None:
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
        runs: Set[str] = set()
        for datasetTypeName in self.standardCuratedDatasetTypes:
            # We need to define the dataset types.
            if datasetTypeName not in StandardCuratedCalibrationDatasetTypes:
                raise ValueError(
                    f"DatasetType {datasetTypeName} not in understood list"
                    f" [{'.'.join(StandardCuratedCalibrationDatasetTypes)}]"
                )
            definition = StandardCuratedCalibrationDatasetTypes[datasetTypeName]
            datasetType = DatasetType(
                datasetTypeName,
                universe=butler.registry.dimensions,
                isCalibration=True,
                # MyPy should be able to figure out that the kwargs here have
                # the right types, but it can't.
                **definition,  # type: ignore
            )
            self._writeSpecificCuratedCalibrationDatasets(
                butler, datasetType, collection, runs=runs, labels=labels
            )

    @classmethod
    def _getSpecificCuratedCalibrationPath(cls, datasetTypeName: str) -> Optional[str]:
        """Return the path of the curated calibration directory.

        Parameters
        ----------
        datasetTypeName : `str`
            The name of the standard dataset type to find.

        Returns
        -------
        path : `str` or `None`
            The path to the standard curated data directory.  `None` if the
            dataset type is not found or the obs data package is not
            available.
        """
        data_package_dir = cls.getObsDataPackageDir()
        if data_package_dir is None:
            # if there is no data package then there can't be datasets
            return None

        if cls.policyName is None:
            raise TypeError(f"Instrument {cls.getName()} has an obs data package but no policy name.")

        calibPath = os.path.join(data_package_dir, cls.policyName, datasetTypeName)

        if os.path.exists(calibPath):
            return calibPath

        return None

    def _writeSpecificCuratedCalibrationDatasets(
        self, butler: Butler, datasetType: DatasetType, collection: str, runs: Set[str], labels: Sequence[str]
    ) -> None:
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
                run = self.makeCuratedCalibrationRunName(md["CALIBDATE"], *labels)
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

    def load_detector(
        self,
        arg: Union[str, int, DatasetRef, Detector, DetectorBuilder],
        /,
        butler: Butler,
        camera: Union[str, DatasetType, DatasetRef, Camera, CameraBuilder] = "camera",
        data_id: Optional[DataId] = None,
        collections: Optional[Iterable[str]] = None,
        **kwargs: Any,
    ) -> Detector:
        """Load a detector from a butler repository and update it with all
        appropriate versioned calibrations for this instrument.

        Parameters
        ----------
        arg : `str`, `int`, `DatasetRef`, `Detector`, `DetectorBuilder`
            If a `DatasetRef`, `Detector`, or `DetectorBuilder`, an object that
            can be used or loaded to provide the detector directly (``camera``
            will be ignored).
            If an `int` or `str`, a detector identifier to use with ``camera``.
        butler : `Butler`
            Butler client to read from.  If not initialized with the desired
            collection search path, the ``collections`` argument must be
            provided.
        camera : `str`, `DatasetType`, `DatasetRef`, `Camera`, \
                `CameraBuilder`, optional
            Used to obtain a `Camera` or `CameraBuilder` as a way to get a
            `DetectorBuilder`.  If a `Camera` is given, the detector is
            extracted and then a `Detector.rebuild` is called, instead of
            rebuilding the full camera.  Defaults to "camera" (as the nominal
            dataset type name).
        data_id : `dict` or `lsst.daf.butler.DataCoordinate`, optional
            Data ID values that identify the ``exposure`` and/or ``visit``
            dimension(s).  May include any extended keys supported by
            `Butler.get`, such as ``day_obs``.  Should not identify the
            ``detector`` dimension.
        collections : `Iterable` [ `str` ], optional
            Collections to search, in order.  If not provided, the butler used
            to initialize this loader must already have the right collections.
        **kwargs
            Additional keyword arguments that are in the keys of
            `detector_calibrations` will be interpreted as overrides of the
            calibrations to apply; see `DetectorLoader.update` for details.
            Other keyword arguments are interpreted as data ID key-value
            pairs and forwarded directly to `Butler.get`.

        Returns
        -------
        detector : `Detector`
            The updated detector.

        Examples
        --------
        Read detector 11 from the nominal camera and apply all default
        calibrations valid during the observation of exposure 500::

            detector = my_instrument.load_detector(11, butler, exposure=500)

        The same, but use a string detector name and day_obs/seq_num for
        exposure, and don't attempt to apply gains from PTC:

            detector = my_instrument.load_detector(
                "R12S21",
                butler,
                day_obs=20250129,
                seq_num=1261,
                ptc=False,
            )
        """
        update_kwargs = {}
        for k, v in self.detector_calibrations.items():
            update_kwargs[k] = kwargs.pop(k, v)
        loader = DetectorLoader(
            arg,
            butler,
            camera=camera,
            data_id=data_id,
            collections=collections,
            **kwargs,
            instrument=self.getName(),
        )
        loader.update(collections=collections, **update_kwargs)
        return loader.finish()

    def load_camera(
        self,
        butler: Butler,
        camera: Union[str, DatasetType, DatasetRef, Camera, CameraBuilder] = "camera",
        data_id: Optional[DataId] = None,
        **kwargs: Any,
    ) -> Camera:
        """
        Parameters
        ----------
        butler : `Butler`
            Butler client to read from.  If not initialized with the desired
            collection search path, the ``collections`` argument must be
            provided.
        camera : `str`, `DatasetType`, `DatasetRef`, `Camera`, \
                `CameraBuilder`, optional
            A `CameraBuilder`, a `Camera` to rebuild into one, or a butler
            dataset type or reference that can be used to load a `Camera`.
            Defaults to "camera" (as the nominal dataset type name).
        data_id : `dict` or `lsst.daf.butler.DataCoordinate`, optional
            Data ID values that identify the ``exposure`` and/or ``visit``
            dimension(s).  May include any extended keys supported by
            `Butler.get`, such as ``day_obs``.
        **kwargs
            Additional keyword arguments that are in the keys of
            `detector_calibrations` or `camera_calibations` will be interpreted
            as overrides of the calibrations to apply; see
            `CameraLoader.update` for details.  Other keyword arguments are
            interpreted as data ID key-value pairs and forwarded directly to
            `Butler.get`.

        Returns
        -------
        camera : `Camera`
            The updated camera.

        Examples
        --------
        Read the nominal camera and apply all default calibrations valid during
        the observation of exposure 500::

            camera = my_instrument.load_camera(butler, exposure=500)

        The same, but use day_obs/seq_num for exposure, and don't attempt to
        apply gains from PTC:

            camera = my_instrument.load_camera(
                butler,
                day_obs=20250129,
                seq_num=1261,
                ptc=False,
            )
        """
        update_camera_kwargs = {}
        for k, v in self.camera_calibrations.items():
            update_camera_kwargs[k] = kwargs.pop(k, v)
        update_detectors_kwargs = {}
        for k, v in self.detector_calibrations.items():
            update_detectors_kwargs = kwargs.pop(k, v)
        loader = CameraLoader(butler, camera, data_id=data_id, **kwargs)
        loader.update_camera(**update_camera_kwargs)
        loader.update_detectors(**update_detectors_kwargs)
        return loader.finish()

    @property
    def detector_calibrations(self) -> Mapping[str, str]:
        """A mapping of the calibration dataset types that should be used to
        update `Detector` objects for this instrument.
        """
        return {}

    @property
    def camera_calibrations(self) -> Mapping[str, str]:
        """A mapping of the calibration dataset types that should be used to
        update `Camera` objects for this instrument.

        This should be disjoint with `detector_calibrations`; those detector
        calibrations are also used to update camera objects, but are not part
        of this container.
        """
        return {}


def makeExposureRecordFromObsInfo(obsInfo: ObservationInfo, universe: DimensionUniverse) -> DimensionRecord:
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

    Raises
    ------
    LookupError
        Raised when ``dataId`` does not specify a valid data ID.
    """
    if collections is None:
        collections = butler.collections
    # Registry would do data ID expansion internally if we didn't do it first,
    # but we might want an expanded data ID ourselves later, so we do it here
    # to ensure it only happens once.
    # This will also catch problems with the data ID not having keys we need.
    try:
        dataId = butler.registry.expandDataId(dataId, graph=butler.registry.dimensions["exposure"].graph)
    except DataIdError as exc:
        raise LookupError(str(exc)) from exc
    try:
        cameraRef = butler.get("camera", dataId=dataId, collections=collections)
        return cameraRef, True
    except LookupError:
        pass
    # We know an instrument data ID is a value, but MyPy doesn't.
    instrument = Instrument.fromName(dataId["instrument"], butler.registry)  # type: ignore
    assert isinstance(instrument, Instrument)  # for mypy
    return instrument.getCamera(), False
