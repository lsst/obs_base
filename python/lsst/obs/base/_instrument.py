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

import logging
import os.path
import re
from abc import abstractmethod
from collections import defaultdict
from collections.abc import Sequence, Set
from functools import lru_cache
from typing import TYPE_CHECKING, Any, cast

import astropy.time
from lsst.afw.cameraGeom import Camera
from lsst.daf.butler import (
    Butler,
    CollectionType,
    DataCoordinate,
    DataId,
    DatasetType,
    DimensionRecord,
    DimensionUniverse,
    Timespan,
)
from lsst.daf.butler.registry import DataIdError
from lsst.pipe.base import Instrument as InstrumentBase
from lsst.utils import doImport, getPackageDir

from ._read_curated_calibs import CuratedCalibration, read_all

if TYPE_CHECKING:
    from astro_metadata_translator import MetadataTranslator, ObservationInfo
    from lsst.daf.butler import Registry

    from .filters import FilterDefinitionCollection

_LOG = logging.getLogger(__name__)

# To be a standard text curated calibration means that we use a
# standard definition for the corresponding DatasetType.
StandardCuratedCalibrationDatasetTypes = {
    "defects": {"dimensions": ("instrument", "detector"), "storageClass": "Defects"},
    "qe_curve": {"dimensions": ("instrument", "detector"), "storageClass": "QECurve"},
    "crosstalk": {"dimensions": ("instrument", "detector"), "storageClass": "CrosstalkCalib"},
    "linearizer": {"dimensions": ("instrument", "detector"), "storageClass": "Linearizer"},
    "bfk": {"dimensions": ("instrument", "detector"), "storageClass": "BrighterFatterKernel"},
    "transmission_optics": {"dimensions": ("instrument",), "storageClass": "TransmissionCurve"},
    "transmission_filter": {
        "dimensions": ("instrument", "physical_filter"),
        "storageClass": "TransmissionCurve",
    },
    "transmission_sensor": {"dimensions": ("instrument", "detector"), "storageClass": "TransmissionCurve"},
    "transmission_atmosphere": {"dimensions": ("instrument",), "storageClass": "TransmissionCurve"},
    "transmission_system": {
        "dimensions": ("instrument", "detector", "physical_filter"),
        "storageClass": "TransmissionCurve",
    },
}


class Instrument(InstrumentBase):
    """Rubin-specified base for instrument-specific logic for the Gen3 Butler.

    Parameters
    ----------
    collection_prefix : `str`, optional
        Prefix for collection names to use instead of the instrument's own
        name. This is primarily for use in simulated-data repositories, where
        the instrument name may not be necessary and/or sufficient to
        distinguish between collections.

    Notes
    -----
    Concrete instrument subclasses must have the same construction signature as
    the base class.
    """

    policyName: str | None = None
    """Instrument specific name to use when locating a policy or configuration
    file in the file system."""

    obsDataPackage: str | None = None
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

    translatorClass: MetadataTranslator | None = None
    """Class to use when extracting information from metadata. If `None`
    the metadata extraction system will determine the translator class itself.
    This class can also be used to calculate the observing day offset in some
    scenarios.
    """

    @property
    @abstractmethod
    def filterDefinitions(self) -> FilterDefinitionCollection:
        """`~lsst.obs.base.FilterDefinitionCollection`, defining the filters
        for this instrument.
        """
        raise NotImplementedError()

    def __init__(self, collection_prefix: str | None = None):
        super().__init__(collection_prefix=collection_prefix)

    @classmethod
    @lru_cache
    def getCuratedCalibrationNames(cls) -> frozenset[str]:
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
    @lru_cache
    def getObsDataPackageDir(cls) -> str | None:
        """Return the root of the obs data package that provides
        specializations for this instrument.

        Returns
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
        registry : `lsst.daf.butler.Registry`
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
        self, butler: Butler, collection: str | None = None, labels: Sequence[str] = ()
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
            them with the standard collection name delimiter.  If provided,
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
        and ``writeAdditionalCuratedCalibrations``.
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
        self, butler: Butler, collection: str | None = None, labels: Sequence[str] = ()
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
            them with the standard collection name delimiter.  If provided,
            these are inserted into the names of the `~CollectionType.RUN`
            collections that datasets are inserted directly into, as well the
            `~CollectionType.CALIBRATION` collection if it is generated
            automatically (i.e. if ``collection is None``).  Usually this is
            just the name of the ticket on which the calibration collection is
            being created.
        """
        return

    def writeCameraGeom(
        self, butler: Butler, collection: str | None = None, labels: Sequence[str] = ()
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
            them with the standard collection name delimiter.  If provided,
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
            "camera", ("instrument",), "Camera", isCalibration=True, universe=butler.dimensions
        )
        butler.registry.registerDatasetType(datasetType)
        camera = self.getCamera()
        ref = butler.put(camera, datasetType, {"instrument": self.getName()}, run=run)
        butler.registry.certify(collection, [ref], Timespan(begin=None, end=None))

    def writeStandardTextCuratedCalibrations(
        self, butler: Butler, collection: str | None = None, labels: Sequence[str] = ()
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
            them with the standard collection name delimiter.  If provided,
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
        runs: set[str] = set()
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
                universe=butler.dimensions,
                isCalibration=True,
                # MyPy should be able to figure out that the kwargs here have
                # the right types, but it can't.
                **definition,  # type: ignore
            )
            self._writeSpecificCuratedCalibrationDatasets(
                butler, datasetType, collection, runs=runs, labels=labels
            )

    @classmethod
    def _getSpecificCuratedCalibrationPath(cls, datasetTypeName: str) -> str | None:
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
        self, butler: Butler, datasetType: DatasetType, collection: str, runs: set[str], labels: Sequence[str]
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
        `~lsst.obs.base._read_curated_calibs.read_all` and provide standard
        metadata.
        """
        calibPath = self._getSpecificCuratedCalibrationPath(datasetType.name)
        if calibPath is None:
            return

        # Register the dataset type
        butler.registry.registerDatasetType(datasetType)
        _LOG.info("Processing %r curated calibration", datasetType.name)

        # The class to use to read these calibrations comes from the storage
        # class.
        calib_class: Any
        calib_class = datasetType.storageClass.pytype
        if not hasattr(calib_class, "readText"):
            # Let's try the default calib class.  All curated
            # calibrations should be subclasses of that, and the
            # parent can identify the correct one to use.
            calib_class = doImport("lsst.ip.isr.IsrCalib")

        calib_class = cast(type[CuratedCalibration], calib_class)

        # Read calibs, registering a new run for each CALIBDATE as needed.
        # We try to avoid registering runs multiple times as an optimization
        # by putting them in the ``runs`` set that was passed in.
        camera = self.getCamera()
        filters = set(self.filterDefinitions.physical_to_band.keys())
        calib_dimensions: list[str]
        if datasetType.name in StandardCuratedCalibrationDatasetTypes:
            calib_dimensions = list(StandardCuratedCalibrationDatasetTypes[datasetType.name]["dimensions"])
        else:
            # This should never trigger with real data, but will
            # trigger on the unit tests.
            _LOG.warning(
                "Unknown curated calibration type %s.  Attempting to use supplied definition.",
                datasetType.name,
            )
            calib_dimensions = list(datasetType.dimensions.names)

        calibsDict, calib_type = read_all(calibPath, camera, calib_class, calib_dimensions, filters)

        datasetRecords = []
        for path in calibsDict:
            times = sorted(calibsDict[path])
            calibs = [calibsDict[path][time] for time in times]
            atimes: list[astropy.time.Time | None] = [
                astropy.time.Time(t, format="datetime", scale="utc") for t in times
            ]
            atimes += [None]
            for calib, beginTime, endTime in zip(calibs, atimes[:-1], atimes[1:], strict=True):
                md = calib.getMetadata()
                run = self.makeCuratedCalibrationRunName(md["CALIBDATE"], *labels)
                if run not in runs:
                    butler.registry.registerRun(run)
                    runs.add(run)

                # DETECTOR and FILTER keywords in the calibration
                # metadata must exist if the calibration depends on
                # those dimensions.
                dimension_arguments = {}
                if "DETECTOR" in md:
                    dimension_arguments["detector"] = md["DETECTOR"]
                if "FILTER" in md:
                    dimension_arguments["physical_filter"] = md["FILTER"]

                dataId = DataCoordinate.standardize(
                    universe=butler.dimensions,
                    instrument=self.getName(),
                    **dimension_arguments,
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

    @classmethod
    def group_name_to_group_id(cls, group_name: str) -> int:
        """Translate the exposure group name to an integer.

        Parameters
        ----------
        group_name : `str`
            The name of the exposure group.

        Returns
        -------
        id : `int`
            The exposure group name in integer form. This integer might be
            used as an ID to uniquely identify the group in contexts where
            a string can not be used.

        Notes
        -----
        The default implementation removes all non numeric characters and casts
        to an integer.
        """
        cleaned = re.sub(r"\D", "", group_name)
        return int(cleaned)


def makeExposureRecordFromObsInfo(
    obsInfo: ObservationInfo, universe: DimensionUniverse, **kwargs: Any
) -> DimensionRecord:
    """Construct an exposure DimensionRecord from
    `astro_metadata_translator.ObservationInfo`.

    Parameters
    ----------
    obsInfo : `astro_metadata_translator.ObservationInfo`
        A `~astro_metadata_translator.ObservationInfo` object corresponding to
        the exposure.
    universe : `DimensionUniverse`
        Set of all known dimensions.
    **kwargs
        Additional field values for this record.

    Returns
    -------
    record : `DimensionRecord`
        A record containing exposure metadata, suitable for insertion into
        a `Registry`.
    """
    dimension = universe["exposure"]

    # Some registries support additional items.
    supported = {meta.name for meta in dimension.metadata}

    ra, dec, sky_angle, azimuth, zenith_angle = (None, None, None, None, None)
    if obsInfo.tracking_radec is not None:
        icrs = obsInfo.tracking_radec.icrs
        ra = icrs.ra.degree
        dec = icrs.dec.degree
        if obsInfo.boresight_rotation_coord == "sky":
            sky_angle = obsInfo.boresight_rotation_angle.degree
    if obsInfo.altaz_begin is not None:
        zenith_angle = obsInfo.altaz_begin.zen.degree
        azimuth = obsInfo.altaz_begin.az.degree

    extras: dict[str, Any] = {}
    for meta_key, info_key in (
        ("has_simulated", "has_simulated_content"),
        ("seq_start", "group_counter_start"),
        ("seq_end", "group_counter_end"),
        ("can_see_sky", "can_see_sky"),
    ):
        if meta_key in supported:
            extras[meta_key] = getattr(obsInfo, info_key)

    if (k := "azimuth") in supported:
        extras[k] = azimuth

    if "group" in dimension.implied:
        extras["group"] = obsInfo.exposure_group
    elif "group_name" in supported:
        extras["group_name"] = obsInfo.exposure_group
        extras["group_id"] = obsInfo.visit_id
    else:
        raise RuntimeError(f"Unable to determine where to put group metadata in exposure record: {supported}")

    # In some bad observations, the end time is before the begin time. We
    # can not let that be ingested as-is because it becomes an unbounded
    # timespan that will not work correctly with calibration lookups. Instead
    # force the end time to be the begin time.
    datetime_end = obsInfo.datetime_end
    if datetime_end < obsInfo.datetime_begin:
        datetime_end = obsInfo.datetime_begin
        _LOG.warning(
            "Exposure %s:%s has end time before begin time. Forcing it to use the begin time.",
            obsInfo.instrument,
            obsInfo.observation_id,
        )

    return dimension.RecordClass(
        instrument=obsInfo.instrument,
        id=obsInfo.exposure_id,
        obs_id=obsInfo.observation_id,
        datetime_begin=obsInfo.datetime_begin,
        datetime_end=datetime_end,
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
        **extras,
        **kwargs,
    )


def loadCamera(butler: Butler, dataId: DataId, *, collections: Any = None) -> tuple[Camera, bool]:
    """Attempt to load versioned camera geometry from a butler, but fall back
    to obtaining a nominal camera from the `Instrument` class if that fails.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
        Butler instance to attempt to query for and load a ``camera`` dataset
        from.
    dataId : `dict` or `~lsst.daf.butler.DataCoordinate`
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
        dataId = butler.registry.expandDataId(dataId, dimensions=butler.dimensions["exposure"].minimal_group)
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
