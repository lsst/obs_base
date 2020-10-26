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

__all__ = [
    "DefineVisitsConfig",
    "DefineVisitsTask",
    "GroupExposuresConfig",
    "GroupExposuresTask",
    "VisitDefinitionData",
]

from abc import ABCMeta, abstractmethod
from collections import defaultdict
import itertools
import dataclasses
from typing import Any, Dict, Iterable, List, Optional, Tuple
from multiprocessing import Pool

from lsst.daf.butler import (
    Butler,
    DataCoordinate,
    DataId,
    DimensionGraph,
    DimensionRecord,
    Timespan,
)

import lsst.geom
from lsst.geom import Box2D
from lsst.pex.config import Config, Field, makeRegistry, registerConfigurable
from lsst.afw.cameraGeom import FOCAL_PLANE, PIXELS
from lsst.pipe.base import Task
from lsst.sphgeom import ConvexPolygon, Region, UnitVector3d
from ._instrument import loadCamera, Instrument


@dataclasses.dataclass
class VisitDefinitionData:
    """Struct representing a group of exposures that will be used to define a
    visit.
    """

    instrument: str
    """Name of the instrument this visit will be associated with.
    """

    id: int
    """Integer ID of the visit.

    This must be unique across all visit systems for the instrument.
    """

    name: str
    """String name for the visit.

    This must be unique across all visit systems for the instrument.
    """

    exposures: List[DimensionRecord] = dataclasses.field(default_factory=list)
    """Dimension records for the exposures that are part of this visit.
    """


@dataclasses.dataclass
class _VisitRecords:
    """Struct containing the dimension records associated with a visit.
    """

    visit: DimensionRecord
    """Record for the 'visit' dimension itself.
    """

    visit_definition: List[DimensionRecord]
    """Records for 'visit_definition', which relates 'visit' to 'exposure'.
    """

    visit_detector_region: List[DimensionRecord]
    """Records for 'visit_detector_region', which associates the combination
    of a 'visit' and a 'detector' with a region on the sky.
    """


class GroupExposuresConfig(Config):
    pass


class GroupExposuresTask(Task, metaclass=ABCMeta):
    """Abstract base class for the subtask of `DefineVisitsTask` that is
    responsible for grouping exposures into visits.

    Subclasses should be registered with `GroupExposuresTask.registry` to
    enable use by `DefineVisitsTask`, and should generally correspond to a
    particular 'visit_system' dimension value.  They are also responsible for
    defining visit IDs and names that are unique across all visit systems in
    use by an instrument.

    Parameters
    ----------
    config : `GroupExposuresConfig`
        Configuration information.
    **kwargs
        Additional keyword arguments forwarded to the `Task` constructor.
    """
    def __init__(self, config: GroupExposuresConfig, **kwargs: Any):
        Task.__init__(self, config=config, **kwargs)

    ConfigClass = GroupExposuresConfig

    _DefaultName = "groupExposures"

    registry = makeRegistry(
        doc="Registry of algorithms for grouping exposures into visits.",
        configBaseType=GroupExposuresConfig,
    )

    @abstractmethod
    def group(self, exposures: List[DimensionRecord]) -> Iterable[VisitDefinitionData]:
        """Group the given exposures into visits.

        Parameters
        ----------
        exposures : `list` [ `DimensionRecord` ]
            DimensionRecords (for the 'exposure' dimension) describing the
            exposures to group.

        Returns
        -------
        visits : `Iterable` [ `VisitDefinitionData` ]
            Structs identifying the visits and the exposures associated with
            them.  This may be an iterator or a container.
        """
        raise NotImplementedError()

    @abstractmethod
    def getVisitSystem(self) -> Tuple[int, str]:
        """Return identifiers for the 'visit_system' dimension this
        algorithm implements.

        Returns
        -------
        id : `int`
            Integer ID for the visit system (given an instrument).
        name : `str`
            Unique string identifier for the visit system (given an
            instrument).
        """
        raise NotImplementedError()


class ComputeVisitRegionsConfig(Config):
    padding = Field(
        dtype=int,
        default=0,
        doc=("Pad raw image bounding boxes with specified number of pixels "
             "when calculating their (conservatively large) region on the "
             "sky."),
    )


class ComputeVisitRegionsTask(Task, metaclass=ABCMeta):
    """Abstract base class for the subtask of `DefineVisitsTask` that is
    responsible for extracting spatial regions for visits and visit+detector
    combinations.

    Subclasses should be registered with `ComputeVisitRegionsTask.registry` to
    enable use by `DefineVisitsTask`.

    Parameters
    ----------
    config : `ComputeVisitRegionsConfig`
        Configuration information.
    butler : `lsst.daf.butler.Butler`
        The butler to use.
    **kwargs
        Additional keyword arguments forwarded to the `Task` constructor.
    """
    def __init__(self, config: ComputeVisitRegionsConfig, *, butler: Butler, **kwargs: Any):
        Task.__init__(self, config=config, **kwargs)
        self.butler = butler
        self.instrumentMap = {}

    ConfigClass = ComputeVisitRegionsConfig

    _DefaultName = "computeVisitRegions"

    registry = makeRegistry(
        doc=("Registry of algorithms for computing on-sky regions for visits "
             "and visit+detector combinations."),
        configBaseType=ComputeVisitRegionsConfig,
    )

    def getInstrument(self, instrumentName) -> Instrument:
        """Retrieve an `~lsst.obs.base.Instrument` associated with this
        instrument name.

        Parameters
        ----------
        instrumentName : `str`
            The name of the instrument.

        Returns
        -------
        instrument : `~lsst.obs.base.Instrument`
            The associated instrument object.

        Notes
        -----
        The result is cached.
        """
        instrument = self.instrumentMap.get(instrumentName)
        if instrument is None:
            instrument = Instrument.fromName(instrumentName, self.butler.registry)
            self.instrumentMap[instrumentName] = instrument
        return instrument

    @abstractmethod
    def compute(self, visit: VisitDefinitionData, *, collections: Any = None
                ) -> Tuple[Region, Dict[int, Region]]:
        """Compute regions for the given visit and all detectors in that visit.

        Parameters
        ----------
        visit : `VisitDefinitionData`
            Struct describing the visit and the exposures associated with it.
        collections : Any, optional
            Collections to be searched for raws and camera geometry, overriding
            ``self.butler.collections``.
            Can be any of the types supported by the ``collections`` argument
            to butler construction.

        Returns
        -------
        visitRegion : `lsst.sphgeom.Region`
            Region for the full visit.
        visitDetectorRegions : `dict` [ `int`, `lsst.sphgeom.Region` ]
            Dictionary mapping detector ID to the region for that detector.
            Should include all detectors in the visit.
        """
        raise NotImplementedError()


class DefineVisitsConfig(Config):
    groupExposures = GroupExposuresTask.registry.makeField(
        doc="Algorithm for grouping exposures into visits.",
        default="one-to-one",
    )
    computeVisitRegions = ComputeVisitRegionsTask.registry.makeField(
        doc="Algorithm from computing visit and visit+detector regions.",
        default="single-raw-wcs",
    )
    ignoreNonScienceExposures = Field(
        doc=("If True, silently ignore input exposures that do not have "
             "observation_type=SCIENCE.  If False, raise an exception if one "
             "encountered."),
        dtype=bool,
        optional=False,
        default=True,
    )


class DefineVisitsTask(Task):
    """Driver Task for defining visits (and their spatial regions) in Gen3
    Butler repositories.

    Parameters
    ----------
    config : `DefineVisitsConfig`
        Configuration for the task.
    butler : `~lsst.daf.butler.Butler`
        Writeable butler instance.  Will be used to read `raw.wcs` and `camera`
        datasets and insert/sync dimension data.
    **kwargs
        Additional keyword arguments are forwarded to the `lsst.pipe.base.Task`
        constructor.

    Notes
    -----
    Each instance of `DefineVisitsTask` reads from / writes to the same Butler.
    Each invocation of `DefineVisitsTask.run` processes an independent group of
    exposures into one or more new vists, all belonging to the same visit
    system and instrument.

    The actual work of grouping exposures and computing regions is delegated
    to pluggable subtasks (`GroupExposuresTask` and `ComputeVisitRegionsTask`),
    respectively.  The defaults are to create one visit for every exposure,
    and to use exactly one (arbitrary) detector-level raw dataset's WCS along
    with camera geometry to compute regions for all detectors.  Other
    implementations can be created and configured for instruments for which
    these choices are unsuitable (e.g. because visits and exposures are not
    one-to-one, or because ``raw.wcs`` datasets for different detectors may not
    be consistent with camera geomery).

    It is not necessary in general to ingest all raws for an exposure before
    defining a visit that includes the exposure; this depends entirely on the
    `ComputeVisitRegionTask` subclass used.  For the default configuration,
    a single raw for each exposure is sufficient.

    Defining the same visit the same way multiple times (e.g. via multiple
    invocations of this task on the same exposures, with the same
    configuration) is safe, but it may be inefficient, as most of the work must
    be done before new visits can be compared to existing visits.
    """
    def __init__(self, config: Optional[DefineVisitsConfig] = None, *, butler: Butler, **kwargs: Any):
        config.validate()  # Not a CmdlineTask nor PipelineTask, so have to validate the config here.
        super().__init__(config, **kwargs)
        self.butler = butler
        self.universe = self.butler.registry.dimensions
        self.makeSubtask("groupExposures")
        self.makeSubtask("computeVisitRegions", butler=self.butler)

    def _reduce_kwargs(self):
        # Add extra parameters to pickle
        return dict(**super()._reduce_kwargs(), butler=self.butler)

    ConfigClass = DefineVisitsConfig

    _DefaultName = "defineVisits"

    def _buildVisitRecords(self, definition: VisitDefinitionData, *,
                           collections: Any = None) -> _VisitRecords:
        """Build the DimensionRecords associated with a visit.

        Parameters
        ----------
        definition : `VisitDefinition`
            Struct with identifiers for the visit and records for its
            constituent exposures.
        collections : Any, optional
            Collections to be searched for raws and camera geometry, overriding
            ``self.butler.collections``.
            Can be any of the types supported by the ``collections`` argument
            to butler construction.

        Results
        -------
        records : `_VisitRecords`
            Struct containing DimensionRecords for the visit, including
            associated dimension elements.
        """
        # Compute all regions.
        visitRegion, visitDetectorRegions = self.computeVisitRegions.compute(definition,
                                                                             collections=collections)
        # Aggregate other exposure quantities.
        timespan = Timespan(
            begin=_reduceOrNone(min, (e.timespan.begin for e in definition.exposures)),
            end=_reduceOrNone(max, (e.timespan.end for e in definition.exposures)),
        )
        exposure_time = _reduceOrNone(sum, (e.exposure_time for e in definition.exposures))
        physical_filter = _reduceOrNone(lambda a, b: a if a == b else None,
                                        (e.physical_filter for e in definition.exposures))
        target_name = _reduceOrNone(lambda a, b: a if a == b else None,
                                    (e.target_name for e in definition.exposures))
        science_program = _reduceOrNone(lambda a, b: a if a == b else None,
                                        (e.science_program for e in definition.exposures))
        observation_reason = _reduceOrNone(lambda a, b: a if a == b else None,
                                           (e.observation_reason for e in definition.exposures))
        if observation_reason is None:
            # Be explicit about there being multiple reasons
            observation_reason = "various"

        # Use the mean zenith angle as an approximation
        zenith_angle = _reduceOrNone(sum, (e.zenith_angle for e in definition.exposures))
        if zenith_angle is not None:
            zenith_angle /= len(definition.exposures)

        # Construct the actual DimensionRecords.
        return _VisitRecords(
            visit=self.universe["visit"].RecordClass(
                instrument=definition.instrument,
                id=definition.id,
                name=definition.name,
                physical_filter=physical_filter,
                target_name=target_name,
                science_program=science_program,
                observation_reason=observation_reason,
                zenith_angle=zenith_angle,
                visit_system=self.groupExposures.getVisitSystem()[0],
                exposure_time=exposure_time,
                timespan=timespan,
                region=visitRegion,
                # TODO: no seeing value in exposure dimension records, so we
                # can't set that here.  But there are many other columns that
                # both dimensions should probably have as well.
            ),
            visit_definition=[
                self.universe["visit_definition"].RecordClass(
                    instrument=definition.instrument,
                    visit=definition.id,
                    exposure=exposure.id,
                    visit_system=self.groupExposures.getVisitSystem()[0],
                )
                for exposure in definition.exposures
            ],
            visit_detector_region=[
                self.universe["visit_detector_region"].RecordClass(
                    instrument=definition.instrument,
                    visit=definition.id,
                    detector=detectorId,
                    region=detectorRegion,
                )
                for detectorId, detectorRegion in visitDetectorRegions.items()
            ]
        )

    def _expandExposureId(self, dataId: DataId) -> DataCoordinate:
        """Return the expanded version of an exposure ID.

        A private method to allow ID expansion in a pool without resorting
        to local callables.

        Parameters
        ----------
        dataId : `dict` or `DataCoordinate`
            Exposure-level data ID.

        Returns
        -------
        expanded : `DataCoordinate`
            A data ID that includes full metadata for all exposure dimensions.
        """
        dimensions = DimensionGraph(self.universe, names=["exposure"])
        return self.butler.registry.expandDataId(dataId, graph=dimensions)

    def _buildVisitRecordsSingle(self, args) -> _VisitRecords:
        """Build the DimensionRecords associated with a visit and collection.

        A wrapper for `_buildVisitRecords` to allow it to be run as part of
        a pool without resorting to local callables.

        Parameters
        ----------
        args : `tuple` [`VisitDefinition`, any]
            A tuple consisting of the ``definition`` and ``collections``
            arguments to `_buildVisitRecords`, in that order.

        Results
        -------
        records : `_VisitRecords`
            Struct containing DimensionRecords for the visit, including
            associated dimension elements.
        """
        return self._buildVisitRecords(args[0], collections=args[1])

    def run(self, dataIds: Iterable[DataId], *,
            pool: Optional[Pool] = None,
            processes: int = 1,
            collections: Optional[str] = None):
        """Add visit definitions to the registry for the given exposures.

        Parameters
        ----------
        dataIds : `Iterable` [ `dict` or `DataCoordinate` ]
            Exposure-level data IDs.  These must all correspond to the same
            instrument, and are expected to be on-sky science exposures.
        pool : `multiprocessing.Pool`, optional
            If not `None`, a process pool with which to parallelize some
            operations.
        processes : `int`, optional
            The number of processes to use.  Ignored if ``pool`` is not `None`.
        collections : Any, optional
            Collections to be searched for raws and camera geometry, overriding
            ``self.butler.collections``.
            Can be any of the types supported by the ``collections`` argument
            to butler construction.

        Raises
        ------
        lsst.daf.butler.registry.ConflictingDefinitionError
            Raised if a visit ID conflict is detected and the existing visit
            differs from the new one.
        """
        # Set up multiprocessing, if desired.
        if pool is None and processes > 1:
            pool = Pool(processes)
        mapFunc = map if pool is None else pool.imap_unordered
        # Normalize, expand, and deduplicate data IDs.
        self.log.info("Preprocessing data IDs.")
        dataIds = set(mapFunc(self._expandExposureId, dataIds))
        if not dataIds:
            raise RuntimeError("No exposures given.")
        # Extract exposure DimensionRecords, check that there's only one
        # instrument in play, and check for non-science exposures.
        exposures = []
        instruments = set()
        for dataId in dataIds:
            record = dataId.records["exposure"]
            if record.observation_type != "science":
                if self.config.ignoreNonScienceExposures:
                    continue
                else:
                    raise RuntimeError(f"Input exposure {dataId} has observation_type "
                                       f"{record.observation_type}, not 'science'.")
            instruments.add(dataId["instrument"])
            exposures.append(record)
        if not exposures:
            self.log.info("No science exposures found after filtering.")
            return
        if len(instruments) > 1:
            raise RuntimeError(
                f"All data IDs passed to DefineVisitsTask.run must be "
                f"from the same instrument; got {instruments}."
            )
        instrument, = instruments
        # Ensure the visit_system our grouping algorithm uses is in the
        # registry, if it wasn't already.
        visitSystemId, visitSystemName = self.groupExposures.getVisitSystem()
        self.log.info("Registering visit_system %d: %s.", visitSystemId, visitSystemName)
        self.butler.registry.syncDimensionData(
            "visit_system",
            {"instrument": instrument, "id": visitSystemId, "name": visitSystemName}
        )
        # Group exposures into visits, delegating to subtask.
        self.log.info("Grouping %d exposure(s) into visits.", len(exposures))
        definitions = list(self.groupExposures.group(exposures))
        # Compute regions and build DimensionRecords for each visit.
        # This is the only parallel step, but it _should_ be the most expensive
        # one (unless DB operations are slow).
        self.log.info("Computing regions and other metadata for %d visit(s).", len(definitions))
        allRecords = mapFunc(self._buildVisitRecordsSingle,
                             zip(definitions, itertools.repeat(collections)))
        # Iterate over visits and insert dimension data, one transaction per
        # visit.  If a visit already exists, we skip all other inserts.
        for visitRecords in allRecords:
            with self.butler.registry.transaction():
                if self.butler.registry.syncDimensionData("visit", visitRecords.visit):
                    self.butler.registry.insertDimensionData("visit_definition",
                                                             *visitRecords.visit_definition)
                    self.butler.registry.insertDimensionData("visit_detector_region",
                                                             *visitRecords.visit_detector_region)


def _reduceOrNone(func, iterable):
    """Apply a binary function to pairs of elements in an iterable until a
    single value is returned, but return `None` if any element is `None` or
    there are no elements.
    """
    r = None
    for v in iterable:
        if v is None:
            return None
        if r is None:
            r = v
        else:
            r = func(r, v)
    return r


class _GroupExposuresOneToOneConfig(GroupExposuresConfig):
    visitSystemId = Field(
        doc=("Integer ID of the visit_system implemented by this grouping "
             "algorithm."),
        dtype=int,
        default=0,
    )
    visitSystemName = Field(
        doc=("String name of the visit_system implemented by this grouping "
             "algorithm."),
        dtype=str,
        default="one-to-one",
    )


@registerConfigurable("one-to-one", GroupExposuresTask.registry)
class _GroupExposuresOneToOneTask(GroupExposuresTask, metaclass=ABCMeta):
    """An exposure grouping algorithm that simply defines one visit for each
    exposure, reusing the exposures identifiers for the visit.
    """

    ConfigClass = _GroupExposuresOneToOneConfig

    def group(self, exposures: List[DimensionRecord]) -> Iterable[VisitDefinitionData]:
        # Docstring inherited from GroupExposuresTask.
        for exposure in exposures:
            yield VisitDefinitionData(
                instrument=exposure.instrument,
                id=exposure.id,
                name=exposure.obs_id,
                exposures=[exposure],
            )

    def getVisitSystem(self) -> Tuple[int, str]:
        # Docstring inherited from GroupExposuresTask.
        return (self.config.visitSystemId, self.config.visitSystemName)


class _GroupExposuresByGroupMetadataConfig(GroupExposuresConfig):
    visitSystemId = Field(
        doc=("Integer ID of the visit_system implemented by this grouping "
             "algorithm."),
        dtype=int,
        default=1,
    )
    visitSystemName = Field(
        doc=("String name of the visit_system implemented by this grouping "
             "algorithm."),
        dtype=str,
        default="by-group-metadata",
    )


@registerConfigurable("by-group-metadata", GroupExposuresTask.registry)
class _GroupExposuresByGroupMetadataTask(GroupExposuresTask, metaclass=ABCMeta):
    """An exposure grouping algorithm that uses exposure.group_name and
    exposure.group_id.

    This algorithm _assumes_ exposure.group_id (generally populated from
    `astro_metadata_translator.ObservationInfo.visit_id`) is not just unique,
    but disjoint from all `ObservationInfo.exposure_id` values - if it isn't,
    it will be impossible to ever use both this grouping algorithm and the
    one-to-one algorithm for a particular camera in the same data repository.
    """

    ConfigClass = _GroupExposuresByGroupMetadataConfig

    def group(self, exposures: List[DimensionRecord]) -> Iterable[VisitDefinitionData]:
        # Docstring inherited from GroupExposuresTask.
        groups = defaultdict(list)
        for exposure in exposures:
            groups[exposure.group_name].append(exposure)
        for visitName, exposuresInGroup in groups.items():
            instrument = exposuresInGroup[0].instrument
            visitId = exposuresInGroup[0].group_id
            assert all(e.group_id == visitId for e in exposuresInGroup), \
                "Grouping by exposure.group_name does not yield consistent group IDs"
            yield VisitDefinitionData(instrument=instrument, id=visitId, name=visitName,
                                      exposures=exposuresInGroup)

    def getVisitSystem(self) -> Tuple[int, str]:
        # Docstring inherited from GroupExposuresTask.
        return (self.config.visitSystemId, self.config.visitSystemName)


class _ComputeVisitRegionsFromSingleRawWcsConfig(ComputeVisitRegionsConfig):
    mergeExposures = Field(
        doc=("If True, merge per-detector regions over all exposures in a "
             "visit (via convex hull) instead of using the first exposure and "
             "assuming its regions are valid for all others."),
        dtype=bool,
        default=False,
    )
    detectorId = Field(
        doc=("Load the WCS for the detector with this ID.  If None, use an "
             "arbitrary detector (the first found in a query of the data "
             "repository for each exposure (or all exposures, if "
             "mergeExposures is True)."),
        dtype=int,
        optional=True,
        default=None
    )
    requireVersionedCamera = Field(
        doc=("If True, raise LookupError if version camera geometry cannot be "
             "loaded for an exposure.  If False, use the nominal camera from "
             "the Instrument class instead."),
        dtype=bool,
        optional=False,
        default=False,
    )


@registerConfigurable("single-raw-wcs", ComputeVisitRegionsTask.registry)
class _ComputeVisitRegionsFromSingleRawWcsTask(ComputeVisitRegionsTask):
    """A visit region calculator that uses a single raw WCS and a camera to
    project the bounding boxes of all detectors onto the sky, relating
    different detectors by their positions in focal plane coordinates.

    Notes
    -----
    Most instruments should have their raw WCSs determined from a combination
    of boresight angle, rotator angle, and camera geometry, and hence this
    algorithm should produce stable results regardless of which detector the
    raw corresponds to.  If this is not the case (e.g. because a per-file FITS
    WCS is used instead), either the ID of the detector should be fixed (see
    the ``detectorId`` config parameter) or a different algorithm used.
    """

    ConfigClass = _ComputeVisitRegionsFromSingleRawWcsConfig

    def computeExposureBounds(self, exposure: DimensionRecord, *, collections: Any = None
                              ) -> Dict[int, List[UnitVector3d]]:
        """Compute the lists of unit vectors on the sphere that correspond to
        the sky positions of detector corners.

        Parameters
        ----------
        exposure : `DimensionRecord`
            Dimension record for the exposure.
        collections : Any, optional
            Collections to be searched for raws and camera geometry, overriding
            ``self.butler.collections``.
            Can be any of the types supported by the ``collections`` argument
            to butler construction.

        Returns
        -------
        bounds : `dict`
            Dictionary mapping detector ID to a list of unit vectors on the
            sphere representing that detector's corners projected onto the sky.
        """
        if collections is None:
            collections = self.butler.collections
        camera, versioned = loadCamera(self.butler, exposure.dataId, collections=collections)
        if not versioned and self.config.requireVersionedCamera:
            raise LookupError(f"No versioned camera found for exposure {exposure.dataId}.")

        # Derive WCS from boresight information -- if available in registry
        use_registry = True
        try:
            orientation = lsst.geom.Angle(exposure.sky_angle, lsst.geom.degrees)
            radec = lsst.geom.SpherePoint(lsst.geom.Angle(exposure.tracking_ra, lsst.geom.degrees),
                                          lsst.geom.Angle(exposure.tracking_dec, lsst.geom.degrees))
        except AttributeError:
            use_registry = False

        if use_registry:
            if self.config.detectorId is None:
                detectorId = next(camera.getIdIter())
            else:
                detectorId = self.config.detectorId
            wcsDetector = camera[detectorId]

            # Ask the raw formatter to create the relevant WCS
            # This allows flips to be taken into account
            instrument = self.getInstrument(exposure.instrument)
            rawFormatter = instrument.getRawFormatter({"detector": detectorId})
            wcs = rawFormatter.makeRawSkyWcsFromBoresight(radec, orientation, wcsDetector)

        else:
            if self.config.detectorId is None:
                wcsRefsIter = self.butler.registry.queryDatasets("raw.wcs", dataId=exposure.dataId,
                                                                 collections=collections)
                if not wcsRefsIter:
                    raise LookupError(f"No raw.wcs datasets found for data ID {exposure.dataId} "
                                      f"in collections {collections}.")
                wcsRef = next(iter(wcsRefsIter))
                wcsDetector = camera[wcsRef.dataId["detector"]]
                wcs = self.butler.getDirect(wcsRef)
            else:
                wcsDetector = camera[self.config.detectorId]
                wcs = self.butler.get("raw.wcs", dataId=exposure.dataId, detector=self.config.detectorId,
                                      collections=collections)
        fpToSky = wcsDetector.getTransform(FOCAL_PLANE, PIXELS).then(wcs.getTransform())
        bounds = {}
        for detector in camera:
            pixelsToSky = detector.getTransform(PIXELS, FOCAL_PLANE).then(fpToSky)
            pixCorners = Box2D(detector.getBBox().dilatedBy(self.config.padding)).getCorners()
            bounds[detector.getId()] = [
                skyCorner.getVector() for skyCorner in pixelsToSky.applyForward(pixCorners)
            ]
        return bounds

    def compute(self, visit: VisitDefinitionData, *, collections: Any = None
                ) -> Tuple[Region, Dict[int, Region]]:
        # Docstring inherited from ComputeVisitRegionsTask.
        if self.config.mergeExposures:
            detectorBounds = defaultdict(list)
            for exposure in visit.exposures:
                exposureDetectorBounds = self.computeExposureBounds(exposure, collections=collections)
                for detectorId, bounds in exposureDetectorBounds.items():
                    detectorBounds[detectorId].extend(bounds)
        else:
            detectorBounds = self.computeExposureBounds(visit.exposures[0], collections=collections)
        visitBounds = []
        detectorRegions = {}
        for detectorId, bounds in detectorBounds.items():
            detectorRegions[detectorId] = ConvexPolygon.convexHull(bounds)
            visitBounds.extend(bounds)
        return ConvexPolygon.convexHull(visitBounds), detectorRegions
