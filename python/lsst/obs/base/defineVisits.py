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
    "VisitSystem",
]

import cmath
import dataclasses
import enum
import math
import operator
from abc import ABCMeta, abstractmethod
from collections import defaultdict
from typing import (
    Any,
    Callable,
    ClassVar,
    Dict,
    FrozenSet,
    Iterable,
    List,
    Optional,
    Set,
    Tuple,
    Type,
    TypeVar,
)

import lsst.geom
from lsst.afw.cameraGeom import FOCAL_PLANE, PIXELS
from lsst.daf.butler import (
    Butler,
    DataCoordinate,
    DataId,
    DimensionGraph,
    DimensionRecord,
    Progress,
    Timespan,
)
from lsst.geom import Box2D
from lsst.pex.config import Config, Field, makeRegistry, registerConfigurable
from lsst.pipe.base import Instrument, Task
from lsst.sphgeom import ConvexPolygon, Region, UnitVector3d
from lsst.utils.introspection import get_full_type_name

from ._instrument import loadCamera


class VisitSystem(enum.Enum):
    """Enumeration used to label different visit systems."""

    ONE_TO_ONE = 0
    """Each exposure is assigned to its own visit."""

    BY_GROUP_METADATA = 1
    """Visit membership is defined by the value of the exposure.group_id."""

    BY_SEQ_START_END = 2
    """Visit membership is defined by the values of the ``exposure.day_obs``,
    ``exposure.seq_start``, and ``exposure.seq_end`` values.
    """

    @classmethod
    def all(cls) -> FrozenSet[VisitSystem]:
        """Return a `frozenset` containing all members."""
        return frozenset(cls.__members__.values())

    @classmethod
    def from_name(cls, external_name: str) -> VisitSystem:
        """Construct the enumeration from given name."""
        name = external_name.upper()
        name = name.replace("-", "_")
        try:
            return cls.__members__[name]
        except KeyError:
            raise KeyError(f"Visit system named '{external_name}' not known.") from None

    @classmethod
    def from_names(cls, names: Optional[Iterable[str]]) -> FrozenSet[VisitSystem]:
        """Return a `frozenset` of all the visit systems matching the supplied
        names.

        Parameters
        ----------
        names : iterable of `str`, or `None`
            Names of visit systems. Case insensitive. If `None` or empty, all
            the visit systems are returned.

        Returns
        -------
        systems : `frozenset` of `VisitSystem`
            The matching visit systems.
        """
        if not names:
            return cls.all()

        return frozenset({cls.from_name(name) for name in names})

    def __str__(self) -> str:
        name = self.name.lower()
        name = name.replace("_", "-")
        return name


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

    visit_systems: Set[VisitSystem]
    """All the visit systems associated with this visit."""

    exposures: List[DimensionRecord] = dataclasses.field(default_factory=list)
    """Dimension records for the exposures that are part of this visit.
    """


@dataclasses.dataclass
class _VisitRecords:
    """Struct containing the dimension records associated with a visit."""

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

    visit_system_membership: List[DimensionRecord]
    """Records relating visits to an associated visit system."""


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

    def getVisitSystems(self) -> Set[VisitSystem]:
        """Return identifiers for the 'visit_system' dimension this
        algorithm implements.

        Returns
        -------
        visit_systems : `Set` [`VisitSystem`]
            The visit systems used by this algorithm.
        """
        raise NotImplementedError()


class ComputeVisitRegionsConfig(Config):
    padding: Field[int] = Field(
        dtype=int,
        default=250,
        doc=(
            "Pad raw image bounding boxes with specified number of pixels "
            "when calculating their (conservatively large) region on the "
            "sky.  Note that the config value for pixelMargin of the "
            "reference object loaders in meas_algorithms should be <= "
            "the value set here."
        ),
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
        self.instrumentMap: Dict[str, Instrument] = {}

    ConfigClass = ComputeVisitRegionsConfig

    _DefaultName = "computeVisitRegions"

    registry = makeRegistry(
        doc=(
            "Registry of algorithms for computing on-sky regions for visits "
            "and visit+detector combinations."
        ),
        configBaseType=ComputeVisitRegionsConfig,
    )

    def getInstrument(self, instrumentName: str) -> Instrument:
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
    def compute(
        self, visit: VisitDefinitionData, *, collections: Any = None
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
        default="one-to-one-and-by-counter",
    )
    computeVisitRegions = ComputeVisitRegionsTask.registry.makeField(
        doc="Algorithm from computing visit and visit+detector regions.",
        default="single-raw-wcs",
    )
    ignoreNonScienceExposures: Field[bool] = Field(
        doc=(
            "If True, silently ignore input exposures that do not have "
            "observation_type=SCIENCE.  If False, raise an exception if one "
            "encountered."
        ),
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

    def __init__(self, config: DefineVisitsConfig, *, butler: Butler, **kwargs: Any):
        config.validate()  # Not a CmdlineTask nor PipelineTask, so have to validate the config here.
        super().__init__(config, **kwargs)
        self.butler = butler
        self.universe = self.butler.registry.dimensions
        self.progress = Progress("obs.base.DefineVisitsTask")
        self.makeSubtask("groupExposures")
        self.makeSubtask("computeVisitRegions", butler=self.butler)

    def _reduce_kwargs(self) -> dict:
        # Add extra parameters to pickle
        return dict(**super()._reduce_kwargs(), butler=self.butler)

    ConfigClass: ClassVar[Type[Config]] = DefineVisitsConfig

    _DefaultName: ClassVar[str] = "defineVisits"

    config: DefineVisitsConfig
    groupExposures: GroupExposuresTask
    computeVisitRegions: ComputeVisitRegionsTask

    def _buildVisitRecords(
        self, definition: VisitDefinitionData, *, collections: Any = None
    ) -> _VisitRecords:
        """Build the DimensionRecords associated with a visit.

        Parameters
        ----------
        definition : `VisitDefinitionData`
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
        dimension = self.universe["visit"]

        # Some registries support additional items.
        supported = {meta.name for meta in dimension.metadata}

        # Compute all regions.
        visitRegion, visitDetectorRegions = self.computeVisitRegions.compute(
            definition, collections=collections
        )
        # Aggregate other exposure quantities.
        timespan = Timespan(
            begin=_reduceOrNone(min, (e.timespan.begin for e in definition.exposures)),
            end=_reduceOrNone(max, (e.timespan.end for e in definition.exposures)),
        )
        exposure_time = _reduceOrNone(operator.add, (e.exposure_time for e in definition.exposures))
        physical_filter = _reduceOrNone(_value_if_equal, (e.physical_filter for e in definition.exposures))
        target_name = _reduceOrNone(_value_if_equal, (e.target_name for e in definition.exposures))
        science_program = _reduceOrNone(_value_if_equal, (e.science_program for e in definition.exposures))

        # observing day for a visit is defined by the earliest observation
        # of the visit
        observing_day = _reduceOrNone(min, (e.day_obs for e in definition.exposures))
        observation_reason = _reduceOrNone(
            _value_if_equal, (e.observation_reason for e in definition.exposures)
        )
        if observation_reason is None:
            # Be explicit about there being multiple reasons
            # MyPy can't really handle DimensionRecord fields as
            # DimensionRecord classes are dynamically defined; easiest to just
            # shush it when it complains.
            observation_reason = "various"  # type: ignore

        # Use the mean zenith angle as an approximation
        zenith_angle = _reduceOrNone(operator.add, (e.zenith_angle for e in definition.exposures))
        if zenith_angle is not None:
            zenith_angle /= len(definition.exposures)

        # New records that may not be supported.
        extras: Dict[str, Any] = {}
        if "seq_num" in supported:
            extras["seq_num"] = _reduceOrNone(min, (e.seq_num for e in definition.exposures))
        if "azimuth" in supported:
            # Must take into account 0/360 problem.
            extras["azimuth"] = _calc_mean_angle([e.azimuth for e in definition.exposures])

        # visit_system handling changed. This is the logic for visit/exposure
        # that has support for seq_start/seq_end.
        if "seq_num" in supported:
            # Map visit to exposure.
            visit_definition = [
                self.universe["visit_definition"].RecordClass(
                    instrument=definition.instrument,
                    visit=definition.id,
                    exposure=exposure.id,
                )
                for exposure in definition.exposures
            ]

            # Map visit to visit system.
            visit_system_membership = []
            for visit_system in self.groupExposures.getVisitSystems():
                if visit_system in definition.visit_systems:
                    record = self.universe["visit_system_membership"].RecordClass(
                        instrument=definition.instrument,
                        visit=definition.id,
                        visit_system=visit_system.value,
                    )
                    visit_system_membership.append(record)

        else:
            # The old approach can only handle one visit system at a time.
            # If we have been configured with multiple options, prefer the
            # one-to-one.
            visit_systems = self.groupExposures.getVisitSystems()
            if len(visit_systems) > 1:
                one_to_one = VisitSystem.from_name("one-to-one")
                if one_to_one not in visit_systems:
                    raise ValueError(
                        f"Multiple visit systems specified ({visit_systems}) for use with old"
                        " dimension universe but unable to find one-to-one."
                    )
                visit_system = one_to_one
            else:
                visit_system = visit_systems.pop()

            extras["visit_system"] = visit_system.value

            # The old visit_definition included visit system.
            visit_definition = [
                self.universe["visit_definition"].RecordClass(
                    instrument=definition.instrument,
                    visit=definition.id,
                    exposure=exposure.id,
                    visit_system=visit_system.value,
                )
                for exposure in definition.exposures
            ]

            # This concept does not exist in old schema.
            visit_system_membership = []

        # Construct the actual DimensionRecords.
        return _VisitRecords(
            visit=dimension.RecordClass(
                instrument=definition.instrument,
                id=definition.id,
                name=definition.name,
                physical_filter=physical_filter,
                target_name=target_name,
                science_program=science_program,
                observation_reason=observation_reason,
                day_obs=observing_day,
                zenith_angle=zenith_angle,
                exposure_time=exposure_time,
                timespan=timespan,
                region=visitRegion,
                # TODO: no seeing value in exposure dimension records, so we
                # can't set that here.  But there are many other columns that
                # both dimensions should probably have as well.
                **extras,
            ),
            visit_definition=visit_definition,
            visit_system_membership=visit_system_membership,
            visit_detector_region=[
                self.universe["visit_detector_region"].RecordClass(
                    instrument=definition.instrument,
                    visit=definition.id,
                    detector=detectorId,
                    region=detectorRegion,
                )
                for detectorId, detectorRegion in visitDetectorRegions.items()
            ],
        )

    def run(
        self,
        dataIds: Iterable[DataId],
        *,
        collections: Optional[str] = None,
        update_records: bool = False,
    ) -> None:
        """Add visit definitions to the registry for the given exposures.

        Parameters
        ----------
        dataIds : `Iterable` [ `dict` or `DataCoordinate` ]
            Exposure-level data IDs.  These must all correspond to the same
            instrument, and are expected to be on-sky science exposures.
        collections : Any, optional
            Collections to be searched for raws and camera geometry, overriding
            ``self.butler.collections``.
            Can be any of the types supported by the ``collections`` argument
            to butler construction.
        update_records : `bool`, optional
            If `True` (`False` is default), update existing visit records that
            conflict with the new ones instead of rejecting them (and when this
            occurs, update visit_detector_region as well).  THIS IS AN ADVANCED
            OPTION THAT SHOULD ONLY BE USED TO FIX REGIONS AND/OR METADATA THAT
            ARE KNOWN TO BE BAD, AND IT CANNOT BE USED TO REMOVE EXPOSURES OR
            DETECTORS FROM A VISIT.

        Raises
        ------
        lsst.daf.butler.registry.ConflictingDefinitionError
            Raised if a visit ID conflict is detected and the existing visit
            differs from the new one.
        """
        # Normalize, expand, and deduplicate data IDs.
        self.log.info("Preprocessing data IDs.")
        dimensions = DimensionGraph(self.universe, names=["exposure"])
        data_id_set: Set[DataCoordinate] = {
            self.butler.registry.expandDataId(d, graph=dimensions) for d in dataIds
        }
        if not data_id_set:
            raise RuntimeError("No exposures given.")
        # Extract exposure DimensionRecords, check that there's only one
        # instrument in play, and check for non-science exposures.
        exposures = []
        instruments = set()
        for dataId in data_id_set:
            record = dataId.records["exposure"]
            assert record is not None, "Guaranteed by expandDataIds call earlier."
            if record.tracking_ra is None or record.tracking_dec is None or record.sky_angle is None:
                if self.config.ignoreNonScienceExposures:
                    continue
                else:
                    raise RuntimeError(
                        f"Input exposure {dataId} has observation_type "
                        f"{record.observation_type}, but is not on sky."
                    )
            instruments.add(dataId["instrument"])
            exposures.append(record)
        if not exposures:
            self.log.info("No on-sky exposures found after filtering.")
            return
        if len(instruments) > 1:
            raise RuntimeError(
                f"All data IDs passed to DefineVisitsTask.run must be "
                f"from the same instrument; got {instruments}."
            )
        (instrument,) = instruments
        # Ensure the visit_system our grouping algorithm uses is in the
        # registry, if it wasn't already.
        visitSystems = self.groupExposures.getVisitSystems()
        for visitSystem in visitSystems:
            self.log.info("Registering visit_system %d: %s.", visitSystem.value, visitSystem)
            self.butler.registry.syncDimensionData(
                "visit_system", {"instrument": instrument, "id": visitSystem.value, "name": str(visitSystem)}
            )
        # Group exposures into visits, delegating to subtask.
        self.log.info("Grouping %d exposure(s) into visits.", len(exposures))
        definitions = list(self.groupExposures.group(exposures))
        # Iterate over visits, compute regions, and insert dimension data, one
        # transaction per visit.  If a visit already exists, we skip all other
        # inserts.
        self.log.info("Computing regions and other metadata for %d visit(s).", len(definitions))
        for visitDefinition in self.progress.wrap(
            definitions, total=len(definitions), desc="Computing regions and inserting visits"
        ):
            visitRecords = self._buildVisitRecords(visitDefinition, collections=collections)
            with self.butler.registry.transaction():
                inserted_or_updated = self.butler.registry.syncDimensionData(
                    "visit",
                    visitRecords.visit,
                    update=update_records,
                )
                if inserted_or_updated:
                    if inserted_or_updated is True:
                        # This is a new visit, not an update to an existing
                        # one, so insert visit definition.
                        # We don't allow visit definitions to change even when
                        # asked to update, because we'd have to delete the old
                        # visit_definitions first and also worry about what
                        # this does to datasets that already use the visit.
                        self.butler.registry.insertDimensionData(
                            "visit_definition", *visitRecords.visit_definition
                        )
                        if visitRecords.visit_system_membership:
                            self.butler.registry.insertDimensionData(
                                "visit_system_membership", *visitRecords.visit_system_membership
                            )
                    # [Re]Insert visit_detector_region records for both inserts
                    # and updates, because we do allow updating to affect the
                    # region calculations.
                    self.butler.registry.insertDimensionData(
                        "visit_detector_region", *visitRecords.visit_detector_region, replace=update_records
                    )


_T = TypeVar("_T")


def _reduceOrNone(func: Callable[[_T, _T], Optional[_T]], iterable: Iterable[Optional[_T]]) -> Optional[_T]:
    """Apply a binary function to pairs of elements in an iterable until a
    single value is returned, but return `None` if any element is `None` or
    there are no elements.
    """
    r: Optional[_T] = None
    for v in iterable:
        if v is None:
            return None
        if r is None:
            r = v
        else:
            r = func(r, v)
    return r


def _value_if_equal(a: _T, b: _T) -> Optional[_T]:
    """Return either argument if they are equal, or `None` if they are not."""
    return a if a == b else None


def _calc_mean_angle(angles: List[float]) -> float:
    """Calculate the mean angle, taking into account 0/360 wrapping.

    Parameters
    ----------
    angles : `list` [`float`]
        Angles to average together, in degrees.

    Returns
    -------
    average : `float`
        Average angle in degrees.
    """
    # Save on all the math if we only have one value.
    if len(angles) == 1:
        return angles[0]

    # Convert polar coordinates of unit circle to complex values.
    # Average the complex values.
    # Convert back to a phase angle.
    return math.degrees(cmath.phase(sum(cmath.rect(1.0, math.radians(d)) for d in angles) / len(angles)))


class _GroupExposuresOneToOneConfig(GroupExposuresConfig):
    visitSystemId: Field[int] = Field(
        doc="Integer ID of the visit_system implemented by this grouping algorithm.",
        dtype=int,
        default=0,
        deprecated="No longer used. Replaced by enum.",
    )
    visitSystemName: Field[str] = Field(
        doc="String name of the visit_system implemented by this grouping algorithm.",
        dtype=str,
        default="one-to-one",
        deprecated="No longer used. Replaced by enum.",
    )


@registerConfigurable("one-to-one", GroupExposuresTask.registry)
class _GroupExposuresOneToOneTask(GroupExposuresTask, metaclass=ABCMeta):
    """An exposure grouping algorithm that simply defines one visit for each
    exposure, reusing the exposures identifiers for the visit.
    """

    ConfigClass = _GroupExposuresOneToOneConfig

    def group(self, exposures: List[DimensionRecord]) -> Iterable[VisitDefinitionData]:
        # Docstring inherited from GroupExposuresTask.
        visit_systems = {VisitSystem.from_name("one-to-one")}
        for exposure in exposures:
            yield VisitDefinitionData(
                instrument=exposure.instrument,
                id=exposure.id,
                name=exposure.obs_id,
                exposures=[exposure],
                visit_systems=visit_systems,
            )

    def getVisitSystems(self) -> Set[VisitSystem]:
        # Docstring inherited from GroupExposuresTask.
        return set(VisitSystem.from_names(["one-to-one"]))


class _GroupExposuresByGroupMetadataConfig(GroupExposuresConfig):
    visitSystemId: Field[int] = Field(
        doc="Integer ID of the visit_system implemented by this grouping algorithm.",
        dtype=int,
        default=1,
        deprecated="No longer used. Replaced by enum.",
    )
    visitSystemName: Field[str] = Field(
        doc="String name of the visit_system implemented by this grouping algorithm.",
        dtype=str,
        default="by-group-metadata",
        deprecated="No longer used. Replaced by enum.",
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
        visit_systems = {VisitSystem.from_name("by-group-metadata")}
        groups = defaultdict(list)
        for exposure in exposures:
            groups[exposure.group_name].append(exposure)
        for visitName, exposuresInGroup in groups.items():
            instrument = exposuresInGroup[0].instrument
            visitId = exposuresInGroup[0].group_id
            assert all(
                e.group_id == visitId for e in exposuresInGroup
            ), "Grouping by exposure.group_name does not yield consistent group IDs"
            yield VisitDefinitionData(
                instrument=instrument,
                id=visitId,
                name=visitName,
                exposures=exposuresInGroup,
                visit_systems=visit_systems,
            )

    def getVisitSystems(self) -> Set[VisitSystem]:
        # Docstring inherited from GroupExposuresTask.
        return set(VisitSystem.from_names(["by-group-metadata"]))


class _GroupExposuresByCounterAndExposuresConfig(GroupExposuresConfig):
    visitSystemId: Field[int] = Field(
        doc="Integer ID of the visit_system implemented by this grouping algorithm.",
        dtype=int,
        default=2,
        deprecated="No longer used. Replaced by enum.",
    )
    visitSystemName: Field[str] = Field(
        doc="String name of the visit_system implemented by this grouping algorithm.",
        dtype=str,
        default="by-counter-and-exposures",
        deprecated="No longer used. Replaced by enum.",
    )


@registerConfigurable("one-to-one-and-by-counter", GroupExposuresTask.registry)
class _GroupExposuresByCounterAndExposuresTask(GroupExposuresTask, metaclass=ABCMeta):
    """An exposure grouping algorithm that uses the sequence start and
    sequence end metadata to create multi-exposure visits, but also
    creates one-to-one visits.

    This algorithm uses the exposure.seq_start and
    exposure.seq_end fields to collect related snaps.
    It also groups single exposures.
    """

    ConfigClass = _GroupExposuresByCounterAndExposuresConfig

    def group(self, exposures: List[DimensionRecord]) -> Iterable[VisitDefinitionData]:
        # Docstring inherited from GroupExposuresTask.
        system_one_to_one = VisitSystem.from_name("one-to-one")
        system_seq_start_end = VisitSystem.from_name("by-seq-start-end")

        groups = defaultdict(list)
        for exposure in exposures:
            groups[exposure.day_obs, exposure.seq_start, exposure.seq_end].append(exposure)
        for visit_key, exposures_in_group in groups.items():
            instrument = exposures_in_group[0].instrument

            # It is possible that the first exposure in a visit has not
            # been ingested. This can be determined and if that is the case
            # we can not reliably define the multi-exposure visit.
            skip_multi = False
            sorted_exposures = sorted(exposures_in_group, key=lambda e: e.seq_num)
            first = sorted_exposures.pop(0)
            if first.seq_num != first.seq_start:
                # Special case seq_num == 0 since that implies that the
                # instrument has no counters and therefore no multi-exposure
                # visits.
                if first.seq_num == 0:
                    self.log.warning(
                        "First exposure for visit %s is not present. Skipping the multi-snap definition.",
                        visit_key,
                    )
                skip_multi = True

            # Define the one-to-one visits.
            num_exposures = len(exposures_in_group)
            for exposure in exposures_in_group:
                # Default is to use the exposure ID and name unless
                # this is the first exposure in a multi-exposure visit.
                visit_name = exposure.obs_id
                visit_id = exposure.id
                visit_systems = {system_one_to_one}

                if num_exposures == 1:
                    # This is also a by-counter visit.
                    # It will use the same visit_name and visit_id.
                    visit_systems.add(system_seq_start_end)

                elif num_exposures > 1 and not skip_multi and exposure == first:
                    # This is the first legitimate exposure in a multi-exposure
                    # visit. It therefore needs a modified visit name and ID
                    # so it does not clash with the multi-exposure visit
                    # definition.
                    visit_name = f"{visit_name}_first"
                    visit_id = int(f"9{visit_id}")

                yield VisitDefinitionData(
                    instrument=instrument,
                    id=visit_id,
                    name=visit_name,
                    exposures=[exposure],
                    visit_systems=visit_systems,
                )

            # Multi-exposure visit.
            if not skip_multi and num_exposures > 1:
                # Define the visit using the first exposure
                visit_name = first.obs_id
                visit_id = first.id

                yield VisitDefinitionData(
                    instrument=instrument,
                    id=visit_id,
                    name=visit_name,
                    exposures=exposures_in_group,
                    visit_systems={system_seq_start_end},
                )

    def getVisitSystems(self) -> Set[VisitSystem]:
        # Docstring inherited from GroupExposuresTask.
        # Using a Config for this is difficult because what this grouping
        # algorithm is doing is using two visit systems.
        # One is using metadata (but not by-group) and the other is the
        # one-to-one. For now hard-code in class.
        return set(VisitSystem.from_names(["one-to-one", "by-seq-start-end"]))


class _ComputeVisitRegionsFromSingleRawWcsConfig(ComputeVisitRegionsConfig):
    mergeExposures: Field[bool] = Field(
        doc=(
            "If True, merge per-detector regions over all exposures in a "
            "visit (via convex hull) instead of using the first exposure and "
            "assuming its regions are valid for all others."
        ),
        dtype=bool,
        default=False,
    )
    detectorId: Field[Optional[int]] = Field(
        doc=(
            "Load the WCS for the detector with this ID.  If None, use an "
            "arbitrary detector (the first found in a query of the data "
            "repository for each exposure (or all exposures, if "
            "mergeExposures is True)."
        ),
        dtype=int,
        optional=True,
        default=None,
    )
    requireVersionedCamera: Field[bool] = Field(
        doc=(
            "If True, raise LookupError if version camera geometry cannot be "
            "loaded for an exposure.  If False, use the nominal camera from "
            "the Instrument class instead."
        ),
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
    config: _ComputeVisitRegionsFromSingleRawWcsConfig

    def computeExposureBounds(
        self, exposure: DimensionRecord, *, collections: Any = None
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
            radec = lsst.geom.SpherePoint(
                lsst.geom.Angle(exposure.tracking_ra, lsst.geom.degrees),
                lsst.geom.Angle(exposure.tracking_dec, lsst.geom.degrees),
            )
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

            try:
                wcs = rawFormatter.makeRawSkyWcsFromBoresight(radec, orientation, wcsDetector)  # type: ignore
            except AttributeError:
                raise TypeError(
                    f"Raw formatter is {get_full_type_name(rawFormatter)} but visit"
                    " definition requires it to support 'makeRawSkyWcsFromBoresight'"
                ) from None
        else:
            if self.config.detectorId is None:
                wcsRefsIter = self.butler.registry.queryDatasets(
                    "raw.wcs", dataId=exposure.dataId, collections=collections
                )
                if not wcsRefsIter:
                    raise LookupError(
                        f"No raw.wcs datasets found for data ID {exposure.dataId} "
                        f"in collections {collections}."
                    )
                wcsRef = next(iter(wcsRefsIter))
                wcsDetector = camera[wcsRef.dataId["detector"]]
                wcs = self.butler.getDirect(wcsRef)
            else:
                wcsDetector = camera[self.config.detectorId]
                wcs = self.butler.get(
                    "raw.wcs",
                    dataId=exposure.dataId,
                    detector=self.config.detectorId,
                    collections=collections,
                )
        fpToSky = wcsDetector.getTransform(FOCAL_PLANE, PIXELS).then(wcs.getTransform())
        bounds = {}
        for detector in camera:
            pixelsToSky = detector.getTransform(PIXELS, FOCAL_PLANE).then(fpToSky)
            pixCorners = Box2D(detector.getBBox().dilatedBy(self.config.padding)).getCorners()
            bounds[detector.getId()] = [
                skyCorner.getVector() for skyCorner in pixelsToSky.applyForward(pixCorners)
            ]
        return bounds

    def compute(
        self, visit: VisitDefinitionData, *, collections: Any = None
    ) -> Tuple[Region, Dict[int, Region]]:
        # Docstring inherited from ComputeVisitRegionsTask.
        if self.config.mergeExposures:
            detectorBounds: Dict[int, List[UnitVector3d]] = defaultdict(list)
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
