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
from collections.abc import Callable, Iterable, Sequence
from typing import Any, ClassVar, TypeVar, cast

import lsst.geom
from lsst.afw.cameraGeom import FOCAL_PLANE, PIXELS
from lsst.daf.butler import Butler, DataId, DimensionRecord, Progress, Timespan
from lsst.geom import Box2D
from lsst.pex.config import Config, Field, makeRegistry, registerConfigurable
from lsst.pipe.base import Struct, Task
from lsst.sphgeom import ConvexPolygon, Region, UnitVector3d
from lsst.utils.introspection import get_full_type_name

from ._instrument import Instrument, loadCamera


class VisitSystem(enum.Enum):
    """Enumeration used to label different visit systems."""

    ONE_TO_ONE = 0
    """Each exposure is assigned to its own visit."""

    BY_GROUP_METADATA = 1
    """Visit membership is defined by the value of the group dimension or, for
    older dimension universes, exposure.group_id."""

    BY_SEQ_START_END = 2
    """Visit membership is defined by the values of the ``exposure.day_obs``,
    ``exposure.seq_start``, and ``exposure.seq_end`` values.
    """

    @classmethod
    def all(cls) -> frozenset[VisitSystem]:
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
    def from_names(cls, names: Iterable[str] | None) -> frozenset[VisitSystem]:
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

    visit_systems: set[VisitSystem]
    """All the visit systems associated with this visit."""

    exposures: list[DimensionRecord] = dataclasses.field(default_factory=list)
    """Dimension records for the exposures that are part of this visit.
    """


@dataclasses.dataclass
class _VisitRecords:
    """Struct containing the dimension records associated with a visit."""

    visit: DimensionRecord
    """Record for the 'visit' dimension itself.
    """

    visit_definition: list[DimensionRecord]
    """Records for 'visit_definition', which relates 'visit' to 'exposure'.
    """

    visit_detector_region: list[DimensionRecord]
    """Records for 'visit_detector_region', which associates the combination
    of a 'visit' and a 'detector' with a region on the sky.
    """

    visit_system_membership: list[DimensionRecord]
    """Records relating visits to an associated visit system."""


class GroupExposuresConfig(Config):
    """Configure exposure grouping."""


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
        Additional keyword arguments forwarded to the `lsst.pipe.baseTask`
        constructor.
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
    def find_missing(
        self, exposures: list[DimensionRecord], registry: lsst.daf.butler.Registry
    ) -> list[DimensionRecord]:
        """Determine, if possible, which exposures might be missing.

        Parameters
        ----------
        exposures : `list` of `lsst.daf.butler.DimensionRecord`
            The exposure records to analyze.
        registry : `lsst.daf.butler.Registry`
            A butler registry that contains these exposure records.

        Returns
        -------
        missing : `list` of `lsst.daf.butler.DimensionRecord`
            Any exposure records present in registry that were related to
            the given exposures but were missing from that list and deemed
            to be relevant.

        Notes
        -----
        Only some grouping schemes are able to find missing exposures. It
        is acceptable to return an empty list.
        """
        raise NotImplementedError()

    @abstractmethod
    def group_exposures(self, exposures: list[DimensionRecord]) -> dict[Any, list[DimensionRecord]]:
        """Group the exposures in a way most natural for this visit definition.

        Parameters
        ----------
        exposures : `list` of `lsst.daf.butler.DimensionRecord`
            The exposure records to group.

        Returns
        -------
        groups : `dict` [Any, `list` [ `lsst.daf.butler.DimensionRecord` ] ]
            Groupings of exposure records. The key type is relevant to the
            specific visit definition and could be a string or a tuple.
        """
        raise NotImplementedError()

    @abstractmethod
    def group(
        self, exposures: list[DimensionRecord], instrument: Instrument
    ) -> Iterable[VisitDefinitionData]:
        """Group the given exposures into visits.

        Parameters
        ----------
        exposures : `list` [ `lsst.daf.butler.DimensionRecord` ]
            DimensionRecords (for the 'exposure' dimension) describing the
            exposures to group.
        instrument : `~lsst.obs.base.Instrument`
            Instrument specification that can be used to optionally support
            some visit ID definitions.

        Returns
        -------
        visits : `Iterable` [ `VisitDefinitionData` ]
            Structs identifying the visits and the exposures associated with
            them.  This may be an iterator or a container.
        """
        raise NotImplementedError()

    def getVisitSystems(self) -> set[VisitSystem]:
        """Return identifiers for the 'visit_system' dimension this
        algorithm implements.

        Returns
        -------
        visit_systems : `Set` [`VisitSystem`]
            The visit systems used by this algorithm.
        """
        raise NotImplementedError()


class ComputeVisitRegionsConfig(Config):
    """Configure visit region calculations."""

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
        self.instrumentMap: dict[str, Instrument] = {}

    ConfigClass = ComputeVisitRegionsConfig

    _DefaultName = "computeVisitRegions"

    registry = makeRegistry(
        doc="Registry of algorithms for computing on-sky regions for visits and visit+detector combinations.",
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
        self,
        visit: VisitDefinitionData,
        *,
        collections: Sequence[str] | str | None = None,
    ) -> tuple[Region, dict[int, Region]]:
        """Compute regions for the given visit and all detectors in that visit.

        Parameters
        ----------
        visit : `VisitDefinitionData`
            Struct describing the visit and the exposures associated with it.
        collections : `Sequence` [ `str` ] or `str` or `None`
            Collections to be searched for camera geometry, overriding
            ``self.butler.collections.defaults``. Can be any of the types
            supported by the ``collections`` argument to butler construction.

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
    """Configure visit definition."""

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
    updateObsCoreTable: Field[bool] = Field(
        doc=(
            "If True, update exposure regions in obscore table after visits "
            "are defined.  If False, do not update obscore table."
        ),
        dtype=bool,
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
        Writeable butler instance.  Will be used to read ``camera`` datasets
        and insert/sync dimension data.
    **kwargs
        Additional keyword arguments are forwarded to the `lsst.pipe.base.Task`
        constructor.

    Notes
    -----
    Each instance of `DefineVisitsTask` reads from / writes to the same Butler.
    Each invocation of `DefineVisitsTask.run` processes an independent group of
    exposures into one or more new visits, all belonging to the same visit
    system and instrument.

    The actual work of grouping exposures and computing regions is delegated to
    pluggable subtasks (`GroupExposuresTask` and `ComputeVisitRegionsTask`),
    respectively.  The defaults are to create one visit for every exposure, and
    to use exactly one (arbitrary) detector-level raw dataset's WCS along with
    camera geometry to compute regions for all detectors, but the raw WCS is
    recomputed from the ``exposure`` dimension record's rotation angle and
    boresight rather than by loading the ``raw.wcs`` dataset directly.  Other
    implementations can be created and configured for instruments for which
    these choices are unsuitable (e.g. because visits and exposures are not
    one-to-one, or because ``raw.wcs`` datasets for different detectors may not
    be consistent with camera geometry).

    Defining the same visit the same way multiple times (e.g. via multiple
    invocations of this task on the same exposures, with the same
    configuration) is safe, but it may be inefficient, as most of the work must
    be done before new visits can be compared to existing visits.
    """

    def __init__(self, config: DefineVisitsConfig, *, butler: Butler, **kwargs: Any):
        config.validate()  # Not a CmdlineTask nor PipelineTask, so have to validate the config here.
        super().__init__(config, **kwargs)
        self.butler = butler
        self.universe = self.butler.dimensions
        self.progress = Progress("obs.base.DefineVisitsTask")
        self.makeSubtask("groupExposures")
        self.makeSubtask("computeVisitRegions", butler=self.butler)

    def _reduce_kwargs(self) -> dict:
        # Add extra parameters to pickle
        return dict(**super()._reduce_kwargs(), butler=self.butler)

    ConfigClass: ClassVar[type[Config]] = DefineVisitsConfig

    _DefaultName: ClassVar[str] = "defineVisits"

    config: DefineVisitsConfig
    groupExposures: GroupExposuresTask
    computeVisitRegions: ComputeVisitRegionsTask

    def _buildVisitRecords(
        self, definition: VisitDefinitionData, *, collections: Sequence[str] | str | None = None
    ) -> _VisitRecords:
        """Build the DimensionRecords associated with a visit.

        Parameters
        ----------
        definition : `VisitDefinitionData`
            Struct with identifiers for the visit and records for its
            constituent exposures.
        collections : `Sequence` [ `str` ] or `str` or `None`
            Collections to be searched for camera geometry, overriding
            ``self.butler.collections.defaults``. Can be any of the types
            supported by the ``collections`` argument to butler construction.

        Returns
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
            observation_reason = "various"

        # Use the mean zenith angle as an approximation
        zenith_angle = _reduceOrNone(operator.add, (e.zenith_angle for e in definition.exposures))
        if zenith_angle is not None:
            zenith_angle /= len(definition.exposures)

        # New records that may not be supported.
        extras: dict[str, Any] = {}
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
        dataIds_or_records: Iterable[DataId | DimensionRecord],
        *,
        collections: Sequence[str] | str | None = None,
        update_records: bool = False,
        incremental: bool = False,
    ) -> Struct:
        """Add visit definitions to the registry for the given exposures.

        Parameters
        ----------
        dataIds_or_records : `Iterable` [ `dict` or \
              `~lsst.daf.butler.DataCoordinate` or \
              `~lsst.daf.butler.DimensionRecord` ]
            Exposure-level data IDs or explicit exposure records.  These must
            all correspond to the same instrument, and are expected to be
            on-sky science exposures.
        collections : `Sequence` [ `str` ] or `str` or `None`
            Collections to be searched for camera geometry, overriding
            ``self.butler.collections.defaults``. Can be any of the types
            supported by the ``collections`` argument to butler construction.
        update_records : `bool`, optional
            If `True` (`False` is default), update existing ``visit`` records
            and ``visit_detector_region`` records.  THIS IS AN ADVANCED OPTION
            THAT SHOULD ONLY BE USED TO FIX REGIONS AND/OR METADATA THAT ARE
            KNOWN TO BE BAD, AND IT CANNOT BE USED TO REMOVE EXPOSURES OR
            DETECTORS FROM A VISIT.
        incremental : `bool`, optional
            If `True` indicate that exposures are being ingested incrementally
            and visit definition will be run on partial visits.  This will
            allow the ``visit`` record to be updated if it already exists, but
            (unlike ``update_records=True``) it will only update the
            ``visit_detector_region`` records if the ``visit`` record's region
            changes. If there is any risk that files are being ingested
            incrementally it is critical that this parameter is set to `True`
            and not to rely on ``update_records``.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Structure with the following attributes (all `int`):

            - n_visits: total number of visits defined
            - n_skipped: number of visits that were already present left alone
            - n_new: number of new visits inserted
            - n_fully_updated: number of existing visits fully updated
            - n_partially_updated: number of visits with non-geometry updates.

        Raises
        ------
        lsst.daf.butler.registry.ConflictingDefinitionError
            Raised if a visit ID conflict is detected and the existing visit
            differs from the new one.
        """
        # Normalize, expand, and deduplicate data IDs.
        self.log.info("Preprocessing data IDs.")
        dimensions = self.universe.conform(["exposure"])

        exposures: list[DimensionRecord] = []
        instruments: set[str] = set()
        instrument_cls_name: str | None = None
        instrument_record: DimensionRecord | None = None

        # Go through the supplied dataset extracting records.
        # Check that only a single instrument is being used.
        for external in dataIds_or_records:
            if isinstance(external, DimensionRecord):
                record = external
                if str(record.definition) != "exposure":
                    raise ValueError(f"Can only define visits from exposure records, not {record}.")
            else:
                data_id = self.butler.registry.expandDataId(external, dimensions=dimensions)
                exp_record = data_id.records["exposure"]
                assert exp_record is not None, "Guaranteed by expandDataIds call earlier."
                record = exp_record
                instrument_record = data_id.records["instrument"]

            # LSSTCam data can assign ra/dec to flats, and dome-closed
            # engineering tests. Do not assign a visit if we know that
            # can_see_sky is False. Treat None as True for this test.
            can_see_sky = getattr(record, "can_see_sky", True)
            if (
                record.tracking_ra is None
                or record.tracking_dec is None
                or record.sky_angle is None
                or can_see_sky is False
            ):
                if self.config.ignoreNonScienceExposures:
                    continue
                else:
                    raise RuntimeError(
                        f"Input exposure {external} has observation_type "
                        f"{record.observation_type}, but is not on sky."
                    )
            instrument_name = record.instrument
            instruments.add(instrument_name)
            exposures.append(record)
        if not exposures:
            self.log.info("No on-sky exposures found after filtering.")
            return Struct(n_visits=0, n_skipped=0, n_new=0, n_partially_updated=0, n_fully_updated=0)
        if len(instruments) > 1:
            raise RuntimeError(
                "All data IDs passed to DefineVisitsTask.run must be "
                f"from the same instrument; got {instruments}."
            )
        (instrument,) = instruments

        # Might need the instrument class for later depending on universe
        # and grouping scheme.
        if instrument_cls_name is None:
            if instrument_record is None:
                # We were given a DimensionRecord instead of a DataCoordinate.

                instrument_records = self.butler.query_dimension_records(
                    "instrument", instrument=instrument, limit=1
                )
                if len(instrument_records) != 1:
                    raise RuntimeError(
                        f"Instrument {instrument} found in dimension record but unknown to butler."
                    )
                instrument_record = instrument_records[0]
            instrument_cls_name = instrument_record.class_name
        assert instrument_cls_name is not None, "Instrument must be defined by this point"
        instrument_helper = Instrument.from_string(instrument_cls_name)

        # Ensure the visit_system our grouping algorithm uses is in the
        # registry, if it wasn't already.
        visitSystems = self.groupExposures.getVisitSystems()
        for visitSystem in visitSystems:
            self.log.info("Registering visit_system %d: %s.", visitSystem.value, visitSystem)
            self.butler.registry.syncDimensionData(
                "visit_system",
                {"instrument": instrument, "id": visitSystem.value, "name": str(visitSystem)},
            )

        # In true incremental we will be given the second snap on its
        # own on the assumption that the previous snap was already handled.
        # For correct grouping we need access to the other exposures in the
        # visit.
        if incremental:
            exposures.extend(self.groupExposures.find_missing(exposures, self.butler.registry))

        # Group exposures into visits, delegating to subtask.
        self.log.info("Grouping %d exposure(s) into visits.", len(exposures))
        definitions = list(self.groupExposures.group(exposures, instrument_helper))
        # Iterate over visits, compute regions, and insert dimension data, one
        # transaction per visit.  If a visit already exists, we skip all other
        # inserts.
        self.log.info("Computing regions and other metadata for %d visit(s).", len(definitions))
        n_skipped: int = 0
        n_new: int = 0
        n_fully_updated: int = 0
        n_partially_updated: int = 0
        for visitDefinition in self.progress.wrap(
            definitions, total=len(definitions), desc="Computing regions and inserting visits"
        ):
            visitRecords = self._buildVisitRecords(visitDefinition, collections=collections)
            with self.butler.registry.transaction():
                inserted_or_updated = self.butler.registry.syncDimensionData(
                    "visit",
                    visitRecords.visit,
                    update=(update_records or incremental),
                )
                if inserted_or_updated or update_records:
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
                    elif incremental and len(visitRecords.visit_definition) > 1:
                        # The visit record was modified. This could happen
                        # if a multi-snap visit was redefined with an
                        # additional snap so play it safe and allow for the
                        # visit definition to be updated. We use update=False
                        # here since there should not be any rows updated,
                        # just additional rows added. update=True does not work
                        # correctly with multiple records. In incremental mode
                        # we assume that the caller wants the visit definition
                        # to be updated and has no worries about provenance
                        # with the previous definition.
                        for definition in visitRecords.visit_definition:
                            self.butler.registry.syncDimensionData("visit_definition", definition)
                    if inserted_or_updated is True:
                        # Insert visit-detector regions if the visit is new.
                        self.butler.registry.insertDimensionData(
                            "visit_detector_region",
                            *visitRecords.visit_detector_region,
                            replace=False,
                        )
                        self.log.verbose(
                            "Inserted %s visit_detector_region records for new visit %s.",
                            len(visitRecords.visit_detector_region),
                            visitRecords.visit.id,
                        )
                        n_new += 1
                    # Cast below is because MyPy can't determine that
                    # inserted_or_updated can only be False if update_records
                    # is True.
                    elif update_records or "region" in cast(dict, inserted_or_updated):
                        # Replace visit-detector regions if we were told to
                        # update records explicitly, or if the visit region
                        # changed in an incremental=True update.
                        self.butler.registry.insertDimensionData(
                            "visit_detector_region",
                            *visitRecords.visit_detector_region,
                            replace=True,
                        )
                        self.log.verbose(
                            "Re-inserted %s visit_detector_region records for updated visit %s.",
                            len(visitRecords.visit_detector_region),
                            visitRecords.visit.id,
                        )
                        n_fully_updated += 1
                    else:
                        self.log.verbose(
                            "Updated visit %s without modifying visit_detector_region records.",
                            visitRecords.visit.id,
                        )
                        n_partially_updated += 1

                    # Update obscore exposure records with region information
                    # from corresponding visits.
                    if self.config.updateObsCoreTable:
                        if obscore_manager := self.butler.registry.obsCoreTableManager:
                            obscore_updates: list[tuple[int, int, Region]] = []
                            exposure_ids = [rec.exposure for rec in visitRecords.visit_definition]
                            for record in visitRecords.visit_detector_region:
                                obscore_updates += [
                                    (exposure, record.detector, record.region) for exposure in exposure_ids
                                ]
                            if obscore_updates:
                                obscore_manager.update_exposure_regions(instrument, obscore_updates)
                else:
                    self.log.verbose("Skipped already-existing visit %s.", visitRecords.visit.id)
                    n_skipped += 1
        self.log.info(
            "Finished writing database records for %d visit(s): %s left unchanged, %s new, "
            "%s updated with new detector regions, %s updated without new detector regions.",
            len(definitions),
            n_skipped,
            n_new,
            n_fully_updated,
            n_partially_updated,
        )
        return Struct(
            n_visits=len(definitions),
            n_skipped=n_skipped,
            n_new=n_new,
            n_fully_updated=n_fully_updated,
            n_partially_updated=n_partially_updated,
        )


_T = TypeVar("_T")


def _reduceOrNone(func: Callable[[_T, _T], _T | None], iterable: Iterable[_T | None]) -> _T | None:
    """Apply a binary function to pairs of elements in an iterable until a
    single value is returned, but return `None` if any element is `None` or
    there are no elements.
    """
    r: _T | None = None
    for v in iterable:
        if v is None:
            return None
        if r is None:
            r = v
        else:
            r = func(r, v)
    return r


def _value_if_equal(a: _T, b: _T) -> _T | None:
    """Return either argument if they are equal, or `None` if they are not."""
    return a if a == b else None


def _calc_mean_angle(angles: list[float]) -> float:
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

    def find_missing(
        self, exposures: list[DimensionRecord], registry: lsst.daf.butler.Registry
    ) -> list[DimensionRecord]:
        # By definition no exposures can be missing.
        return []

    def group_exposures(self, exposures: list[DimensionRecord]) -> dict[Any, list[DimensionRecord]]:
        # No grouping.
        return {exposure.id: [exposure] for exposure in exposures}

    def group(
        self, exposures: list[DimensionRecord], instrument: Instrument
    ) -> Iterable[VisitDefinitionData]:
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

    def getVisitSystems(self) -> set[VisitSystem]:
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
    """An exposure grouping algorithm that uses the exposure group.

    This algorithm uses the ``group`` dimension for modern universes and the
    ``exposure.group_id`` for older universes.

    This algorithm *assumes* group ID (generally populated from
    `astro_metadata_translator.ObservationInfo.visit_id`) is not just unique,
    but disjoint from all `ObservationInfo.exposure_id` values - if it isn't,
    it will be impossible to ever use both this grouping algorithm and the
    one-to-one algorithm for a particular camera in the same data repository.
    """

    ConfigClass = _GroupExposuresByGroupMetadataConfig

    def find_missing(
        self, exposures: list[DimensionRecord], registry: lsst.daf.butler.Registry
    ) -> list[DimensionRecord]:
        groups = self.group_exposures(exposures)
        # Determine which group implementation we are using.
        if "group" in registry.dimensions["exposure"].implied:
            group_key = "group"
        else:
            group_key = "group_name"
        missing_exposures: list[DimensionRecord] = []
        for exposures_in_group in groups.values():
            # We can not tell how many exposures are expected to be in the
            # visit so we have to query every time.
            first = exposures_in_group[0]
            records = set(
                registry.queryDimensionRecords(
                    "exposure",
                    where=f"exposure.{group_key} = groupnam",
                    bind={"groupnam": getattr(first, group_key)},
                    instrument=first.instrument,
                )
            )
            records.difference_update(set(exposures_in_group))
            missing_exposures.extend(list(records))
        return missing_exposures

    def group_exposures(self, exposures: list[DimensionRecord]) -> dict[Any, list[DimensionRecord]]:
        groups = defaultdict(list)
        group_key = "group"
        if exposures and hasattr(exposures[0], "group_name"):
            group_key = "group_name"
        for exposure in exposures:
            groups[getattr(exposure, group_key)].append(exposure)
        return groups

    def group(
        self, exposures: list[DimensionRecord], instrument: Instrument
    ) -> Iterable[VisitDefinitionData]:
        # Docstring inherited from GroupExposuresTask.
        visit_systems = {VisitSystem.from_name("by-group-metadata")}
        groups = self.group_exposures(exposures)
        has_group_dimension: bool | None = None
        for visitName, exposuresInGroup in groups.items():
            instrument_name = exposuresInGroup[0].instrument
            assert instrument_name == instrument.getName(), "Inconsistency in instrument name"
            visit_ids: set[int] = set()
            if has_group_dimension is None:
                has_group_dimension = hasattr(exposuresInGroup[0], "group")
            if has_group_dimension:
                visit_ids = {instrument.group_name_to_group_id(e.group) for e in exposuresInGroup}
            else:
                visit_ids = {e.group_id for e in exposuresInGroup}
            assert len(visit_ids) == 1, "Grouping by exposure group does not yield consistent group IDs"
            yield VisitDefinitionData(
                instrument=instrument_name,
                id=visit_ids.pop(),
                name=visitName,
                exposures=exposuresInGroup,
                visit_systems=visit_systems,
            )

    def getVisitSystems(self) -> set[VisitSystem]:
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

    def find_missing(
        self, exposures: list[DimensionRecord], registry: lsst.daf.butler.Registry
    ) -> list[DimensionRecord]:
        """Analyze the exposures and return relevant exposures known to
        registry.
        """
        groups = self.group_exposures(exposures)
        missing_exposures: list[DimensionRecord] = []
        for exposures_in_group in groups.values():
            sorted_exposures = sorted(exposures_in_group, key=lambda e: e.seq_num)
            first = sorted_exposures[0]

            # Only need to look for the seq_nums that we don't already have.
            seq_nums = set(range(first.seq_start, first.seq_end + 1))
            seq_nums.difference_update({exp.seq_num for exp in sorted_exposures})

            if seq_nums:
                # Missing something. Check registry.
                records = list(
                    registry.queryDimensionRecords(
                        "exposure",
                        where="exposure.seq_start = seq_start AND exposure.seq_end = seq_end AND "
                        "exposure.seq_num IN (seq_nums)",
                        bind={"seq_start": first.seq_start, "seq_end": first.seq_end, "seq_nums": seq_nums},
                        instrument=first.instrument,
                    )
                )
                missing_exposures.extend(records)

        return missing_exposures

    def group_exposures(self, exposures: list[DimensionRecord]) -> dict[Any, list[DimensionRecord]]:
        groups = defaultdict(list)
        for exposure in exposures:
            groups[exposure.day_obs, exposure.seq_start, exposure.seq_end].append(exposure)
        return groups

    def group(
        self, exposures: list[DimensionRecord], instrument: Instrument
    ) -> Iterable[VisitDefinitionData]:
        # Docstring inherited from GroupExposuresTask.
        system_one_to_one = VisitSystem.from_name("one-to-one")
        system_seq_start_end = VisitSystem.from_name("by-seq-start-end")

        groups = self.group_exposures(exposures)
        for visit_key, exposures_in_group in groups.items():
            instrument_name = exposures_in_group[0].instrument

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
                if first.seq_num != 0:
                    self.log.warning(
                        "First exposure for visit %s is not present. Skipping the multi-snap definition.",
                        visit_key,
                    )
                skip_multi = True

            multi_exposure = False
            if first.seq_start != first.seq_end:
                # This is a multi-exposure visit regardless of the number
                # of exposures present.
                multi_exposure = True

            # Define the one-to-one visits.
            for exposure in exposures_in_group:
                # Default is to use the exposure ID and name unless
                # this is the first exposure in a multi-exposure visit.
                visit_name = exposure.obs_id
                visit_id = exposure.id
                visit_systems = {system_one_to_one}

                if not multi_exposure:
                    # This is also a by-counter visit.
                    # It will use the same visit_name and visit_id.
                    visit_systems.add(system_seq_start_end)

                elif not skip_multi and exposure == first:
                    # This is the first legitimate exposure in a multi-exposure
                    # visit. It therefore needs a modified visit name and ID
                    # so it does not clash with the multi-exposure visit
                    # definition.
                    visit_name = f"{visit_name}_first"
                    visit_id = int(f"9{visit_id}")

                yield VisitDefinitionData(
                    instrument=instrument_name,
                    id=visit_id,
                    name=visit_name,
                    exposures=[exposure],
                    visit_systems=visit_systems,
                )

            # Multi-exposure visit.
            if not skip_multi and multi_exposure:
                # Define the visit using the first exposure
                visit_name = first.obs_id
                visit_id = first.id

                yield VisitDefinitionData(
                    instrument=instrument_name,
                    id=visit_id,
                    name=visit_name,
                    exposures=exposures_in_group,
                    visit_systems={system_seq_start_end},
                )

    def getVisitSystems(self) -> set[VisitSystem]:
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
    detectorId: Field[int | None] = Field(
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
    """A visit region calculator that uses a single raw WCS (recomputed from
    the ``exposure`` dimension record) and a camera to project the bounding
    boxes of all detectors onto the sky, relating different detectors by their
    positions in focal plane coordinates.

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
    ) -> dict[int, list[UnitVector3d]]:
        """Compute the lists of unit vectors on the sphere that correspond to
        the sky positions of detector corners.

        Parameters
        ----------
        exposure : `DimensionRecord`
            Dimension record for the exposure.
        collections : Any, optional
            Collections to be searched for raws and camera geometry, overriding
            ``self.butler.collections.defaults``.
            Can be any of the types supported by the ``collections`` argument
            to butler construction.

        Returns
        -------
        bounds : `dict`
            Dictionary mapping detector ID to a list of unit vectors on the
            sphere representing that detector's corners projected onto the sky.
        """
        if collections is None:
            collections = list(self.butler.collections.defaults)
        camera, versioned = loadCamera(self.butler, exposure.dataId, collections=collections)
        if not versioned and self.config.requireVersionedCamera:
            raise LookupError(f"No versioned camera found for exposure {exposure.dataId}.")

        orientation = lsst.geom.Angle(exposure.sky_angle, lsst.geom.degrees)
        radec = lsst.geom.SpherePoint(
            lsst.geom.Angle(exposure.tracking_ra, lsst.geom.degrees),
            lsst.geom.Angle(exposure.tracking_dec, lsst.geom.degrees),
        )

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
    ) -> tuple[Region, dict[int, Region]]:
        # Docstring inherited from ComputeVisitRegionsTask.
        if self.config.mergeExposures:
            detectorBounds: dict[int, list[UnitVector3d]] = defaultdict(list)
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
