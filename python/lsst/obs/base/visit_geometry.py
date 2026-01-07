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

__all__ = ("VisitGeometry",)

import itertools
from typing import Any

import pydantic

from lsst.daf.butler import Butler, DimensionRecord
from lsst.daf.butler.pydantic_utils import SerializableRegion
from lsst.utils.logging import getLogger

_LOG = getLogger(__name__)


class VisitGeometry(pydantic.BaseModel):
    """A serializable struct that holds geometry information for a visit.

    This is intended to be used as a butler output dataset intermediary between
    a tasks that fit coordinate mappings from stars to reference catalogs and
    tool that update butler dimension record regions.
    """

    boresight_ra: float
    """Re-fit boresight right ascension in degrees."""

    boresight_dec: float
    """Re-fit boresight declination in degrees."""

    orientation: float
    """Re-fit rotation angle in degrees."""

    visit_region: SerializableRegion
    """Updated region for the visit."""

    detector_regions: dict[int, SerializableRegion] = pydantic.Field(default_factory=dict)
    """Updated region for each detector in this visit."""

    @classmethod
    def update_dimension_records(
        cls,
        butler: Butler,
        instrument: str,
        *where: Any,
        dataset_type: str = "visit_geometry",
        batch_size: int = 100,
        **kwargs: Any,
    ) -> int:
        """Update the dimension records in a data repository by querying for
        and reading visit geometry datasets.

        Parameters
        ----------
        butler : `lsst.daf.butler.Butler`
            Butler client.  Must be writable and have its default collections
            set to the location of the visit geometry datasets.
        instrument : `str`
            Name of the instrument whose records will be updated.
        *where
            Positional argument query constraint terms, forwarded to
            `lsst.daf.butler.queries.Query.where` (e.g. an expression string or
            data ID).
        dataset_type : `str`, optional
            Dataset type name for the visit geometry datasets.
        batch_size : `int`, optional
            Number of visits to process at once.
        **kwargs
            Keyword-only query constraint terms, forwarded to
            `lsst.daf.butler.queries.Query.where` (e.g. data ID key-value
            pairs).

        Returns
        -------
        n_visits : `int`
            The number of visits updated.
        """
        with butler.query() as query:
            query = query.join_dataset_search(dataset_type)
            query = query.where(instrument=instrument, *where, **kwargs)
            _LOG.verbose("Querying for %s datasets.", dataset_type)
            all_refs = list(query.datasets(dataset_type))
            all_refs.sort(key=lambda ref: ref.dataId["visit"])
            _LOG.verbose("Querying for visit records (%d expected).", len(all_refs))
            old_visit_records = {r.id: r for r in query.dimension_records("visit")}
            _LOG.verbose("Querying for exposure records (%d expected).", len(all_refs))
            old_exposure_records = {r.id: r for r in query.dimension_records("exposure")}
        vdr_record_cls = butler.dimensions.elements["visit_detector_region"].RecordClass
        for ref_batch in itertools.batched(all_refs, batch_size):
            _LOG.verbose(
                "Updating records for %d visits with IDs between %d and %d.",
                len(ref_batch),
                ref_batch[0].dataId["visit"],
                ref_batch[-1].dataId["visit"],
            )
            new_exposure_records: list[DimensionRecord] = []
            new_visit_records: list[DimensionRecord] = []
            new_vdr_records: list[DimensionRecord] = []
            for ref in ref_batch:
                visit_id = ref.dataId["visit"]
                visit_geometry: VisitGeometry = butler.get(ref)
                new_exposure_records.append(
                    cls._make_updated_dimension_record(
                        old_exposure_records[visit_id],
                        tracking_ra=visit_geometry.boresight_ra,
                        tracking_dec=visit_geometry.boresight_dec,
                        sky_angle=visit_geometry.orientation,
                    )
                )
                new_visit_records.append(
                    cls._make_updated_dimension_record(
                        old_visit_records[visit_id], region=visit_geometry.visit_region
                    )
                )
                for detector_id, vdr_region in visit_geometry.detector_regions.items():
                    new_vdr_records.append(
                        vdr_record_cls(
                            instrument=instrument, visit=visit_id, detector=detector_id, region=vdr_region
                        )
                    )
            with butler.transaction():
                butler.registry.insertDimensionData("exposure", *new_exposure_records, replace=True)
                butler.registry.insertDimensionData("visit", *new_visit_records, replace=True)
                butler.registry.insertDimensionData("visit_detector_region", *new_vdr_records, replace=True)
        _LOG.info("Records for %d visits updated.", len(all_refs))
        return len(all_refs)

    @staticmethod
    def _make_updated_dimension_record(original: DimensionRecord, **kwargs: Any) -> DimensionRecord:
        as_dict = original.toDict()
        as_dict.update(**kwargs)
        return type(original)(**as_dict)
