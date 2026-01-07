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

import tempfile
import unittest

import numpy as np

from lsst.daf.butler import Butler, DatasetType, DimensionRecord
from lsst.obs.base.visit_geometry import VisitGeometry
from lsst.sphgeom import Angle, ConvexPolygon, LonLat, UnitVector3d


class VisitGeometryTestCase(unittest.TestCase):
    """Tests for VisitGeometry."""

    def setUp(self):
        root = self.enterContext(tempfile.TemporaryDirectory(ignore_cleanup_errors=True))
        Butler.makeRepo(root)
        self.butler: Butler = self.enterContext(Butler.from_config(root, run="geometry"))
        self.rng = np.random.default_rng(500)
        self.rot_about = UnitVector3d(*self.rng.uniform(0.0, 1.0, size=3).tolist())
        self.rot_amount = Angle.fromDegrees(self.rng.uniform(0.0, 5.0))

    def test_round_trip(self):
        """Test round-tripping through a JSON string."""
        self.butler.import_(filename="resource://lsst.daf.butler/tests/registry_data/base.yaml")
        self.butler.import_(filename="resource://lsst.daf.butler/tests/registry_data/spatial.yaml")
        (visit_record,) = self.butler.query_dimension_records("visit", visit=1, instrument="Cam1")
        visit_detector_region_records = self.butler.query_dimension_records(
            "visit_detector_region", visit=1, instrument="Cam1"
        )
        vg1 = VisitGeometry(
            boresight_ra=45.0,
            boresight_dec=60.0,
            orientation=30.0,
            visit_region=visit_record.region,
            detector_regions={r.detector: r.region for r in visit_detector_region_records},
        )
        json_data = vg1.model_dump_json()
        vg2 = VisitGeometry.model_validate_json(json_data)
        self.assertEqual(vg1.boresight_ra, vg2.boresight_ra)
        self.assertEqual(vg1.boresight_dec, vg2.boresight_dec)
        self.assertEqual(vg1.orientation, vg2.orientation)
        self.assertEqual(vg1.visit_region, vg2.visit_region)
        self.assertEqual(vg1.detector_regions, vg2.detector_regions)

    def test_update_dimension_records(self):
        """Test the update_dimension_records method."""
        self.butler.import_(filename="resource://lsst.daf.butler/tests/registry_data/base.yaml")
        self.butler.import_(filename="resource://lsst.daf.butler/tests/registry_data/spatial.yaml")
        old_visit_records = {r.id: r for r in self.butler.query_dimension_records("visit")}
        old_exposure_records = {
            r.id: self._invent_and_insert_exposure_record(r) for r in old_visit_records.values()
        }
        self.butler.registry.registerDatasetType(
            DatasetType("visit_geometry", self.butler.dimensions.conform(["visit"]), "VisitGeometry")
        )
        visit_geometries: dict[int, VisitGeometry] = {}
        self.assertEqual(len(old_visit_records), 2)
        for visit_id, visit_record in old_visit_records.items():
            exposure_record = old_exposure_records[visit_id]
            old_boresight = LonLat.fromDegrees(exposure_record.tracking_ra, exposure_record.tracking_dec)
            new_boresight = LonLat(self._modify_point(UnitVector3d(old_boresight)))
            visit_geometry = VisitGeometry(
                boresight_ra=new_boresight.getLon().asDegrees(),
                boresight_dec=new_boresight.getLat().asDegrees(),
                # We *mostly* try to modify the geometry self-consistently, but
                # the orientation angle is just arbitrary.  Nothing in the test
                # right now requires any kind of consistency between the
                # regions and the exposure angle fields.
                orientation=self.rng.uniform(0.0, 5.0),
                visit_region=self._modify_polygon(visit_record.region),
                detector_regions={
                    r.detector: self._modify_polygon(r.region)
                    for r in self.butler.query_dimension_records(
                        "visit_detector_region", data_id=visit_record.dataId
                    )
                },
            )
            visit_geometries[visit_id] = visit_geometry
            self.butler.put(visit_geometry, "visit_geometry", visit_record.dataId)
        n_updated = VisitGeometry.update_dimension_records(self.butler, instrument="Cam1")
        self.assertEqual(n_updated, 2)
        for new_visit_record in self.butler.query_dimension_records("visit"):
            visit_id = new_visit_record.id
            old_visit_record = old_visit_records[visit_id]
            visit_geometry = visit_geometries[visit_id]
            self.assertNotEqual(new_visit_record.region, old_visit_record.region)
            self.assertEqual(new_visit_record.region, visit_geometry.visit_region)
            for new_vdr_record in self.butler.query_dimension_records(
                "visit_detector_region", data_id=new_visit_record.dataId
            ):
                self.assertEqual(
                    new_vdr_record.region, visit_geometry.detector_regions[new_vdr_record.detector]
                )
        for new_exposure_record in self.butler.query_dimension_records("exposure"):
            visit_id = new_exposure_record.id
            visit_geometry = visit_geometries[visit_id]
            self.assertEqual(new_exposure_record.tracking_ra, visit_geometry.boresight_ra)
            self.assertEqual(new_exposure_record.tracking_dec, visit_geometry.boresight_dec)
            self.assertEqual(new_exposure_record.sky_angle, visit_geometry.orientation)

    def _invent_and_insert_exposure_record(self, visit_record: DimensionRecord) -> DimensionRecord:
        group_record = self.butler.dimensions["group"].RecordClass(
            instrument=visit_record.instrument, name=str(visit_record.id)
        )
        self.butler.registry.insertDimensionData("group", group_record, skip_existing=True)
        center = LonLat(visit_record.region.getCentroid())
        exposure_record = self.butler.dimensions["exposure"].RecordClass(
            instrument=visit_record.instrument,
            id=visit_record.id,
            day_obs=visit_record.day_obs,
            group=group_record.name,
            physical_filter=visit_record.physical_filter,
            tracking_ra=center.getLon().asDegrees(),
            tracking_dec=center.getLat().asDegrees(),
            sky_angle=float(self.rng.uniform(0.0, 5.0)),
        )
        self.butler.registry.insertDimensionData("exposure", exposure_record)
        visit_definition_record = self.butler.dimensions["visit_definition"].RecordClass(
            instrument=visit_record.instrument,
            visit=visit_record.id,
            exposure=exposure_record.id,
        )
        self.butler.registry.insertDimensionData("visit_definition", visit_definition_record)
        return exposure_record

    def _modify_polygon(self, polygon: ConvexPolygon) -> ConvexPolygon:
        return ConvexPolygon([self._modify_point(v) for v in polygon.getVertices()])

    def _modify_point(self, u: UnitVector3d) -> UnitVector3d:
        return u.rotatedAround(self.rot_about, self.rot_amount)


if __name__ == "__main__":
    unittest.main()
