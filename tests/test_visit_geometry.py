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

import unittest

from lsst.obs.base.visit_geometry import VisitGeometry
from lsst.pipe.base.tests.mocks import InMemoryRepo


class VisitGeometryTestCase(unittest.TestCase):
    """Tests for VisitGeometry."""

    def setUp(self):
        self.repo: InMemoryRepo = self.enterContext(InMemoryRepo("base.yaml", "spatial.yaml"))

    def test_round_trip(self):
        """Test round-tripping through a JSON string."""
        (visit_record,) = self.repo.butler.query_dimension_records("visit", visit=1, instrument="Cam1")
        visit_detector_region_records = self.repo.butler.query_dimension_records(
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


if __name__ == "__main__":
    unittest.main()
