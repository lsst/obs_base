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
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import unittest

import lsst.geom as geom
import lsst.obs.base as obsBase
import lsst.utils.tests


class BboxFromIrafTestCase(lsst.utils.tests.TestCase):
    """Demonstrate that we can correctly parse IRAF-style BBOXes"""

    def testValid(self):
        test_data = {
            "[1:1084,1:1024]": geom.BoxI(geom.PointI(0, 0), geom.PointI(1083, 1023)),
            "[0:0,0:0]": geom.BoxI(geom.PointI(-1, -1), geom.PointI(-1, -1)),
        }
        for val, result in test_data.items():
            self.assertEqual(obsBase.bboxFromIraf(val), result)

    def testInvalid(self):
        test_data = {
            "1:1084,1:1024": RuntimeError,
            "(1:1084,1:1024)": RuntimeError,
            ("1:1084", "1:1024"): TypeError,
        }
        for val, err in test_data.items():
            self.assertRaises(err, obsBase.bboxFromIraf, val)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
