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

from lsst.obs.base import ExposureIdInfo


class ExposureIdInfoTestCase(unittest.TestCase):
    def testExposureIdInfo(self):
        idInfo = ExposureIdInfo(12345, 32, 128)
        self.assertEqual(idInfo.unusedBits, 96)

        self.assertIn("=12345", str(idInfo))

        with self.assertRaises(RuntimeError):
            ExposureIdInfo.fromDataId({})

        with self.assertRaises(RuntimeError):
            ExposureIdInfo(1e6, 4)

        with self.assertRaises(RuntimeError):
            ExposureIdInfo(12345, 64, 32)


if __name__ == "__main__":
    unittest.main()
