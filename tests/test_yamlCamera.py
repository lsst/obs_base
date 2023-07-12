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

from lsst.obs.base.yamlCamera import makeCamera
from lsst.resources import ResourcePath


class YamlCameraTestCase(unittest.TestCase):
    """Test the YAML camera geometry."""

    def setUp(self):
        self.cameraFile = ResourcePath("resource://lsst.obs.base/test/dummycam.yaml")

    def test_basics(self):
        """Basic test of yaml camera construction."""
        with self.cameraFile.as_local() as local_file:
            camera = makeCamera(local_file.ospath)

        self.assertEqual(len(camera), 2)
        self.assertEqual(camera[0].getName(), "RXX_S00")
        self.assertEqual(camera[1].getName(), "RXX_S01")


if __name__ == "__main__":
    unittest.main()
