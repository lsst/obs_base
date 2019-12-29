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
import os

import lsst.utils.tests
import lsst.afw.image
import lsst.daf.persistence as dafPersist
import lsst.obs.base


ROOT = os.path.abspath(os.path.dirname(__file__))


class MinMapper2(lsst.obs.base.CameraMapper):
    packageName = 'larry'

    def __init__(self):
        policy = dafPersist.Policy(os.path.join(ROOT, 'MinMapper2.yaml'))
        lsst.obs.base.CameraMapper.__init__(self,
                                            policy=policy,
                                            repositoryDir=ROOT,
                                            root=ROOT,
                                            registry=os.path.join(ROOT, 'cfhtls.sqlite3'))
        return

    def _transformId(self, dataId):
        return dataId

    def _extractDetectorName(self, dataId):
        return "Detector"

    @classmethod
    def getPackageDir(cls):
        return "/path/to/nowhere"


class DM329TestCase(unittest.TestCase):

    def testHdu(self):
        mapper = MinMapper2()
        butler = dafPersist.ButlerFactory(mapper=mapper).create()
        # HDU INT_MIN returns primary array (skipping empty PDU)
        # HDU 0 returns (header-only) PDU
        # HDU 1 returns first image plane
        # HDU 2 returns mask plane
        # HDU 3 returns variance plane
        for i in (1, 2, 3):
            loc = mapper.map("other", dict(ccd=35, hdu=i))
            expectedLocations = ["bar-35.fits[%d]" % (i,)]
            self.assertEqual(loc.getStorage().root, ROOT)
            self.assertEqual(loc.getLocations(), expectedLocations)
            image = butler.get("other", ccd=35, hdu=i, immediate=True)
            self.assertIsInstance(image, lsst.afw.image.ImageF)
            self.assertEqual(image.getHeight(), 2024)
            self.assertEqual(image.getWidth(), 2248)
            self.assertEqual(image[200, 25, lsst.afw.image.LOCAL], (0.0, 20.0, 0.0)[i-1])
            self.assertAlmostEqual(image[200, 26, lsst.afw.image.LOCAL], (1.20544, 0.0, 5.82185)[i-1],
                                   places=5)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
