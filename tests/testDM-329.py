#!/usr/bin/env python

#
# LSST Data Management System
# Copyright 2014 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#


import unittest
import lsst.utils.tests

import lsst.afw.image
import lsst.daf.persistence as dafPersist
import lsst.obs.base
import lsst.pex.policy as pexPolicy
from lsst.utils import getPackageDir

import os

testDir = os.path.relpath(os.path.join(getPackageDir('obs_base'), 'tests'))


class MinMapper2(lsst.obs.base.CameraMapper):
    packageName = 'larry'

    def __init__(self):
        policy = pexPolicy.Policy.createPolicy(os.path.join(testDir, 'MinMapper2.paf'))
        lsst.obs.base.CameraMapper.__init__(self,
                                            policy=policy,
                                            repositoryDir=testDir,
                                            root=testDir,
                                            registry=os.path.join(testDir, 'cfhtls.sqlite3'))
        return

    def _transformId(self, dataId):
        return dataId

    def _extractDetectorName(self, dataId):
        return "Detector"


class DM329TestCase(unittest.TestCase):

    def testHdu(self):
        mapper = MinMapper2()
        butler = dafPersist.ButlerFactory(mapper=mapper).create()
        # HDU 0 returns (header-only) primary array
        # Note difference from afw Image constructor with hdu=0
        # HDU 1 returns first image plane
        # HDU 2 returns mask plane
        # HDU 3 returns variance plane
        for i in (1, 2, 3):
            loc = mapper.map("other", dict(ccd=35, hdu=i))
            expectedLocations = ["bar-35.fits[%d]" % (i,)]
            self.assertEqual(loc.getStorage().root, testDir)
            self.assertEqual(loc.getLocations(), expectedLocations)
            image = butler.get("other", ccd=35, hdu=i, immediate=True)
            self.assertIsInstance(image, lsst.afw.image.ImageF)
            self.assertEqual(image.getHeight(), 2024)
            self.assertEqual(image.getWidth(), 2248)
            self.assertEqual(image.get(200, 25), (0.0, 20.0, 0.0)[i-1])
            self.assertAlmostEqual(image.get(200, 26), (1.20544, 0.0, 5.82185)[i-1], places=5)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
