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
import lsst.utils.tests as utilsTests

import lsst.daf.persistence as dafPersist
import lsst.daf.butlerUtils as butlerUtils
import lsst.pex.policy as pexPolicy
import lsst.afw.image

class MinMapper2(butlerUtils.CameraMapper):
    def __init__(self):
        policy = pexPolicy.Policy.createPolicy("tests/MinMapper2.paf")
        butlerUtils.CameraMapper.__init__(self, policy=policy,
                repositoryDir="tests", root="tests",
                registry="tests/cfhtls.sqlite3")
        return

    def _transformId(self, dataId):
        return dataId

    def _extractDetectorName(self, dataId):
        return "Detector"
    
class DM329TestCase(unittest.TestCase):

    def testHdu(self):
        mapper = MinMapper2()
        # Load just the mask plane into an ImageF
        loc = mapper.map("other", dict(ccd=35, hdu=3))
        self.assertEqual(loc.getLocations(), ["tests/bar-35.fits[3]"])
        butler = dafPersist.ButlerFactory(mapper=mapper).create()
        image = butler.get("other", ccd=35, hdu=3, immediate=True)
        self.assertEqual(type(image), lsst.afw.image.ImageF)
        self.assertEqual(image.getHeight(), 2024)
        self.assertEqual(image.getWidth(), 2248)
        self.assertEqual(image.get(200, 25), 20.0)
        self.assertEqual(image.get(200, 26), 0.0)

def suite():
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(DM329TestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(), shouldExit)

if __name__ == '__main__':
    run(True)
