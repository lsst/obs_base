#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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

import lsst.pex.policy as pexPolicy
import lsst.daf.persistence as dafPersist
import lsst.daf.butlerUtils as butlerUtils

class MinMapper1(butlerUtils.CameraMapper):
    def __init__(self):
        policy = pexPolicy.Policy.createPolicy("tests/MinMapper1.paf")
        butlerUtils.CameraMapper.__init__(self,
                policy=policy, repositoryDir="tests", root="tests")
        return

    def std_x(self, item, dataId):
        return float(item)

class MinMapper2(butlerUtils.CameraMapper):
    # CalibRoot in policy
    # needCalibRegistry
    def __init__(self):
        policy = pexPolicy.Policy.createPolicy("tests/MinMapper2.paf")
        butlerUtils.CameraMapper.__init__(self,
                policy=policy, repositoryDir="tests", root="tests",
                registry="tests/cfhtls.sqlite3")
        return

    def _transformId(self, dataId):
        return dataId

    def _extractDetectorName(self, dataId):
        return "Detector"

    def std_x(self, item, dataId):
        return float(item)

class Mapper1TestCase(unittest.TestCase):
    """A test case for the mapper used by the data butler."""

    def setUp(self):
        self.mapper = MinMapper1()

    def tearDown(self):
        del self.mapper

    def testGetDatasetTypes(self):
        self.assertEqual(set(self.mapper.getDatasetTypes()),
                set(["x", "badSourceHist", "camera"]))

    def testMap(self):
        loc = self.mapper.map("x", {"sensor": "1,1"})
        self.assertEqual(loc.getPythonType(), "lsst.afw.image.BBox")
        self.assertEqual(loc.getCppType(), "lsst::afw::image::BBox")
        self.assertEqual(loc.getStorageName(), "PickleStorage")
        self.assertEqual(loc.getLocations(), ["tests/foo-1,1.pickle"])
        self.assertEqual(loc.getAdditionalData().toString(),
                "sensor = \"1,1\"\n")

    def testQueryMetadata(self):
        self.assertEqual(
                self.mapper.queryMetadata("x", "sensor", ["sensor"], None),
                [("1,1",)])

    def testStandardize(self):
        self.assertEqual(self.mapper.canStandardize("x"), True)
        self.assertEqual(self.mapper.canStandardize("badSourceHist"), True)
        self.assertEqual(self.mapper.canStandardize("notPresent"), False)
        result = self.mapper.standardize("x", 3, None)
        self.assertEqual(isinstance(result, float), True)
        self.assertEqual(result, 3.0)
        result = self.mapper.standardize("x", 3.14, None)
        self.assertEqual(isinstance(result, float), True)
        self.assertEqual(result, 3.14)
        result = self.mapper.standardize("x", "3.14", None)
        self.assertEqual(isinstance(result, float), True)
        self.assertEqual(result, 3.14)

class Mapper2TestCase(unittest.TestCase):
    """A test case for the mapper used by the data butler."""

    def setUp(self):
        self.mapper = MinMapper2()

    def tearDown(self):
        del self.mapper

    def testGetDatasetTypes(self):
        self.assertEqual(set(self.mapper.getDatasetTypes()),
                set(["flat", "raw", "camera"]))

    def testMap(self):
        loc = self.mapper.map("raw", {"ccd": 13})
        self.assertEqual(loc.getPythonType(), "lsst.afw.image.Exposure")
        self.assertEqual(loc.getCppType(), "lsst::afw::image::Exposure")
        self.assertEqual(loc.getStorageName(), "FitsStorage")
        self.assertEqual(loc.getLocations(), ["tests/foo-13.fits"])
        self.assertEqual(loc.getAdditionalData().toString(),
                "ccd = 13\n")

    def testQueryMetadata(self):
        self.assertEqual(
                self.mapper.queryMetadata("raw", "ccd", ["ccd"], None),
                [(x,) for x in xrange(36) if x != 3])

    def testStandardize(self):
        self.assertEqual(self.mapper.canStandardize("raw"), True)
        self.assertEqual(self.mapper.canStandardize("notPresent"), False)

    def testCalib(self):
        loc = self.mapper.map("flat", {"visit": 787650})
        self.assertEqual(loc.getPythonType(), "lsst.afw.image.Exposure")
        self.assertEqual(loc.getCppType(), "lsst::afw::image::Exposure")
        self.assertEqual(loc.getStorageName(), "FitsStorage")
        self.assertEqual(loc.getLocations(), ["tests/flat-fr-c13-a1.fits"])
        self.assertEqual(loc.getAdditionalData().toString(),
                'filter = "r"\n')

def suite():
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(Mapper1TestCase)
    suites += unittest.makeSuite(Mapper2TestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(), shouldExit)

if __name__ == '__main__':
    run(True)
