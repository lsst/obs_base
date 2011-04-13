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

import lsst.afw.geom as afwGeom
import lsst.pex.policy as pexPolicy
import lsst.daf.base as dafBase
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
                         set(["x", "x_filename", "badSourceHist",
                             "badSourceHist_filename", "camera", "skypolicy"]))

    def testMap(self):
        loc = self.mapper.map("x", {"sensor": "1,1"})
        self.assertEqual(loc.getPythonType(), "lsst.afw.image.BBox")
        self.assertEqual(loc.getCppType(), "BBox")
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
        self.assertEqual(self.mapper.canStandardize("badSourceHist"), False)
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
                         set(["flat", "flat_filename", "raw", "raw_md",
                             "raw_filename", "raw_sub",
                              "camera", "src", "src_filename", "skypolicy"]))

    def testMap(self):
        loc = self.mapper.map("raw", {"ccd": 13})
        self.assertEqual(loc.getPythonType(), "lsst.afw.image.ExposureU")
        self.assertEqual(loc.getCppType(), "ImageU")
        self.assertEqual(loc.getStorageName(), "FitsStorage")
        self.assertEqual(loc.getLocations(), ["tests/foo-13.fits"])
        self.assertEqual(loc.getAdditionalData().toString(),
                "ccd = 13\n")

    def testSubMap(self):
        bbox = afwGeom.Box2I(afwGeom.Point2I(200, 100), afwGeom.Extent2I(300, 400))
        loc = self.mapper.map("raw_sub", {"ccd": 13, "bbox": bbox})
        self.assertEqual(loc.getPythonType(), "lsst.afw.image.ExposureU")
        self.assertEqual(loc.getCppType(), "ImageU")
        self.assertEqual(loc.getStorageName(), "FitsStorage")
        self.assertEqual(loc.getLocations(), ["tests/foo-13.fits"])
        self.assertEqual(loc.getAdditionalData().toString(),
                'ccd = 13\nheight = 400\nllcX = 200\nllcY = 100\nwidth = 300\n')

    def testQueryMetadata(self):
        self.assertEqual(
                self.mapper.queryMetadata("raw", "ccd", ["ccd"], None),
                [(x,) for x in xrange(36) if x != 3])

    def testStandardize(self):
        self.assertEqual(self.mapper.canStandardize("raw"), True)
        self.assertEqual(self.mapper.canStandardize("notPresent"), False)

    def testCalib(self):
        loc = self.mapper.map("flat", {"visit": 787650, "ccd": 13})
        self.assertEqual(loc.getPythonType(), "lsst.afw.image.ExposureF")
        self.assertEqual(loc.getCppType(), "ExposureF")
        self.assertEqual(loc.getStorageName(), "FitsStorage")
        self.assertEqual(loc.getLocations(), ["tests/flat-05Am03-fi.fits"])
        self.assertEqual(loc.getAdditionalData().toString(),
                'ccd = 13\nderivedRunId = "05Am03"\nfilter = "i"\nvisit = 787650\n')

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
