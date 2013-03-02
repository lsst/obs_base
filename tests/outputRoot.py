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

import cPickle
import glob
import os
import subprocess

import lsst.afw.geom as afwGeom
import lsst.pex.policy as pexPolicy
import lsst.daf.base as dafBase
import lsst.daf.persistence as dafPersist
import lsst.daf.butlerUtils as butlerUtils

class MinMapper1(butlerUtils.CameraMapper):
    def __init__(self, root="tests", outputRoot=None):
        policy = pexPolicy.Policy.createPolicy("tests/MinMapper1.paf")
        butlerUtils.CameraMapper.__init__(self,
                policy=policy, repositoryDir="tests", root=root,
                outputRoot=outputRoot)
        return

    def std_x(self, item, dataId):
        return float(item)

class OutputRootTestCase(unittest.TestCase):
    """A test case for output roots."""

    def setUp(self):
        # Note: these will succeed even if the directories don't exist.
        subprocess.call(["rm", "-rf", "testOutput"])
        subprocess.call(["rm", "-rf", "testOutput2"])
        subprocess.call(["rm", "-rf", "testInput1"])
        subprocess.call(["rm", "-rf", "testInput2"])

    def tearDown(self):
        subprocess.call(["rm", "-rf", "testOutput", "testOutput2",
            "testOutput3", "testInput1", "testInput2"])

    def testCreateOutputRoot(self):
        mapper = MinMapper1(outputRoot="testOutput")
        self.assert_(os.path.exists("testOutput"))
        self.assert_(os.path.isdir("testOutput"))
        self.assert_(os.path.islink("testOutput/_parent"))
        self.assert_(os.path.exists("testOutput/_parent/MinMapper1.paf"))
        self.assert_(os.path.exists("testOutput/_parent/outputRoot.py"))

    def testParentNormal(self):
        mapper1 = MinMapper1(outputRoot="testOutput")
        mapper2 = MinMapper1(root="testOutput", outputRoot="testOutput2")
        loc = mapper1.map("x", dict(sensor="1,1"), write=True)
        self.assertEqual(loc.getPythonType(), "lsst.afw.geom.BoxI")
        self.assertEqual(loc.getCppType(), "BoxI")
        self.assertEqual(loc.getStorageName(), "PickleStorage")
        self.assertEqual(loc.getLocations(), ["testOutput/foo-1,1.pickle"])
        self.assertEqual(loc.getAdditionalData().toString(), "sensor = \"1,1\"\n")
        box = afwGeom.BoxI(afwGeom.PointI(0, 1), afwGeom.PointI(2, 3))
        with open(loc.getLocations()[0], "w") as f:
            cPickle.dump(box, f)

        loc = mapper2.map("x", dict(sensor="1,1"))
        self.assertEqual(loc.getPythonType(), "lsst.afw.geom.BoxI")
        self.assertEqual(loc.getCppType(), "BoxI")
        self.assertEqual(loc.getStorageName(), "PickleStorage")
        self.assertEqual(loc.getLocations(), ["testOutput2/_parent/foo-1,1.pickle"])
        self.assertEqual(loc.getAdditionalData().toString(), "sensor = \"1,1\"\n")

    def testParentTrailingSlash2527(self):
        mapper1 = MinMapper1(outputRoot="testOutput/")
        mapper2 = MinMapper1(root="testOutput", outputRoot="testOutput2/")
        os.symlink("testOutput2", "testOutput3")

        loc = mapper1.map("x", dict(sensor="1,1"), write=True)
        self.assertEqual(loc.getPythonType(), "lsst.afw.geom.BoxI")
        self.assertEqual(loc.getCppType(), "BoxI")
        self.assertEqual(loc.getStorageName(), "PickleStorage")
        self.assertEqual(loc.getLocations(), ["testOutput/foo-1,1.pickle"])
        self.assertEqual(loc.getAdditionalData().toString(), "sensor = \"1,1\"\n")
        box = afwGeom.BoxI(afwGeom.PointI(0, 1), afwGeom.PointI(2, 3))
        with open(loc.getLocations()[0], "w") as f:
            cPickle.dump(box, f)

        parent = mapper2._parentSearch("testOutput3/foo-1,1.pickle")
        self.assertEqual(parent, "testOutput3/_parent/foo-1,1.pickle")

        loc = mapper2.map("x", dict(sensor="1,1"))
        self.assertEqual(loc.getPythonType(), "lsst.afw.geom.BoxI")
        self.assertEqual(loc.getCppType(), "BoxI")
        self.assertEqual(loc.getStorageName(), "PickleStorage")
        self.assertEqual(loc.getLocations(), ["testOutput2/_parent/foo-1,1.pickle"])
        self.assertEqual(loc.getAdditionalData().toString(), "sensor = \"1,1\"\n")

    def testReuseOutputRoot(self):
        mapper = MinMapper1(outputRoot="testOutput")
        self.assert_(os.path.exists("testOutput"))
        self.assert_(os.path.isdir("testOutput"))
        self.assert_(os.path.islink("testOutput/_parent"))
        self.assert_(os.path.exists("testOutput/_parent/MinMapper1.paf"))
        self.assert_(os.path.exists("testOutput/_parent/outputRoot.py"))

        self.assertRaises(RuntimeError, MinMapper1,
                root="testOutput", outputRoot="testOutput")

        mapper = MinMapper1(root="testOutput", outputRoot="testOutput2")
        self.assert_(os.path.exists("testOutput2"))
        self.assert_(os.path.isdir("testOutput2"))
        self.assert_(os.path.islink("testOutput2/_parent"))
        self.assert_(os.path.exists("testOutput2/_parent/_parent/MinMapper1.paf"))
        self.assert_(os.path.exists("testOutput2/_parent/_parent/outputRoot.py"))

    def testDiffInput(self):
        os.mkdir("testInput1")
        with open("testInput1/foo", "w") as f:
            pass
        os.mkdir("testInput2")
        with open("testInput2/foo", "w") as f:
            pass
        mapper = MinMapper1(root="testInput1", outputRoot="testOutput")
        self.assert_(os.path.exists("testOutput"))
        self.assert_(os.path.isdir("testOutput"))
        self.assert_(os.path.islink("testOutput/_parent"))
        self.assert_(os.path.exists("testOutput/_parent/foo"))
        self.assertRaises(RuntimeError, MinMapper1,
                root="testInput2", outputRoot="testOutput")
        os.unlink("testInput1/foo")
        os.unlink("testInput2/foo")
        os.rmdir("testInput1")
        os.rmdir("testInput2")


def suite():
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(OutputRootTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(), shouldExit)

if __name__ == '__main__':
    run(True)
