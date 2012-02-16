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
        try:
            subprocess.call(["rm", "-rf", "testOutput"])
            subprocess.call(["rm", "-rf", "testInput1"])
            subprocess.call(["rm", "-rf", "testInput2"])
        except:
            pass

    def tearDown(self):
        subprocess.call(["rm", "-rf", "testOutput"])
        subprocess.call(["rm", "-rf", "testInput1"])
        subprocess.call(["rm", "-rf", "testInput2"])

    def testCreateOutputRoot(self):
        mapper = MinMapper1(outputRoot="testOutput")
        self.assert_(os.path.exists("testOutput"))
        self.assert_(os.path.isdir("testOutput"))
        self.assert_(os.path.islink("testOutput/MinMapper1.paf"))
        self.assert_(os.path.islink("testOutput/outputRoot.py"))

    def testReuseOutputRoot(self):
        mapper = MinMapper1(outputRoot="testOutput")
        self.assert_(os.path.exists("testOutput"))
        self.assert_(os.path.isdir("testOutput"))
        self.assert_(os.path.islink("testOutput/MinMapper1.paf"))
        self.assert_(os.path.islink("testOutput/outputRoot.py"))

        mapper = MinMapper1(root="testOutput", outputRoot="testOutput")
        self.assert_(os.path.exists("testOutput"))
        self.assert_(os.path.isdir("testOutput"))
        self.assert_(os.path.islink("testOutput/MinMapper1.paf"))
        self.assert_(os.path.islink("testOutput/outputRoot.py"))

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
        self.assert_(os.path.islink("testOutput/foo"))
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
