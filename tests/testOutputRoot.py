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


from future import standard_library
standard_library.install_aliases()
import unittest
import lsst.utils.tests
import pickle
import os
import subprocess
import lsst.afw.geom as afwGeom
import lsst.pex.policy as pexPolicy
import lsst.daf.persistence as dafPersist
import lsst.obs.base

# Define paths used for testing
testPath = os.path.abspath(os.path.dirname(__file__))
testOutput = os.path.join(testPath, "testOutput")
testOutput2 = os.path.join(testPath, "testOutput2")
testOutput3 = os.path.join(testPath, "testOutput3")
testInput1 = os.path.join(testPath, "testInput1")
testInput2 = os.path.join(testPath, "testInput2")


def setup_module(module):
    lsst.utils.tests.init()


class MinMapper1(lsst.obs.base.CameraMapper):
    packageName = 'larry'

    def __init__(self, root=testPath, outputRoot=None):
        policy = pexPolicy.Policy.createPolicy(os.path.join(testPath, "MinMapper1.paf"))
        lsst.obs.base.CameraMapper.__init__(self,
                                            policy=policy, repositoryDir=testPath, root=root,
                                            outputRoot=outputRoot)
        return

    def std_x(self, item, dataId):
        return float(item)


class OutputRootTestCase(unittest.TestCase):
    """A test case for output roots."""

    def setUp(self):
        # Note: these will succeed even if the directories don't exist.
        subprocess.call(["rm", "-rf", testOutput])
        subprocess.call(["rm", "-rf", testOutput2])
        subprocess.call(["rm", "-rf", testInput1])
        subprocess.call(["rm", "-rf", testInput2])

    def tearDown(self):
        subprocess.call(["rm", "-rf", testOutput, testOutput2,
                         testOutput3, testInput1, testInput2])

    def testCreateOutputRoot(self):
        MinMapper1(outputRoot=testOutput)
        self.assertTrue(os.path.exists(testOutput))
        self.assertTrue(os.path.isdir(testOutput))
        self.assertTrue(os.path.islink(os.path.join(testOutput, "_parent")))
        self.assertTrue(os.path.exists(os.path.join(testOutput, "_parent", "MinMapper1.paf")))
        self.assertTrue(os.path.exists(os.path.join(testOutput, "_parent", "testOutputRoot.py")))

    def testParentNormal(self):
        mapper1 = MinMapper1(outputRoot=testOutput)
        mapper2 = MinMapper1(root=testOutput, outputRoot=testOutput2)
        loc = mapper1.map("x", dict(sensor="1,1"), write=True)
        self.assertEqual(loc.getPythonType(), "lsst.afw.geom.BoxI")
        self.assertEqual(loc.getCppType(), "BoxI")
        self.assertEqual(loc.getStorageName(), "PickleStorage")
        self.assertEqual(loc.getLocations(), ["foo-1,1.pickle"])
        self.assertEqual(loc.getAdditionalData().toString(), "sensor = \"1,1\"\n")
        box = afwGeom.BoxI(afwGeom.PointI(0, 1), afwGeom.PointI(2, 3))
        with open(os.path.join(testOutput, loc.getLocations()[0]), "wb") as f:
            pickle.dump(box, f)

        loc = mapper2.map("x", dict(sensor="1,1"))
        self.assertEqual(loc.getPythonType(), "lsst.afw.geom.BoxI")
        self.assertEqual(loc.getCppType(), "BoxI")
        self.assertEqual(loc.getStorageName(), "PickleStorage")
        self.assertEqual(loc.getLocations(), [os.path.join("_parent", "foo-1,1.pickle")])
        self.assertEqual(loc.getStorage().getRoot(), testOutput2)
        self.assertEqual(loc.getAdditionalData().toString(), "sensor = \"1,1\"\n")

    def testParentTrailingSlash2527(self):
        # todo these shouldn't be commented out, I think the test wants the trailing slash.
        #mapper1 = MinMapper1(outputRoot="testOutput/")
        mapper1 = MinMapper1(outputRoot=testOutput)
        #mapper2 = MinMapper1(root="testOutput", outputRoot="testOutput2/")
        mapper2 = MinMapper1(root=testOutput, outputRoot=testOutput2)
        os.symlink('testOutput2', testOutput3)

        loc = mapper1.map("x", dict(sensor="1,1"), write=True)
        self.assertEqual(loc.getPythonType(), "lsst.afw.geom.BoxI")
        self.assertEqual(loc.getCppType(), "BoxI")
        self.assertEqual(loc.getStorageName(), "PickleStorage")
        self.assertEqual(loc.getLocations(), ["foo-1,1.pickle"])
        self.assertEqual(loc.getStorage().getRoot(), testOutput)
        self.assertEqual(loc.getAdditionalData().toString(), "sensor = \"1,1\"\n")
        box = afwGeom.BoxI(afwGeom.PointI(0, 1), afwGeom.PointI(2, 3))
        with open(os.path.join(loc.getStorage().getRoot(), loc.getLocations()[0]), "wb") as f:
            pickle.dump(box, f)

        parent = mapper2._parentSearch(os.path.join(testOutput3, "foo-1,1.pickle"))
        self.assertEqual(parent, [os.path.join("_parent", "foo-1,1.pickle")])

        loc = mapper2.map("x", dict(sensor="1,1"))
        self.assertEqual(loc.getPythonType(), "lsst.afw.geom.BoxI")
        self.assertEqual(loc.getCppType(), "BoxI")
        self.assertEqual(loc.getStorageName(), "PickleStorage")
        self.assertEqual(loc.getLocations(), [os.path.join("_parent", "foo-1,1.pickle")])
        self.assertEqual(loc.getStorage().getRoot(), testOutput2)
        self.assertEqual(loc.getAdditionalData().toString(), "sensor = \"1,1\"\n")

    def testReuseOutputRoot(self):
        MinMapper1(outputRoot=testOutput)
        self.assertTrue(os.path.exists(testOutput))
        self.assertTrue(os.path.isdir(testOutput))
        self.assertTrue(os.path.islink(os.path.join(testOutput, "_parent")))
        self.assertTrue(os.path.exists(os.path.join(testOutput, "_parent", "MinMapper1.paf")))
        self.assertTrue(os.path.exists(os.path.join(testOutput, "_parent", "testOutputRoot.py")))

        MinMapper1(root=testOutput, outputRoot=testOutput2)
        self.assertTrue(os.path.exists(testOutput2))
        self.assertTrue(os.path.isdir(testOutput2))
        self.assertTrue(os.path.islink(os.path.join(testOutput2, "_parent")))
        self.assertTrue(os.path.exists(os.path.join(testOutput2, "_parent", "_parent", "MinMapper1.paf")))
        self.assertTrue(os.path.exists(os.path.join(testOutput2, "_parent", "_parent", "testOutputRoot.py")))

    def testDiffInput(self):
        os.mkdir(testInput1)
        with open(os.path.join(testInput1, "foo"), "w"):
            pass
        os.mkdir(testInput2)
        with open(os.path.join(testInput2, "foo"), "w"):
            pass
        MinMapper1(root=testInput1, outputRoot=testOutput)
        self.assertTrue(os.path.exists(testOutput))
        self.assertTrue(os.path.isdir(testOutput))
        self.assertTrue(os.path.islink(os.path.join(testOutput, "_parent")))
        self.assertTrue(os.path.exists(os.path.join(testOutput, "_parent", "foo")))
        self.assertRaises(RuntimeError, MinMapper1,
                          root=testInput2, outputRoot=testOutput)
        os.unlink(os.path.join(testInput1, "foo"))
        os.unlink(os.path.join(testInput2, "foo"))
        os.rmdir(testInput1)
        os.rmdir(testInput2)

    @unittest.expectedFailure
    def testBackup(self):
        mapper1 = MinMapper1(outputRoot=testOutput)
        butler1 = dafPersist.Butler(outputs=dafPersist.RepositoryArgs(mode='w',
                                                                      root=testOutput,
                                                                      mapper=mapper1))
        b1 = afwGeom.Box2I(afwGeom.Point2I(3, 4), afwGeom.Point2I(7, 6))
        butler1.put(b1, "x")
        self.assertTrue(os.path.exists(os.path.join(testOutput, "foo-1,1.pickle")))
        b2 = afwGeom.Box2I(b1)
        b2.grow(1)
        butler1.put(b2, "x", doBackup=True)
        self.assertTrue(os.path.exists(os.path.join(testOutput, "foo-1,1.pickle")))
        self.assertTrue(os.path.exists(os.path.join(testOutput, "foo-1,1.pickle~1")))
        mapper2 = MinMapper1(root=testOutput, outputRoot=testOutput2)
        butler2 = dafPersist.Butler(
            # MinMapper is a little unconventional in that it takes its root and output root as separate
            # arguments, meaning that in effect it's a mapper for 2 different repositories
            outputs=dafPersist.RepositoryArgs(
                mode='rw',
                root=testOutput2,
                mapper=mapper2),
        )
        b3 = afwGeom.Box2I(b2)
        b3.grow(1)
        butler2.put(b3, "x", doBackup=True)

        self.assertTrue(os.path.exists(os.path.join(testOutput2, "foo-1,1.pickle")))
        self.assertTrue(os.path.exists(os.path.join(testOutput2, "foo-1,1.pickle~1")))
        self.assertTrue(os.path.exists(os.path.join(testOutput2, "foo-1,1.pickle~2")))
        with open(os.path.join(testOutput2, "foo-1,1.pickle~2"), 'rb') as f:
            b1Check = pickle.load(f)
            self.assertEqual(b1Check, b1)
        with open(os.path.join(testOutput2, "foo-1,1.pickle~1"), 'rb') as f:
            b2Check = pickle.load(f)
            self.assertEqual(b2Check, b2)
        with open(os.path.join(testOutput2, "foo-1,1.pickle"), 'rb') as f:
            b3Check = pickle.load(f)
            self.assertEqual(b3Check, b3)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == '__main__':
    lsst.utils.tests.init()
    unittest.main()
