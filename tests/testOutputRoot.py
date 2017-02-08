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

    def __init__(self, **kwargs):
        policy = pexPolicy.Policy.createPolicy(os.path.join(testPath, "MinMapper1.paf"))
        lsst.obs.base.CameraMapper.__init__(self,
                                            policy=policy, repositoryDir=testPath, **kwargs)
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
        """Create an input repository and a related output repo, and verify there is a parent relationship
        from the output repo to the input repo."""
        butler = dafPersist.Butler(inputs={'root': testPath, 'mapper': MinMapper1},
                                   outputs=testOutput)
        self.assertTrue(butler)
        self.assertTrue(os.path.exists(testOutput))
        self.assertTrue(os.path.isdir(testOutput))
        self.assertTrue(os.path.exists(os.path.join(testOutput, "repositoryCfg.yaml")))
        cfg = dafPersist.PosixStorage.getRepositoryCfg(testOutput)
        self.assertEqual(len(cfg.parents), 1)
        self.assertEqual(cfg.parents[0], testPath)

    def testParentNormal(self):
        """Test that an object can be found at root location and put into an output location.
        Then test that when the output locaiton is used as an input location, and with a new output location,
        that the object is found in the first output location.
        """
        butler = dafPersist.Butler(inputs={'root': testPath, 'mapper': MinMapper1}, outputs=testOutput)
        mapper1 = butler._repos.inputs()[0].repo._mapper
        loc = mapper1.map("x", dict(sensor="1,1"), write=True)
        self.assertEqual(loc.getPythonType(), "lsst.afw.geom.BoxI")
        self.assertEqual(loc.getCppType(), "BoxI")
        self.assertEqual(loc.getStorageName(), "PickleStorage")
        self.assertEqual(loc.getLocations(), ["foo-1,1.pickle"])
        self.assertEqual(loc.getAdditionalData().toString(), "sensor = \"1,1\"\n")
        box = afwGeom.BoxI(afwGeom.PointI(0, 1), afwGeom.PointI(2, 3))
        butler.put(box, "x", sensor="1,1")
        self.assertTrue(os.path.exists(os.path.join(testOutput, loc.getLocations()[0])))
        del butler

        butler = dafPersist.Butler(inputs={'root': testOutput, 'mapper': MinMapper1}, outputs=testOutput2)
        mapper2 = butler._repos.inputs()[0].repo._mapper
        loc = mapper2.map("x", dict(sensor="1,1"))
        self.assertEqual(loc.getPythonType(), "lsst.afw.geom.BoxI")
        self.assertEqual(loc.getCppType(), "BoxI")
        self.assertEqual(loc.getStorageName(), "PickleStorage")
        self.assertEqual(loc.getLocations(), ["foo-1,1.pickle"])
        self.assertEqual(loc.getStorage().getRoot(), testOutput)
        self.assertEqual(loc.getAdditionalData().toString(), "sensor = \"1,1\"\n")

    def testParentTrailingSlash2527(self):
        """Just like testParentNormal, but put a trailing slash on the root paths.
        Test that an object can be found at root location and put into an output location.
        Then test that when the output locaiton is used as an input location, and with a new output location,
        that the object is found in the first output location."""
        # todo these shouldn't be commented out, I think the test wants the trailing slash.
        butler = dafPersist.Butler(inputs={'root': testPath + '/', 'mapper': MinMapper1},
                                   outputs=testOutput + '/')
        mapper1 = butler._repos.inputs()[0].repo._mapper
        loc = mapper1.map("x", dict(sensor="1,1"), write=True)
        self.assertEqual(loc.getPythonType(), "lsst.afw.geom.BoxI")
        self.assertEqual(loc.getCppType(), "BoxI")
        self.assertEqual(loc.getStorageName(), "PickleStorage")
        self.assertEqual(loc.getLocations(), ["foo-1,1.pickle"])
        self.assertEqual(loc.getStorage().getRoot(), testPath)
        self.assertEqual(loc.getAdditionalData().toString(), "sensor = \"1,1\"\n")
        box = afwGeom.BoxI(afwGeom.PointI(0, 1), afwGeom.PointI(2, 3))
        butler.put(box, "x", sensor="1,1")
        self.assertTrue(os.path.exists(os.path.join(testOutput, loc.getLocations()[0])))
        del butler
        del mapper1

        butler = dafPersist.Butler(inputs={'root': testOutput, 'mapper': MinMapper1},
                                   outputs=testOutput2 + '/')
        mapper2 = butler._repos.inputs()[0].repo._mapper
        loc = mapper2.map("x", dict(sensor="1,1"))
        self.assertEqual(loc.getPythonType(), "lsst.afw.geom.BoxI")
        self.assertEqual(loc.getCppType(), "BoxI")
        self.assertEqual(loc.getStorageName(), "PickleStorage")
        self.assertEqual(loc.getLocations(), ["foo-1,1.pickle"])
        self.assertEqual(loc.getStorage().getRoot(), testOutput)
        self.assertEqual(loc.getAdditionalData().toString(), "sensor = \"1,1\"\n")

    def testReuseOutputRoot(self):
        """Set up an output repositoriy and verify its parent relationship to the input repository.
        Then set up an output repository with the first output as an input, and verify the parent
        relationships."""
        butler = dafPersist.Butler(inputs={'root': testPath, 'mapper': MinMapper1},
                                   outputs=testOutput)
        self.assertTrue(os.path.exists(testOutput))
        self.assertTrue(os.path.isdir(testOutput))
        cfg = dafPersist.Storage.getRepositoryCfg(testOutput)
        self.assertEqual(cfg.parents, [testPath])
        del butler

        butler = dafPersist.Butler(inputs={'root': testOutput, 'mapper': MinMapper1},
                                   outputs=testOutput2)
        self.assertTrue(os.path.exists(testOutput2))
        self.assertTrue(os.path.isdir(testOutput2))
        cfg = dafPersist.Storage.getRepositoryCfg(testOutput2)
        self.assertEqual(cfg.parents, [testOutput])
        del butler

    def testDiffInput(self):
        """Verify that if an output repository is loaded/created twice, and the second time it has a different
        parent than the first time, then the second instantiation should raise an exception."""
        butler = dafPersist.Butler(outputs={'root': testInput1, 'mapper': MinMapper1})
        del butler
        butler = dafPersist.Butler(outputs={'root': testInput2, 'mapper': MinMapper1})
        del butler
        butler = dafPersist.Butler(inputs=testInput1, outputs=testOutput)
        del butler
        # should raise:
        with self.assertRaises(RuntimeError):
            butler = dafPersist.Butler(inputs=testInput2, outputs=testOutput)
            del butler

    @unittest.expectedFailure  # this is flagged to be fixed in DM-9048
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
