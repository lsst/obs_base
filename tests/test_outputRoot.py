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
import tempfile
import shutil
import pickle
import os

import lsst.utils.tests
import lsst.afw.geom as afwGeom
import lsst.pex.policy as pexPolicy
import lsst.daf.persistence as dafPersist
import lsst.obs.base

# Define paths used for testing
ROOT = os.path.abspath(os.path.dirname(__file__))


def setup_module(module):
    lsst.utils.tests.init()


class MinMapper1(lsst.obs.base.CameraMapper):
    packageName = 'larry'

    def __init__(self, **kwargs):
        policy = pexPolicy.Policy.createPolicy(os.path.join(ROOT, "MinMapper1.paf"))
        lsst.obs.base.CameraMapper.__init__(self,
                                            policy=policy, repositoryDir=ROOT, **kwargs)
        return

    def std_x(self, item, dataId):
        return float(item)

    @classmethod
    def getPackageDir(cls):
        return "/path/to/nowhere"


class OutputRootTestCase(unittest.TestCase):
    """A test case for output roots."""

    def setUp(self):
        self.tempDirs = {}

    def tearDown(self):
        for d in self.tempDirs:
            shutil.rmtree(self.tempDirs[d])

    def mkdtemp(self, prefix):
        """Create a temporary directory and return its name.
        The resulting path is stored in a dict (self.tempDirs) with
        the supplied prefix as key. This allows the name to be retrieved
        later. The directory path is returned.
        A '-' is appended to the prefix if not there."""
        dprefix = prefix
        if not dprefix.endswith('-'):
            dprefix = dprefix + '-'
        tempDir = tempfile.mkdtemp(dir=ROOT, prefix=dprefix)
        self.tempDirs[prefix] = tempDir
        return tempDir

    def testCreateOutputRoot(self):
        """Create an input repository and a related output repo, and verify there is a parent relationship
        from the output repo to the input repo."""
        testOutput = self.mkdtemp("testOutput")
        butler = dafPersist.Butler(inputs={'root': ROOT, 'mapper': MinMapper1},
                                   outputs=testOutput)
        self.assertTrue(butler)
        self.assertTrue(os.path.exists(testOutput))
        self.assertTrue(os.path.isdir(testOutput))
        self.assertTrue(os.path.exists(os.path.join(testOutput, "repositoryCfg.yaml")))
        cfg = dafPersist.PosixStorage.getRepositoryCfg(testOutput)
        expectedCfg = dafPersist.RepositoryCfg(root=ROOT,
                                               mapper=MinMapper1,
                                               mapperArgs=None,
                                               parents=None,
                                               policy=None)
        self.assertEqual(len(cfg.parents), 1)
        self.assertEqual(cfg.parents[0], expectedCfg)

    def testParentNormal(self):
        """Test that an object can be found at root location and put into an output location.
        Then test that when the output locaiton is used as an input location, and with a new output location,
        that the object is found in the first output location.
        """
        testOutput = self.mkdtemp("testOutput")
        butler = dafPersist.Butler(inputs={'root': ROOT, 'mapper': MinMapper1}, outputs=testOutput)
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

        testOutput2 = self.mkdtemp("testOutput2")
        butler = dafPersist.Butler(inputs={'root': testOutput, 'mapper': MinMapper1}, outputs=testOutput2)
        mapper2 = butler._repos.inputs()[0].repo._mapper
        loc = mapper2.map("x", dict(sensor="1,1"))
        self.assertEqual(loc.getPythonType(), "lsst.afw.geom.BoxI")
        self.assertEqual(loc.getCppType(), "BoxI")
        self.assertEqual(loc.getStorageName(), "PickleStorage")
        self.assertEqual(loc.getLocations(), ["foo-1,1.pickle"])
        self.assertEqual(loc.getStorage().root, testOutput)
        self.assertEqual(loc.getAdditionalData().toString(), "sensor = \"1,1\"\n")

    def testParentTrailingSlash2527(self):
        """Just like testParentNormal, but put a trailing slash on the root paths.
        Test that an object can be found at root location and put into an output location.
        Then test that when the output locaiton is used as an input location, and with a new output location,
        that the object is found in the first output location."""
        # todo these shouldn't be commented out, I think the test wants the trailing slash.
        testOutput = self.mkdtemp("testOutput") + '/'
        butler = dafPersist.Butler(inputs={'root': ROOT + '/', 'mapper': MinMapper1},
                                   outputs=testOutput)
        mapper1 = butler._repos.inputs()[0].repo._mapper
        loc = mapper1.map("x", dict(sensor="1,1"), write=True)
        self.assertEqual(loc.getPythonType(), "lsst.afw.geom.BoxI")
        self.assertEqual(loc.getCppType(), "BoxI")
        self.assertEqual(loc.getStorageName(), "PickleStorage")
        self.assertEqual(loc.getLocations(), ["foo-1,1.pickle"])
        self.assertEqual(loc.getStorage().root, ROOT)
        self.assertEqual(loc.getAdditionalData().toString(), "sensor = \"1,1\"\n")
        box = afwGeom.BoxI(afwGeom.PointI(0, 1), afwGeom.PointI(2, 3))
        butler.put(box, "x", sensor="1,1")
        self.assertTrue(os.path.exists(os.path.join(testOutput, loc.getLocations()[0])))
        del butler
        del mapper1

        testOutput2 = self.mkdtemp("testOutput2") + '/'
        butler = dafPersist.Butler(inputs={'root': testOutput, 'mapper': MinMapper1},
                                   outputs=testOutput2)
        mapper2 = butler._repos.inputs()[0].repo._mapper
        loc = mapper2.map("x", dict(sensor="1,1"))
        self.assertEqual(loc.getPythonType(), "lsst.afw.geom.BoxI")
        self.assertEqual(loc.getCppType(), "BoxI")
        self.assertEqual(loc.getStorageName(), "PickleStorage")
        self.assertEqual(loc.getLocations(), ["foo-1,1.pickle"])
        self.assertEqual(os.path.normpath(loc.getStorage().root),
                         os.path.normpath(testOutput))
        self.assertEqual(loc.getAdditionalData().toString(), "sensor = \"1,1\"\n")

    def testReuseOutputRoot(self):
        """Set up an output repositoriy and verify its parent relationship to the input repository.
        Then set up an output repository with the first output as an input, and verify the parent
        relationships."""
        testOutput = self.mkdtemp("testOutput")
        butler = dafPersist.Butler(inputs={'root': ROOT, 'mapper': MinMapper1},
                                   outputs=testOutput)
        self.assertTrue(os.path.exists(testOutput))
        self.assertTrue(os.path.isdir(testOutput))
        cfg = dafPersist.Storage().getRepositoryCfg(testOutput)
        expectedCfg = dafPersist.RepositoryCfg(root=ROOT,
                                               mapper=MinMapper1,
                                               mapperArgs=None,
                                               parents=None,
                                               policy=None)
        self.assertEqual(cfg.parents, [expectedCfg])
        del butler

        testOutput2 = self.mkdtemp("testOutput2")
        butler = dafPersist.Butler(inputs={'root': testOutput, 'mapper': MinMapper1},
                                   outputs=testOutput2)
        self.assertTrue(os.path.exists(testOutput2))
        self.assertTrue(os.path.isdir(testOutput2))
        cfg = dafPersist.Storage().getRepositoryCfg(testOutput2)
        self.assertEqual(cfg.parents, [testOutput])
        del butler

    def testDiffInput(self):
        """Verify that if an output repository is loaded/created twice, and the second time it has a different
        parent than the first time, then the second instantiation should raise an exception."""
        testInput1 = self.mkdtemp("testInput1")
        butler = dafPersist.Butler(outputs={'root': testInput1, 'mapper': MinMapper1})
        del butler
        testInput2 = self.mkdtemp("testInput2")
        butler = dafPersist.Butler(outputs={'root': testInput2, 'mapper': MinMapper1})
        del butler
        testOutput = self.mkdtemp("testOutput")
        butler = dafPersist.Butler(inputs=testInput1, outputs=testOutput)
        del butler
        # should raise:
        with self.assertRaises(RuntimeError):
            butler = dafPersist.Butler(inputs=testInput2, outputs=testOutput)
            del butler

    @unittest.expectedFailure  # this is flagged to be fixed in DM-9048
    def testBackup(self):
        testOutput = self.mkdtemp("testOutput")
        mapper1 = MinMapper1(outputRoot=testOutput)
        butler1 = dafPersist.Butler(outputs=dafPersist.RepositoryArgs(mode='w',
                                                                      root=testOutput,
                                                                      mapper=mapper1))
        b1 = afwGeom.Box2I(afwGeom.Point2I(3, 4), afwGeom.Point2I(7, 6), invert=False)
        butler1.put(b1, "x")
        self.assertTrue(os.path.exists(os.path.join(testOutput, "foo-1,1.pickle")))
        b2 = afwGeom.Box2I(b1)
        b2.grow(1)
        butler1.put(b2, "x", doBackup=True)
        self.assertTrue(os.path.exists(os.path.join(testOutput, "foo-1,1.pickle")))
        self.assertTrue(os.path.exists(os.path.join(testOutput, "foo-1,1.pickle~1")))
        testOutput2 = self.mkdtemp("testOutput2")
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
