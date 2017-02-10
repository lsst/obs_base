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

from builtins import range
import collections
import os
import unittest

import lsst.utils.tests
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.daf.persistence as dafPersist
import lsst.obs.base
from lsst.utils import getPackageDir
import shutil

testDir = os.path.relpath(os.path.join(getPackageDir('obs_base'), 'tests'))


def setup_module(module):
    lsst.utils.tests.init()


class BaseMapper(lsst.obs.base.CameraMapper):
    packageName = 'base'

    def __init__(self):
        policy = dafPersist.Policy(os.path.join(testDir, "BaseMapper.paf"))
        lsst.obs.base.CameraMapper.__init__(self, policy=policy, repositoryDir=testDir, root=testDir)
        return


class MinMapper1(lsst.obs.base.CameraMapper):
    packageName = 'larry'

    def __init__(self, **kwargs):
        policy = dafPersist.Policy(os.path.join(testDir, "MinMapper1.paf"))
        lsst.obs.base.CameraMapper.__init__(self, policy=policy, repositoryDir=testDir, **kwargs)
        return

    def std_x(self, item, dataId):
        return float(item)

    @classmethod
    def getCameraName(cls):
        """Return the name of the camera that this CameraMapper is for."""
        return "min"


class MinMapper2(lsst.obs.base.CameraMapper):
    packageName = 'moe'

    # CalibRoot in policy
    # needCalibRegistry
    def __init__(self, **kwargs):
        policy = dafPersist.Policy(os.path.join(testDir, "MinMapper2.paf"))
        lsst.obs.base.CameraMapper.__init__(self, policy=policy, repositoryDir=testDir,
                                            registry=os.path.join(testDir, "cfhtls.sqlite3"), **kwargs)
        return

    def _transformId(self, dataId):
        return dataId

    def _extractDetectorName(self, dataId):
        return "Detector"

    def std_x(self, item, dataId):
        return float(item)

    @classmethod
    def getCameraName(cls):
        """Return the name of the camera that this CameraMapper is for."""
        return "min"


# does not assign packageName
class MinMapper3(lsst.obs.base.CameraMapper):

    def __init__(self, **kwargs):
        policy = dafPersist.Policy(os.path.join(testDir, "MinMapper1.paf"))
        lsst.obs.base.CameraMapper.__init__(self, policy=policy, repositoryDir=testDir, root=testDir)
        return


class Mapper1TestCase(unittest.TestCase):
    """A test case for the mapper used by the data butler."""

    def setUp(self):
        self.mapper = MinMapper1(root=testDir)

    def tearDown(self):
        del self.mapper

    def testGetDatasetTypes(self):
        expectedTypes = BaseMapper().getDatasetTypes()
        #   Add the expected additional types to what the base class provides
        expectedTypes.extend(["x", "x_filename",
                              "badSourceHist", "badSourceHist_filename", ])
        self.assertEqual(set(self.mapper.getDatasetTypes()), set(expectedTypes))

    def testMap(self):
        loc = self.mapper.map("x", {"sensor": "1,1"}, write=True)
        self.assertEqual(loc.getPythonType(), "lsst.afw.geom.BoxI")
        self.assertEqual(loc.getCppType(), "BoxI")
        self.assertEqual(loc.getStorageName(), "PickleStorage")
        expectedRoot = testDir
        expectedLocations = ["foo-1,1.pickle"]
        self.assertEqual(loc.getStorage().getRoot(), expectedRoot)
        self.assertEqual(loc.getLocations(), expectedLocations)
        self.assertEqual(loc.getAdditionalData().toString(),
                         "sensor = \"1,1\"\n")

    def testQueryMetadata(self):
        self.assertEqual(self.mapper.queryMetadata("x", ["sensor"], None), [("1,1",)])

    def testStandardize(self):
        self.assertTrue(self.mapper.canStandardize("x"))
        self.assertFalse(self.mapper.canStandardize("badSourceHist"))
        self.assertFalse(self.mapper.canStandardize("notPresent"))
        result = self.mapper.standardize("x", 3, None)
        self.assertIsInstance(result, float)
        self.assertEqual(result, 3.0)
        result = self.mapper.standardize("x", 3.14, None)
        self.assertIsInstance(result, float)
        self.assertEqual(result, 3.14)
        result = self.mapper.standardize("x", "3.14", None)
        self.assertIsInstance(result, float)
        self.assertEqual(result, 3.14)

    def testNames(self):
        name = MinMapper1.getCameraName()
        self.assertEqual(MinMapper1.getCameraName(), "min")
        self.assertEqual(MinMapper1.getPackageName(), "larry")


class Mapper2TestCase(unittest.TestCase):
    """A test case for the mapper used by the data butler."""

    def testGetDatasetTypes(self):
        mapper = MinMapper2(root=testDir)
        expectedTypes = BaseMapper().getDatasetTypes()
        #   Add the expected additional types to what the base class provides
        expectedTypes.extend(["flat", "flat_md", "flat_filename", "flat_sub",
                              "raw", "raw_md", "raw_filename", "raw_sub",
                              "some", "some_filename", "some_md", "some_sub",
                              "someCatalog", "someCatalog_md", "someCatalog_filename",
                              "someCatalog_len", "someCatalog_schema",
                              "forced_src", "forced_src_md", "forced_src_filename",
                              "forced_src_len","forced_src_schema",
                              "other_sub", "other_filename", "other_md", "other",
                              "someGz", "someGz_filename", "someFz", "someFz_filename", "someGz_md",
                              "someFz_sub", "someFz_md", "someGz_sub"
                              ])
        self.assertEqual(set(mapper.getDatasetTypes()),
                         set(expectedTypes))

    def testMap(self):
        mapper = MinMapper2(root=testDir)
        loc = mapper.map("raw", {"ccd": 13}, write=True)
        self.assertEqual(loc.getPythonType(), "lsst.afw.image.ExposureU")
        self.assertEqual(loc.getCppType(), "ImageU")
        self.assertEqual(loc.getStorageName(), "FitsStorage")
        self.assertEqual(loc.getLocations(), ["foo-13.fits"])
        self.assertEqual(loc.getStorage().getRoot(), testDir)
        self.assertEqual(loc.getAdditionalData().toString(), "ccd = 13\n")

    def testSubMap(self):
        if hasattr(afwGeom, 'makePointI'):
            # old afw (pre-#1556) interface
            bbox = afwGeom.BoxI(afwGeom.makePointI(200, 100),
                                afwGeom.makeExtentI(300, 400))
        else:
            # new afw (post-#1556) interface
            bbox = afwGeom.BoxI(afwGeom.Point2I(200, 100),
                                afwGeom.Extent2I(300, 400))
        mapper = MinMapper2(root=testDir)
        loc = mapper.map("raw_sub", {"ccd": 13, "bbox": bbox}, write=True)
        self.assertEqual(loc.getPythonType(), "lsst.afw.image.ExposureU")
        self.assertEqual(loc.getCppType(), "ImageU")
        self.assertEqual(loc.getStorageName(), "FitsStorage")
        self.assertEqual(loc.getLocations(), ["foo-13.fits"])
        self.assertEqual(loc.getStorage().getRoot(), testDir)
        self.assertEqual(loc.getAdditionalData().toString(),
                         'ccd = 13\nheight = 400\nllcX = 200\nllcY = 100\nwidth = 300\n')

        loc = mapper.map("raw_sub", {"ccd": 13, "bbox": bbox, "imageOrigin": "PARENT"}, write=True)
        self.assertEqual(loc.getPythonType(), "lsst.afw.image.ExposureU")
        self.assertEqual(loc.getCppType(), "ImageU")
        self.assertEqual(loc.getStorageName(), "FitsStorage")
        self.assertEqual(loc.getLocations(), ["foo-13.fits"])
        self.assertEqual(loc.getStorage().getRoot(), testDir)
        self.assertEqual(loc.getAdditionalData().toString(),
                         'ccd = 13\nheight = 400\nimageOrigin = "PARENT"\n'
                         'llcX = 200\nllcY = 100\nwidth = 300\n')

    def testCatalogExtras(self):
        butler = dafPersist.Butler(root=testDir, mapper=MinMapper2)
        schema = afwTable.Schema()
        aa = schema.addField("a", type=int, doc="a")
        bb = schema.addField("b", type=float, doc="b")
        catalog = lsst.afw.table.BaseCatalog(schema)
        row = catalog.addNew()
        row.set(aa, 12345)
        row.set(bb, 1.2345)
        size = len(catalog)
        dataId = dict(visit=123, ccd=45)
        butler.put(catalog, "someCatalog", dataId)
        filename = butler.get("someCatalog_filename", dataId)[0]
        try:
            self.assertTrue(os.path.exists(filename))
            self.assertEqual(butler.get("someCatalog_schema", dataId), schema)
            self.assertEqual(butler.get("someCatalog_len", dataId), size)
            header = butler.get("someCatalog_md", dataId)
            self.assertEqual(header.get("NAXIS2"), size)
        finally:
            try:
                os.remove(filename)
            except OSError as exc:
                print("Warning: could not remove file %r: %s" % (filename, exc))

    def testImage(self):
        mapper = MinMapper2(root=testDir)
        loc = mapper.map("some", dict(ccd=35))
        expectedLocations = ["bar-35.fits"]
        self.assertEqual(loc.getStorage().getRoot(), testDir)
        self.assertEqual(loc.getLocations(), expectedLocations)

        butler = dafPersist.ButlerFactory(mapper=mapper).create()
        image = butler.get("some", ccd=35)
        self.assertEqual(image.getFilter().getName(), "r")

        bbox = afwGeom.BoxI(afwGeom.Point2I(200, 100),
                            afwGeom.Extent2I(300, 400))
        image = butler.get("some_sub", ccd=35, bbox=bbox, imageOrigin="LOCAL", immediate=True)
        self.assertEqual(image.getHeight(), 400)
        self.assertEqual(image.getWidth(), 300)

    def testGzImage(self):
        mapper = MinMapper2(root=testDir)
        loc = mapper.map("someGz", dict(ccd=35))
        expectedLocations = [os.path.join("gz", "bar-35.fits.gz")]
        self.assertEqual(loc.getStorage().getRoot(), testDir)
        self.assertEqual(loc.getLocations(), expectedLocations)

        butler = dafPersist.ButlerFactory(mapper=mapper).create()
        image = butler.get("someGz", ccd=35)
        self.assertEqual(image.getFilter().getName(), "r")

        bbox = afwGeom.BoxI(afwGeom.Point2I(200, 100),
                            afwGeom.Extent2I(300, 400))
        image = butler.get("someGz_sub", ccd=35, bbox=bbox, imageOrigin="LOCAL", immediate=True)
        self.assertEqual(image.getHeight(), 400)
        self.assertEqual(image.getWidth(), 300)

    def testFzImage(self):
        mapper = MinMapper2(root=testDir)
        loc = mapper.map("someFz", dict(ccd=35))
        expectedRoot = testDir
        expectedLocations = [os.path.join("fz", "bar-35.fits.fz")]
        self.assertEqual(loc.getStorage().getRoot(), expectedRoot)
        self.assertEqual(loc.getLocations(), expectedLocations)

        butler = dafPersist.ButlerFactory(mapper=mapper).create()
        image = butler.get("someFz", ccd=35)
        self.assertEqual(image.getFilter().getName(), "r")

        bbox = afwGeom.BoxI(afwGeom.Point2I(200, 100),
                            afwGeom.Extent2I(300, 400))
        image = butler.get("someFz_sub", ccd=35, bbox=bbox, imageOrigin="LOCAL", immediate=True)
        self.assertEqual(image.getHeight(), 400)
        self.assertEqual(image.getWidth(), 300)

    def testButlerQueryMetadata(self):
        mapper = MinMapper2(root=testDir)
        butler = dafPersist.ButlerFactory(mapper=mapper).create()
        kwargs = {"ccd": 35, "filter": "r", "visit": 787731,
                  "taiObs": "2005-04-02T09:24:49.933440000"}
        self.assertEqual(butler.queryMetadata("other", "visit", **kwargs), [787731])
        self.assertEqual(butler.queryMetadata("other", "visit",
                                              visit=kwargs["visit"], ccd=kwargs["ccd"],
                                              taiObs=kwargs["taiObs"], filter=kwargs["filter"]),
                         [787731])
        # now test we get no matches if ccd is out of range
        self.assertEqual(butler.queryMetadata("raw", "ccd", ccd=36, filter="r", visit=787731), [])

    def testQueryMetadata(self):
        mapper = MinMapper2(root=testDir)
        self.assertEqual(mapper.queryMetadata("raw", ["ccd"], None),
                         [(x,) for x in range(36) if x != 3])

    def testStandardize(self):
        mapper = MinMapper2(root=testDir)
        self.assertEqual(mapper.canStandardize("raw"), True)
        self.assertEqual(mapper.canStandardize("notPresent"), False)

    def testCalib(self):
        mapper = MinMapper2(root=testDir)
        loc = mapper.map("flat", {"visit": 787650, "ccd": 13}, write=True)
        self.assertEqual(loc.getPythonType(), "lsst.afw.image.ExposureF")
        self.assertEqual(loc.getCppType(), "ExposureF")
        self.assertEqual(loc.getStorageName(), "FitsStorage")
        expectedRoot = testDir
        expectedLocations = ["flat-05Am03-fi.fits"]
        self.assertEqual(loc.getStorage().getRoot(), expectedRoot)
        self.assertEqual(loc.getLocations(), expectedLocations)
        self.assertEqual(loc.getAdditionalData().toString(),
                         'ccd = 13\nderivedRunId = "05Am03"\nfilter = "i"\nvisit = 787650\n')

    def testNames(self):
        self.assertEqual(MinMapper2.getCameraName(), "min")
        self.assertEqual(MinMapper2.getPackageName(), "moe")

    @unittest.expectedFailure
    def testParentSearch(self):
        mapper = MinMapper2(root=testDir)
        paths = mapper.parentSearch(os.path.join(testDir, 'testParentSearch'),
                                    os.path.join(testDir, os.path.join('testParentSearch', 'bar.fits')))
        self.assertEqual(paths, [os.path.join(testDir, os.path.join('testParentSearch', 'bar.fits'))])
        paths = mapper.parentSearch(os.path.join(testDir, 'testParentSearch'),
                                    os.path.join(testDir,
                                                 os.path.join('testParentSearch', 'bar.fits[1]')))
        self.assertEqual(paths, [os.path.join(testDir, os.path.join('testParentSearch', 'bar.fits[1]'))])

        paths = mapper.parentSearch(os.path.join(testDir, 'testParentSearch'),
                                    os.path.join(testDir, os.path.join('testParentSearch', 'baz.fits')))
        self.assertEqual(paths, [os.path.join(testDir,
                                              os.path.join('testParentSearch', '_parent', 'baz.fits'))])
        paths = mapper.parentSearch(os.path.join(testDir, 'testParentSearch'),
                                    os.path.join(testDir,
                                                 os.path.join('testParentSearch', 'baz.fits[1]')))
        self.assertEqual(paths, [os.path.join(testDir,
                                              os.path.join('testParentSearch', '_parent', 'baz.fits[1]'))])

    def testSkymapLookups(self):
        """Test that metadata lookups don't try to get skymap data ID values from the registry.
        """
        mapper = MinMapper2(root=testDir)
        butler = dafPersist.Butler(mapper=mapper)
        with self.assertRaises(RuntimeError) as manager:
            butler.dataRef("forced_src", visit=787650, ccd=13)
            self.assertIn("Cannot lookup skymap key 'tract'", str(manager.exception))
        # We're mostly concerned that the statements below will raise an exception;
        # if they don't, it's not likely the following tests will fail.
        subset = butler.subset("forced_src", visit=787650, ccd=13, tract=0)
        self.assertEqual(len(subset), 1)
        dataRef = butler.dataRef("forced_src", visit=787650, ccd=13, tract=0)
        self.assertFalse(dataRef.datasetExists("forced_src"))


class Mapper3TestCase(unittest.TestCase):
    """A test case for a mapper subclass which does not assign packageName."""

    def testPackageName(self):
        with self.assertRaises(ValueError):
            MinMapper3()
        with self.assertRaises(ValueError):
            MinMapper3.getPackageName()

class ParentRegistryTestCase(unittest.TestCase):

    def setUp(self):
        self.testDir = os.path.join(testDir, "ParentRegistryTestCase")
        self.tearDown()
        self.repoARoot = os.path.join(self.testDir, 'a')
        args = dafPersist.RepositoryArgs(root=self.repoARoot, mapper=MinMapper1)
        butler = dafPersist.Butler(outputs=args)
        with open(os.path.join(self.repoARoot, 'registry.sqlite3'), 'w') as f:
            f.write('123')

    def tearDown(self):
        if os.path.exists(self.testDir):
            shutil.rmtree(self.testDir)

    def test(self):
        """Verify that when the child repo does not have a registry it is assigned the registry from the
        parent."""
        repoBRoot = os.path.join(self.testDir, 'b')
        butler = dafPersist.Butler(inputs=self.repoARoot, outputs=repoBRoot)
        # This way of getting the registry from the mapping is obviously going way into private members and
        # the python lambda implementation code. It is very brittle and should not be duplicated in user code
        # or any location that is not trivial to fix along with changes to the CameraMapper or Mapping.
        registryA = butler._repos.inputs()[0].repo._mapper.registry
        registryB = butler._repos.outputs()[0].repo._mapper.registry
        self.assertEqual(id(registryA), id(registryB))
        del butler

        with open(os.path.join(repoBRoot, 'registry.sqlite3'), 'w') as f:
            f.write('123')
        butler = dafPersist.Butler(inputs=self.repoARoot, outputs=repoBRoot)
        # see above; don't copy this way of getting the registry.
        registryA = butler._repos.inputs()[0].repo._mapper.registry
        registryB = butler._repos.outputs()[0].repo._mapper.registry
        self.assertNotEqual(id(registryA), id(registryB))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == '__main__':
    lsst.utils.tests.init()
    unittest.main()
