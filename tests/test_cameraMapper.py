# This file is part of obs_base.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import gc
import os
import shutil
import sqlite3
import tempfile
import unittest

import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.daf.persistence as dafPersist
import lsst.geom as geom
import lsst.obs.base
import lsst.utils.tests
import numpy as np
from lsst.obs.base.test import BaseMapper

ROOT = os.path.abspath(os.path.dirname(__file__))


def setup_module(module):
    lsst.utils.tests.init()


class MinCam(lsst.obs.base.Instrument):
    @property
    def filterDefinitions(self):
        return lsst.obs.base.FilterDefinitionCollection(
            lsst.obs.base.FilterDefinition(physical_filter="u.MP9301", band="u"),
            lsst.obs.base.FilterDefinition(physical_filter="g.MP9401", band="g"),
            lsst.obs.base.FilterDefinition(physical_filter="r.MP9601", band="r", alias={"old-r"}),
            lsst.obs.base.FilterDefinition(physical_filter="i.MP9701", band="i", alias={"old-i"}),
            lsst.obs.base.FilterDefinition(physical_filter="z.MP9801", band="z"),
            # afw_name is so special-cased that only a real example will work
            lsst.obs.base.FilterDefinition(physical_filter="HSC-I2", band="i", afw_name="i2"),
        )

    @classmethod
    def getName(cls):
        return "min"

    def getCamera(self):
        raise NotImplementedError()

    def register(self, registry):
        raise NotImplementedError()

    def getRawFormatter(self, dataId):
        raise NotImplementedError()

    def makeDataIdTranslatorFactory(self):
        raise NotImplementedError()


class MinMapper1(lsst.obs.base.CameraMapper):
    packageName = "larry"

    def __init__(self, **kwargs):
        policy = dafPersist.Policy(os.path.join(ROOT, "MinMapper1.yaml"))
        lsst.obs.base.CameraMapper.__init__(self, policy=policy, repositoryDir=ROOT, **kwargs)
        return

    def std_x(self, item, dataId):
        return float(item)

    @classmethod
    def getCameraName(cls):
        """Return the name of the camera that this CameraMapper is for."""
        return "min"

    @classmethod
    def getPackageDir(cls):
        return "/path/to/nowhere"


class MinMapper2(lsst.obs.base.CameraMapper):
    packageName = "moe"
    _gen3instrument = MinCam

    # CalibRoot in policy
    # needCalibRegistry
    def __init__(self, **kwargs):
        policy = dafPersist.Policy(os.path.join(ROOT, "MinMapper2.yaml"))
        lsst.obs.base.CameraMapper.__init__(
            self, policy=policy, repositoryDir=ROOT, registry="cfhtls.sqlite3", **kwargs
        )
        return

    def _transformId(self, dataId):
        return dataId

    def _extractDetectorName(self, dataId):
        return "ccd00"

    def std_x(self, item, dataId):
        return float(item)

    @classmethod
    def getCameraName(cls):
        """Return the name of the camera that this CameraMapper is for."""
        return "min"

    @classmethod
    def getPackageDir(cls):
        return "/path/to/nowhere"


# does not assign packageName
class MinMapper3(lsst.obs.base.CameraMapper):
    def __init__(self, **kwargs):
        policy = dafPersist.Policy(os.path.join(ROOT, "MinMapper1.yaml"))
        lsst.obs.base.CameraMapper.__init__(self, policy=policy, repositoryDir=ROOT, root=ROOT)
        return

    @classmethod
    def getPackageDir(cls):
        return "/path/to/nowhere"


def checkCompression(testCase, additionalData):
    """Check that compression settings are present

    We check that we can access the required settings, and that
    the seed is non-zero (zero causes lsst.afw.math.Random to fail).
    """
    for plane in ("image", "mask", "variance"):
        for entry in (
            "compression.algorithm",
            "compression.columns",
            "compression.rows",
            "compression.quantizeLevel",
            "scaling.algorithm",
            "scaling.bitpix",
            "scaling.maskPlanes",
            "scaling.seed",
            "scaling.quantizeLevel",
            "scaling.quantizePad",
            "scaling.fuzz",
            "scaling.bscale",
            "scaling.bzero",
        ):
            additionalData.getScalar(plane + "." + entry)
        testCase.assertNotEqual(additionalData.getScalar(plane + ".scaling.seed"), 0)


class Mapper1TestCase(unittest.TestCase):
    """A test case for the mapper used by the data butler."""

    def setUp(self):
        self.mapper = MinMapper1(root=ROOT)

    def tearDown(self):
        del self.mapper

    def testGetDatasetTypes(self):
        expectedTypes = BaseMapper(ROOT).getDatasetTypes()
        #   Add the expected additional types to what the base class provides
        expectedTypes.extend(
            [
                "x",
                "x_filename",
                "badSourceHist",
                "badSourceHist_filename",
            ]
        )
        self.assertEqual(set(self.mapper.getDatasetTypes()), set(expectedTypes))

    def testMap(self):
        loc = self.mapper.map("x", {"sensor": "1,1"}, write=True)
        self.assertEqual(loc.getPythonType(), "lsst.afw.geom.BoxI")
        self.assertEqual(loc.getCppType(), "BoxI")
        self.assertEqual(loc.getStorageName(), "PickleStorage")
        expectedRoot = ROOT
        expectedLocations = ["foo-1,1.pickle"]
        self.assertEqual(loc.getStorage().root, expectedRoot)
        self.assertEqual(loc.getLocations(), expectedLocations)
        self.assertEqual(loc.getAdditionalData().toString(), 'sensor = "1,1"\n')

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
        self.assertEqual(MinMapper1.getCameraName(), "min")
        self.assertEqual(MinMapper1.getPackageName(), "larry")


class Mapper2TestCase(unittest.TestCase):
    """A test case for the mapper used by the data butler."""

    def setUp(self):
        super().setUp()
        # Force a standard set of filters even for tests that don't use
        # MinCam directly.
        MinCam()

    def testGetDatasetTypes(self):
        mapper = MinMapper2(root=ROOT)
        expectedTypes = BaseMapper(ROOT).getDatasetTypes()
        #   Add the expected additional types to what the base class provides
        expectedTypes.extend(
            [
                "flat",
                "flat_md",
                "flat_filename",
                "flat_sub",
                "raw",
                "raw_md",
                "raw_filename",
                "raw_sub",
                "some",
                "some_filename",
                "some_md",
                "some_sub",
                "someCatalog",
                "someCatalog_md",
                "someCatalog_filename",
                "someCatalog_len",
                "someCatalog_schema",
                "forced_src",
                "forced_src_md",
                "forced_src_filename",
                "forced_src_len",
                "forced_src_schema",
                "other_sub",
                "other_filename",
                "other_md",
                "other",
                "someGz",
                "someGz_filename",
                "someFz",
                "someFz_filename",
                "someGz_md",
                "someFz_sub",
                "someFz_md",
                "someGz_sub",
                "someGz_bbox",
                "someFz_bbox",
                "some_bbox",
                "other_bbox",
                "someExp",
                "someExp_filename",
                "someExp_md",
                "someExp_sub",
                "someExp_bbox",
                "someExp_filterLabel",
                "someExp_photoCalib",
                "someExp_visitInfo",
                "someExp_detector",
                "someExp_filter",
                "someExp_header_wcs",
                "someExp_wcs",
            ]
        )
        self.assertEqual(set(mapper.getDatasetTypes()), set(expectedTypes))

    def testMap(self):
        mapper = MinMapper2(root=ROOT)
        loc = mapper.map("raw", {"ccd": 13}, write=True)
        self.assertEqual(loc.getPythonType(), "lsst.afw.image.ExposureU")
        self.assertEqual(loc.getCppType(), "ImageU")
        self.assertEqual(loc.getStorageName(), "FitsStorage")
        self.assertEqual(loc.getLocations(), ["foo-13.fits"])
        self.assertEqual(loc.getStorage().root, ROOT)
        self.assertEqual(loc.getAdditionalData().getScalar("ccd"), 13)
        checkCompression(self, loc.getAdditionalData())

    def testSubMap(self):
        bbox = geom.BoxI(geom.Point2I(200, 100), geom.Extent2I(300, 400))
        mapper = MinMapper2(root=ROOT)
        loc = mapper.map("raw_sub", {"ccd": 13, "bbox": bbox}, write=True)
        self.assertEqual(loc.getPythonType(), "lsst.afw.image.ExposureU")
        self.assertEqual(loc.getCppType(), "ImageU")
        self.assertEqual(loc.getStorageName(), "FitsStorage")
        self.assertEqual(loc.getLocations(), ["foo-13.fits"])
        self.assertEqual(loc.getStorage().root, ROOT)
        self.assertEqual(loc.getAdditionalData().getScalar("ccd"), 13)
        self.assertEqual(loc.getAdditionalData().getScalar("width"), 300)
        self.assertEqual(loc.getAdditionalData().getScalar("height"), 400)
        self.assertEqual(loc.getAdditionalData().getScalar("llcX"), 200)
        self.assertEqual(loc.getAdditionalData().getScalar("llcY"), 100)
        checkCompression(self, loc.getAdditionalData())

        loc = mapper.map("raw_sub", {"ccd": 13, "bbox": bbox, "imageOrigin": "PARENT"}, write=True)
        self.assertEqual(loc.getPythonType(), "lsst.afw.image.ExposureU")
        self.assertEqual(loc.getCppType(), "ImageU")
        self.assertEqual(loc.getStorageName(), "FitsStorage")
        self.assertEqual(loc.getLocations(), ["foo-13.fits"])
        self.assertEqual(loc.getStorage().root, ROOT)
        self.assertEqual(loc.getAdditionalData().getScalar("ccd"), 13)
        self.assertEqual(loc.getAdditionalData().getScalar("width"), 300)
        self.assertEqual(loc.getAdditionalData().getScalar("height"), 400)
        self.assertEqual(loc.getAdditionalData().getScalar("llcX"), 200)
        self.assertEqual(loc.getAdditionalData().getScalar("llcY"), 100)
        self.assertEqual(loc.getAdditionalData().getScalar("imageOrigin"), "PARENT")
        checkCompression(self, loc.getAdditionalData())

    def testCatalogExtras(self):
        butler = dafPersist.Butler(root=ROOT, mapper=MinMapper2)
        schema = afwTable.Schema()
        aa = schema.addField("a", type=np.int32, doc="a")
        bb = schema.addField("b", type=np.float64, doc="b")
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
            self.assertEqual(header.getScalar("NAXIS2"), size)
        finally:
            try:
                os.remove(filename)
            except OSError as exc:
                print("Warning: could not remove file %r: %s" % (filename, exc))

    def testImage(self):
        mapper = MinMapper2(root=ROOT)
        loc = mapper.map("some", dict(ccd=35))
        expectedLocations = ["bar-35.fits"]
        self.assertEqual(loc.getStorage().root, ROOT)
        self.assertEqual(loc.getLocations(), expectedLocations)

        butler = dafPersist.ButlerFactory(mapper=mapper).create()
        image = butler.get("some", ccd=35)
        self.assertEqual(image.getFilter().bandLabel, "r")

        self.assertEqual(butler.get("some_bbox", ccd=35), image.getBBox())

        bbox = geom.BoxI(geom.Point2I(200, 100), geom.Extent2I(300, 400))
        image = butler.get("some_sub", ccd=35, bbox=bbox, imageOrigin="LOCAL", immediate=True)
        self.assertEqual(image.getHeight(), 400)
        self.assertEqual(image.getWidth(), 300)

    def testFilter(self):
        """Test that the same (patched) filter is returned through all Butler
        retrieval paths.
        """
        mapper = MinMapper2(root=ROOT)

        butler = dafPersist.ButlerFactory(mapper=mapper).create()
        image = butler.get("someExp", ccd=35)
        filter = butler.get("someExp_filter", ccd=35)
        # Test only valid with a complete filter
        self.assertEqual(image.getFilter(), afwImage.FilterLabel(band="r", physical="r.MP9601"))
        # Datasets should give consistent answers
        self.assertEqual(filter, image.getFilter())

    def testDetector(self):
        mapper = MinMapper2(root=ROOT)
        butler = dafPersist.ButlerFactory(mapper=mapper).create()
        detector = butler.get("raw_detector", ccd=0)
        self.assertEqual(detector.getName(), "ccd00")

    def testGzImage(self):
        mapper = MinMapper2(root=ROOT)
        loc = mapper.map("someGz", dict(ccd=35))
        expectedLocations = [os.path.join("gz", "bar-35.fits.gz")]
        self.assertEqual(loc.getStorage().root, ROOT)
        self.assertEqual(loc.getLocations(), expectedLocations)

        butler = dafPersist.ButlerFactory(mapper=mapper).create()
        image = butler.get("someGz", ccd=35)
        self.assertEqual(image.getFilter().bandLabel, "r")

        bbox = geom.BoxI(geom.Point2I(200, 100), geom.Extent2I(300, 400))
        image = butler.get("someGz_sub", ccd=35, bbox=bbox, imageOrigin="LOCAL", immediate=True)
        self.assertEqual(image.getHeight(), 400)
        self.assertEqual(image.getWidth(), 300)

    def testFzImage(self):
        mapper = MinMapper2(root=ROOT)
        loc = mapper.map("someFz", dict(ccd=35))
        expectedRoot = ROOT
        expectedLocations = [os.path.join("fz", "bar-35.fits.fz")]
        self.assertEqual(loc.getStorage().root, expectedRoot)
        self.assertEqual(loc.getLocations(), expectedLocations)

        butler = dafPersist.ButlerFactory(mapper=mapper).create()
        image = butler.get("someFz", ccd=35)
        self.assertEqual(image.getFilter().bandLabel, "r")

        bbox = geom.BoxI(geom.Point2I(200, 100), geom.Extent2I(300, 400))
        image = butler.get("someFz_sub", ccd=35, bbox=bbox, imageOrigin="LOCAL", immediate=True)
        self.assertEqual(image.getHeight(), 400)
        self.assertEqual(image.getWidth(), 300)

    def testButlerQueryMetadata(self):
        mapper = MinMapper2(root=ROOT)
        butler = dafPersist.ButlerFactory(mapper=mapper).create()
        kwargs = {"ccd": 35, "filter": "r", "visit": 787731, "taiObs": "2005-04-02T09:24:49.933440000"}
        self.assertEqual(butler.queryMetadata("other", "visit", **kwargs), [787731])
        self.assertEqual(
            butler.queryMetadata(
                "other",
                "visit",
                visit=kwargs["visit"],
                ccd=kwargs["ccd"],
                taiObs=kwargs["taiObs"],
                filter=kwargs["filter"],
            ),
            [787731],
        )
        # now test we get no matches if ccd is out of range
        self.assertEqual(butler.queryMetadata("raw", "ccd", ccd=36, filter="r", visit=787731), [])

    def testQueryMetadata(self):
        mapper = MinMapper2(root=ROOT)
        self.assertEqual(mapper.queryMetadata("raw", ["ccd"], None), [(x,) for x in range(36) if x != 3])

    def testStandardize(self):
        mapper = MinMapper2(root=ROOT)
        self.assertEqual(mapper.canStandardize("raw"), True)
        self.assertEqual(mapper.canStandardize("notPresent"), False)

    def testStandardizeFiltersFilterDefs(self):
        # tuples are (input, desired output)
        testLabels = [
            (None, None),
            (
                afwImage.FilterLabel(band="i", physical="i.MP9701"),
                afwImage.FilterLabel(band="i", physical="i.MP9701"),
            ),
            (afwImage.FilterLabel(band="i"), afwImage.FilterLabel(band="i", physical="i.MP9701")),
            (afwImage.FilterLabel(physical="i.MP9701"), afwImage.FilterLabel(band="i", physical="i.MP9701")),
            (
                afwImage.FilterLabel(band="i", physical="old-i"),
                afwImage.FilterLabel(band="i", physical="i.MP9701"),
            ),
            (afwImage.FilterLabel(physical="old-i"), afwImage.FilterLabel(band="i", physical="i.MP9701")),
            (afwImage.FilterLabel(physical="i2"), afwImage.FilterLabel(band="i", physical="HSC-I2")),
        ]
        testIds = [
            {"visit": 12345, "ccd": 42, "filter": f}
            for f in {
                "i",
                "i.MP9701",
                "old-i",
                "i2",
            }
        ]
        testData = []
        # Resolve special combinations where the expected output is different
        for input, corrected in testLabels:
            for dataId in testIds:
                if input is None:
                    if dataId["filter"] == "i2":
                        data = (input, dataId, afwImage.FilterLabel(band="i", physical="HSC-I2"))
                    else:
                        data = (input, dataId, afwImage.FilterLabel(band="i", physical="i.MP9701"))
                elif input == afwImage.FilterLabel(band="i") and dataId["filter"] == "i2":
                    data = (input, dataId, afwImage.FilterLabel(band="i", physical="HSC-I2"))
                elif corrected.physicalLabel == "HSC-I2" and dataId["filter"] in ("i.MP9701", "old-i"):
                    # Contradictory inputs, leave as-is
                    data = (input, dataId, input)
                elif corrected.physicalLabel == "i.MP9701" and dataId["filter"] == "i2":
                    # Contradictory inputs, leave as-is
                    data = (input, dataId, input)
                else:
                    data = (input, dataId, corrected)
                testData.append(data)

        mapper = MinMapper2(root=ROOT)
        for label, dataId, corrected in testData:
            exposure = afwImage.ExposureF()
            exposure.setFilter(label)
            mapper._setFilter(mapper.exposures["raw"], exposure, dataId)
            self.assertEqual(exposure.getFilter(), corrected, msg=f"Started from {label} and {dataId}")

    def testStandardizeFiltersFilterNoDefs(self):
        testLabels = [
            None,
            afwImage.FilterLabel(band="i", physical="i.MP9701"),
            afwImage.FilterLabel(band="i"),
            afwImage.FilterLabel(physical="i.MP9701"),
            afwImage.FilterLabel(band="i", physical="old-i"),
            afwImage.FilterLabel(physical="old-i"),
            afwImage.FilterLabel(physical="i2"),
        ]
        testIds = [
            {"visit": 12345, "ccd": 42, "filter": f}
            for f in {
                "i",
                "i.MP9701",
                "old-i",
                "i2",
            }
        ]
        testData = []
        # Resolve special combinations where the expected output is different
        for input in testLabels:
            for dataId in testIds:
                # FilterLabel is only source of truth if no FilterDefinitions.
                data = (input, dataId, input)
                testData.append(data)

        mapper = MinMapper1(root=ROOT)
        for label, dataId, corrected in testData:
            exposure = afwImage.ExposureF()
            exposure.setFilter(label)
            mapper._setFilter(mapper.exposures["raw"], exposure, dataId)
            self.assertEqual(exposure.getFilter(), corrected, msg=f"Started from {label} and {dataId}")

    def testCalib(self):
        mapper = MinMapper2(root=ROOT)
        loc = mapper.map("flat", {"visit": 787650, "ccd": 13}, write=True)
        self.assertEqual(loc.getPythonType(), "lsst.afw.image.ExposureF")
        self.assertEqual(loc.getCppType(), "ExposureF")
        self.assertEqual(loc.getStorageName(), "FitsStorage")
        expectedRoot = ROOT
        expectedLocations = ["flat-05Am03-fi.fits"]
        self.assertEqual(loc.getStorage().root, expectedRoot)
        self.assertEqual(loc.getLocations(), expectedLocations)
        self.assertEqual(loc.getAdditionalData().getScalar("ccd"), 13)
        self.assertEqual(loc.getAdditionalData().getScalar("visit"), 787650)
        self.assertEqual(loc.getAdditionalData().getScalar("derivedRunId"), "05Am03")
        self.assertEqual(loc.getAdditionalData().getScalar("filter"), "i")
        checkCompression(self, loc.getAdditionalData())

    def testNames(self):
        self.assertEqual(MinMapper2.getCameraName(), "min")
        self.assertEqual(MinMapper2.getPackageName(), "moe")

    @unittest.expectedFailure
    def testParentSearch(self):
        mapper = MinMapper2(root=ROOT)
        paths = mapper.parentSearch(
            os.path.join(ROOT, "testParentSearch"),
            os.path.join(ROOT, os.path.join("testParentSearch", "bar.fits")),
        )
        self.assertEqual(paths, [os.path.join(ROOT, os.path.join("testParentSearch", "bar.fits"))])
        paths = mapper.parentSearch(
            os.path.join(ROOT, "testParentSearch"),
            os.path.join(ROOT, os.path.join("testParentSearch", "bar.fits[1]")),
        )
        self.assertEqual(paths, [os.path.join(ROOT, os.path.join("testParentSearch", "bar.fits[1]"))])

        paths = mapper.parentSearch(
            os.path.join(ROOT, "testParentSearch"),
            os.path.join(ROOT, os.path.join("testParentSearch", "baz.fits")),
        )
        self.assertEqual(paths, [os.path.join(ROOT, os.path.join("testParentSearch", "_parent", "baz.fits"))])
        paths = mapper.parentSearch(
            os.path.join(ROOT, "testParentSearch"),
            os.path.join(ROOT, os.path.join("testParentSearch", "baz.fits[1]")),
        )
        self.assertEqual(
            paths, [os.path.join(ROOT, os.path.join("testParentSearch", "_parent", "baz.fits[1]"))]
        )

    def testSkymapLookups(self):
        """Test that metadata lookups don't try to get skymap data ID values
        from the registry.
        """
        mapper = MinMapper2(root=ROOT)
        butler = dafPersist.Butler(mapper=mapper)
        with self.assertRaises(RuntimeError) as manager:
            butler.dataRef("forced_src", visit=787650, ccd=13)
            self.assertIn("Cannot lookup skymap key 'tract'", str(manager.exception))
        # We're mostly concerned that the statements below will raise an
        # exception; if they don't, it's not likely the following tests will
        # fail.
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
    @staticmethod
    def _createRegistry(path):
        cmd = """CREATE TABLE x(
           id INT,
           visit INT,
           filter TEXT,
           snap INT,
           raft TEXT,
           sensor TEXT,
           channel TEXT,
           taiObs TEXT,
           expTime REAL
           );
        """
        conn = sqlite3.connect(path)
        conn.cursor().execute(cmd)
        conn.commit()
        conn.close()

    def setUp(self):
        self.ROOT = tempfile.mkdtemp(dir=ROOT, prefix="ParentRegistryTestCase-")
        self.repoARoot = os.path.join(self.ROOT, "a")
        args = dafPersist.RepositoryArgs(root=self.repoARoot, mapper=MinMapper1)
        butler = dafPersist.Butler(outputs=args)
        self._createRegistry(os.path.join(self.repoARoot, "registry.sqlite3"))
        del butler

    def tearDown(self):
        # the butler sql registry closes its database connection in __del__.
        # To trigger __del__ we explicitly collect the garbage here. If we
        # find having or closing the open database connection is a problem in
        # production code, we may need to add api to butler to explicity
        # release database connections (and maybe other things like in-memory
        # cached objects).
        gc.collect()
        if os.path.exists(self.ROOT):
            shutil.rmtree(self.ROOT)

    def test(self):
        """Verify that when the child repo does not have a registry it is
        assigned the registry from the parent.
        """
        repoBRoot = os.path.join(self.ROOT, "b")
        butler = dafPersist.Butler(inputs=self.repoARoot, outputs=repoBRoot)
        # This way of getting the registry from the mapping is obviously going
        # way into private members and the python lambda implementation code.
        # It is very brittle and should not be duplicated in user code
        # or any location that is not trivial to fix along with changes to the
        # CameraMapper or Mapping.
        registryA = butler._repos.inputs()[0].repo._mapper.registry
        registryB = butler._repos.outputs()[0].repo._mapper.registry
        self.assertEqual(id(registryA), id(registryB))

        self._createRegistry(os.path.join(repoBRoot, "registry.sqlite3"))
        butler = dafPersist.Butler(inputs=self.repoARoot, outputs=repoBRoot)
        # see above; don't copy this way of getting the registry.
        registryA = butler._repos.inputs()[0].repo._mapper.registry
        registryB = butler._repos.outputs()[0].repo._mapper.registry
        self.assertNotEqual(id(registryA), id(registryB))


class MissingPolicyKeyTestCase(unittest.TestCase):
    def testGetRaises(self):
        butler = dafPersist.Butler(inputs={"root": ROOT, "mapper": MinMapper1})
        # MinMapper1 does not specify a template for the raw dataset type so
        # trying to use it for get should raise
        with self.assertRaises(RuntimeError) as contextManager:
            butler.get("raw")
        # This test demonstrates and verifies that simple use of the incomplete
        # dataset type returns a helpful (I hope) error message.
        self.assertEqual(
            str(contextManager.exception),
            "Template is not defined for the raw dataset type, it must be set before it can be used.",
        )
        with self.assertRaises(RuntimeError) as contextManager:
            butler.queryMetadata("raw", "unused", {})

    def testQueryMetadataRaises(self):
        butler = dafPersist.Butler(inputs={"root": ROOT, "mapper": MinMapper1})
        # MinMapper1 does not specify a template for the raw dataset type so
        # trying to use it for queryMetadata should raise
        with self.assertRaises(RuntimeError) as contextManager:
            butler.queryMetadata("raw", "unused", {})
        # This test demonstrates and verifies that simple use of the incomplete
        # dataset type returns a helpful (I hope) error message.
        self.assertEqual(
            str(contextManager.exception),
            "Template is not defined for the raw dataset type, it must be set before it can be used.",
        )

    def testFilenameRaises(self):
        butler = dafPersist.Butler(inputs={"root": ROOT, "mapper": MinMapper1})
        # MinMapper1 does not specify a template for the raw dataset type so
        # trying to use it for <datasetType>_filename should raise
        with self.assertRaises(RuntimeError) as contextManager:
            butler.get("raw_filename")
        # This test demonstrates and verifies that simple use of the
        # incomplete dataset type returns a helpful (I hope) error message.
        self.assertEqual(
            str(contextManager.exception),
            "Template is not defined for the raw dataset type, it must be set before it can be used.",
        )

    def testWcsRaises(self):
        butler = dafPersist.Butler(inputs={"root": ROOT, "mapper": MinMapper1})
        # MinMapper1 does not specify a template for the raw dataset type so
        # trying to use it for <datasetType>_wcs should raise
        with self.assertRaises(RuntimeError) as contextManager:
            butler.get("raw_wcs")
        # This test demonstrates and verifies that simple use of the
        # incomplete dataset type returns a helpful (I hope) error message.
        self.assertEqual(
            str(contextManager.exception),
            "Template is not defined for the raw dataset type, it must be set before it can be used.",
        )

    def testConflictRaises(self):
        policy = dafPersist.Policy(os.path.join(ROOT, "ConflictMapper.yaml"))
        with self.assertRaisesRegex(ValueError, r"Duplicate mapping policy for dataset type packages"):
            mapper = lsst.obs.base.CameraMapper(policy=policy, repositoryDir=ROOT, root=ROOT)  # noqa F841


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
