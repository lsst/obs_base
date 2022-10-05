# This file is part of daf_butler.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import annotations

import os
import shutil
import tempfile
import unittest
from typing import TYPE_CHECKING

import astropy.table
import lsst.afw.cameraGeom.testUtils  # for test asserts injected into TestCase
import lsst.afw.image
import lsst.pex.config
import lsst.utils.tests
from lsst.afw.fits import readMetadata
from lsst.afw.image import LOCAL, ExposureFitsReader
from lsst.afw.math import flipImage
from lsst.daf.base import PropertyList, PropertySet
from lsst.daf.butler import Config, DatasetType, StorageClassFactory
from lsst.daf.butler.tests import DatasetTestHelper, addDatasetType, makeTestCollection, makeTestRepo
from lsst.geom import Box2I, Extent2I, Point2I
from lsst.obs.base.exposureAssembler import ExposureAssembler
from lsst.obs.base.tests import make_ramp_exposure_trimmed, make_ramp_exposure_untrimmed

if TYPE_CHECKING:
    from lsst.daf.butler import DatasetRef

TESTDIR = os.path.dirname(__file__)

BUTLER_CONFIG = """
storageClasses:
  ExposureCompositeF:
    inheritsFrom: ExposureF
datastore:
  # Want to check disassembly so can't use InMemory
  cls: lsst.daf.butler.datastores.fileDatastore.FileDatastore
  formatters:
    ExposureCompositeF: lsst.obs.base.formatters.fitsExposure.FitsExposureFormatter
    lossless:
      formatter: lsst.obs.base.formatters.fitsExposure.FitsExposureFormatter
      parameters:
        recipe: lossless
    uncompressed:
      formatter: lsst.obs.base.formatters.fitsExposure.FitsExposureFormatter
      parameters:
        recipe: noCompression
    lossy:
      formatter: lsst.obs.base.formatters.fitsExposure.FitsExposureFormatter
      parameters:
        recipe: lossyBasic
  composites:
    disassembled:
      ExposureCompositeF: True
"""

# Components present in the test file
COMPONENTS = {
    "wcs",
    "image",
    "mask",
    "coaddInputs",
    "psf",
    "visitInfo",
    "variance",
    "metadata",
    "photoCalib",
    "filter",
    "validPolygon",
    "transmissionCurve",
    "detector",
    "apCorrMap",
    "summaryStats",
    "id",
}
READ_COMPONENTS = {
    "bbox",
    "xy0",
    "dimensions",
}


class ButlerFitsTests(DatasetTestHelper, lsst.utils.tests.TestCase):
    @classmethod
    def setUpClass(cls):
        """Create a new butler once only."""

        cls.storageClassFactory = StorageClassFactory()

        cls.root = tempfile.mkdtemp(dir=TESTDIR)

        dataIds = {
            "instrument": ["DummyCam"],
            "physical_filter": ["d-r"],
            "visit": [42, 43, 44],
        }

        # Ensure that we test in a directory that will include some
        # metacharacters
        subdir = "sub?#dir"
        butlerRoot = os.path.join(cls.root, subdir)

        cls.creatorButler = makeTestRepo(butlerRoot, dataIds, config=Config.fromYaml(BUTLER_CONFIG))

        # Create dataset types used by the tests
        for datasetTypeName, storageClassName in (
            ("calexp", "ExposureF"),
            ("unknown", "ExposureCompositeF"),
            ("testCatalog", "SourceCatalog"),
            ("lossless", "ExposureF"),
            ("uncompressed", "ExposureF"),
            ("lossy", "ExposureF"),
        ):
            storageClass = cls.storageClassFactory.getStorageClass(storageClassName)
            addDatasetType(cls.creatorButler, datasetTypeName, set(dataIds), storageClass)

        # And some dataset types that have no dimensions for easy testing
        for datasetTypeName, storageClassName in (
            ("ps", "PropertySet"),
            ("pl", "PropertyList"),
            ("int_exp_trimmed", "ExposureI"),
            ("int_exp_untrimmed", "ExposureI"),
        ):
            storageClass = cls.storageClassFactory.getStorageClass(storageClassName)
            addDatasetType(cls.creatorButler, datasetTypeName, {}, storageClass)

    @classmethod
    def tearDownClass(cls):
        if cls.root is not None:
            shutil.rmtree(cls.root, ignore_errors=True)

    def setUp(self):
        self.butler = makeTestCollection(self.creatorButler)

    def makeExampleCatalog(self) -> lsst.afw.table.SourceCatalog:
        catalogPath = os.path.join(TESTDIR, "data", "source_catalog.fits")
        return lsst.afw.table.SourceCatalog.readFits(catalogPath)

    def assertCatalogEqual(
        self, inputCatalog: lsst.afw.table.SourceCatalog, outputCatalog: lsst.afw.table.SourceCatalog
    ) -> None:
        self.assertIsInstance(outputCatalog, lsst.afw.table.SourceCatalog)
        inputTable = inputCatalog.getTable()
        inputRecord = inputCatalog[0]
        outputTable = outputCatalog.getTable()
        outputRecord = outputCatalog[0]
        self.assertEqual(inputRecord.getPsfInstFlux(), outputRecord.getPsfInstFlux())
        self.assertEqual(inputRecord.getPsfFluxFlag(), outputRecord.getPsfFluxFlag())
        self.assertEqual(
            inputTable.getSchema().getAliasMap().get("slot_Centroid"),
            outputTable.getSchema().getAliasMap().get("slot_Centroid"),
        )
        self.assertEqual(inputRecord.getCentroid(), outputRecord.getCentroid())
        self.assertFloatsAlmostEqual(
            inputRecord.getCentroidErr()[0, 0], outputRecord.getCentroidErr()[0, 0], rtol=1e-6
        )
        self.assertFloatsAlmostEqual(
            inputRecord.getCentroidErr()[1, 1], outputRecord.getCentroidErr()[1, 1], rtol=1e-6
        )
        self.assertEqual(
            inputTable.getSchema().getAliasMap().get("slot_Shape"),
            outputTable.getSchema().getAliasMap().get("slot_Shape"),
        )
        self.assertFloatsAlmostEqual(
            inputRecord.getShapeErr()[0, 0], outputRecord.getShapeErr()[0, 0], rtol=1e-6
        )
        self.assertFloatsAlmostEqual(
            inputRecord.getShapeErr()[1, 1], outputRecord.getShapeErr()[1, 1], rtol=1e-6
        )
        self.assertFloatsAlmostEqual(
            inputRecord.getShapeErr()[2, 2], outputRecord.getShapeErr()[2, 2], rtol=1e-6
        )

    def runFundamentalTypeTest(self, datasetTypeName, entity):
        """Put and get the supplied entity and compare."""
        ref = self.butler.put(entity, datasetTypeName)
        butler_ps = self.butler.get(ref)
        self.assertEqual(butler_ps, entity)

        # Break the contact by ensuring that we are writing YAML
        uri = self.butler.getURI(ref)
        self.assertEqual(uri.getExtension(), ".yaml", f"Check extension of {uri}")

    def testFundamentalTypes(self) -> None:
        """Ensure that some fundamental stack types round trip."""
        ps = PropertySet()
        ps["a.b"] = 5
        ps["c.d.e"] = "string"
        self.runFundamentalTypeTest("ps", ps)

        pl = PropertyList()
        pl["A"] = 1
        pl.setComment("A", "An int comment")
        pl["B"] = "string"
        pl.setComment("B", "A string comment")
        self.runFundamentalTypeTest("pl", pl)

    def testFitsCatalog(self) -> None:
        """Test reading of a FITS catalog"""
        catalog = self.makeExampleCatalog()
        dataId = {"visit": 42, "instrument": "DummyCam", "physical_filter": "d-r"}
        ref = self.butler.put(catalog, "testCatalog", dataId)
        stored = self.butler.get(ref)
        self.assertCatalogEqual(catalog, stored)

        # Override the storage class.
        astropy_table = self.butler.getDirect(ref, storageClass="AstropyTable")
        self.assertIsInstance(astropy_table, astropy.table.Table)
        self.assertEqual(len(astropy_table), len(stored))

    def testExposureCompositePutGetConcrete(self) -> None:
        """Test composite with no disassembly"""
        ref = self.runExposureCompositePutGetTest("calexp")

        uri = self.butler.getURI(ref)
        self.assertTrue(uri.exists(), f"Checking URI {uri} existence")

    def testExposureCompositePutGetVirtual(self) -> None:
        """Testing composite disassembly"""
        ref = self.runExposureCompositePutGetTest("unknown")

        primary, components = self.butler.getURIs(ref)
        self.assertIsNone(primary)
        self.assertEqual(set(components), COMPONENTS)
        for compName, uri in components.items():
            self.assertTrue(uri.exists(), f"Checking URI {uri} existence for component {compName}")

    def runExposureCompositePutGetTest(self, datasetTypeName: str) -> DatasetRef:
        example = os.path.join(TESTDIR, "data", "calexp.fits")
        exposure = lsst.afw.image.ExposureF(example)

        dataId = {"visit": 42, "instrument": "DummyCam", "physical_filter": "d-r"}
        ref = self.butler.put(exposure, datasetTypeName, dataId)

        # Get the full thing
        composite = self.butler.get(datasetTypeName, dataId)

        # There is no assert for Exposure so just look at maskedImage
        self.assertMaskedImagesEqual(composite.maskedImage, exposure.maskedImage)

        # Helper for extracting components
        assembler = ExposureAssembler(ref.datasetType.storageClass)

        # Check all possible components that can be read
        allComponents = set()
        allComponents.update(COMPONENTS, READ_COMPONENTS)

        # Get each component from butler independently
        for compName in allComponents:
            compTypeName = DatasetType.nameWithComponent(datasetTypeName, compName)
            component = self.butler.get(compTypeName, dataId)

            reference = assembler.getComponent(exposure, compName)

            self.assertIsInstance(component, type(reference), f"Checking type of component {compName}")

            if compName in ("image", "variance"):
                self.assertImagesEqual(component, reference)
            elif compName == "mask":
                self.assertMasksEqual(component, reference)
            elif compName == "wcs":
                self.assertWcsAlmostEqualOverBBox(component, reference, exposure.getBBox())
            elif compName == "coaddInputs":
                self.assertEqual(
                    len(component.visits), len(reference.visits), f"cf visits {component.visits}"
                )
                self.assertEqual(len(component.ccds), len(reference.ccds), f"cf CCDs {component.ccds}")
            elif compName == "psf":
                # Equality for PSF does not work
                pass
            elif compName == "filter":
                self.assertEqual(component, reference)
            elif compName == "id":
                self.assertEqual(component, reference)
            elif compName == "visitInfo":
                self.assertEqual(component, reference, "VisitInfo comparison")
            elif compName == "metadata":
                # The component metadata has extra fields in it so cannot
                # compare directly.
                for k, v in reference.items():
                    self.assertEqual(component[k], v)
            elif compName == "photoCalib":
                # This example has a
                # "spatially constant with mean: inf error: nan" entry
                # which does not compare directly.
                self.assertEqual(str(component), str(reference))
                self.assertIn("spatially constant with mean: 1.99409", str(component), "Checking photoCalib")
            elif compName in ("bbox", "xy0", "dimensions", "validPolygon"):
                self.assertEqual(component, reference)
            elif compName == "apCorrMap":
                self.assertEqual(set(component.keys()), set(reference.keys()))
            elif compName == "transmissionCurve":
                self.assertEqual(component.getThroughputAtBounds(), reference.getThroughputAtBounds())
            elif compName == "detector":
                c_amps = {a.getName() for a in component.getAmplifiers()}
                r_amps = {a.getName() for a in reference.getAmplifiers()}
                self.assertEqual(c_amps, r_amps)
            elif compName == "summaryStats":
                self.assertEqual(component.psfSigma, reference.psfSigma)
            else:
                raise RuntimeError(f"Unexpected component '{compName}' encountered in test")

        # Full Exposure with parameters
        inBBox = Box2I(minimum=Point2I(3, 3), maximum=Point2I(21, 16))
        parameters = dict(bbox=inBBox, origin=LOCAL)
        subset = self.butler.get(datasetTypeName, dataId, parameters=parameters)
        outBBox = subset.getBBox()
        self.assertEqual(inBBox, outBBox)
        self.assertImagesEqual(subset.getImage(), exposure.subset(inBBox, origin=LOCAL).getImage())

        return ref

    def putFits(self, exposure, datasetTypeName, visit):
        """Put different datasetTypes and return information."""
        dataId = {"visit": visit, "instrument": "DummyCam", "physical_filter": "d-r"}
        refC = self.butler.put(exposure, datasetTypeName, dataId)
        uriC = self.butler.getURI(refC)
        stat = os.stat(uriC.ospath)
        size = stat.st_size
        # We can't use butler's Exposure storage class metadata component,
        # because that intentionally strips keywords that science code
        # shouldn't ever read (to at least weaken our assumptions that we write
        # to FITS).  Instead use a lower-level method on the URI for this test.
        meta = readMetadata(uriC.ospath)
        return meta, size

    def testCompression(self):
        """Test that we can write compressed and uncompressed FITS."""
        example = os.path.join(TESTDIR, "data", "small.fits")
        exposure = lsst.afw.image.ExposureF(example)

        # Write a lossless compressed
        metaC, sizeC = self.putFits(exposure, "lossless", 42)
        self.assertEqual(metaC["TTYPE1"], "COMPRESSED_DATA")
        self.assertEqual(metaC["ZCMPTYPE"], "GZIP_2")

        # Write an uncompressed FITS file
        metaN, sizeN = self.putFits(exposure, "uncompressed", 43)
        self.assertNotIn("ZCMPTYPE", metaN)

        # Write an uncompressed FITS file
        metaL, sizeL = self.putFits(exposure, "lossy", 44)
        self.assertEqual(metaL["TTYPE1"], "COMPRESSED_DATA")
        self.assertEqual(metaL["ZCMPTYPE"], "RICE_1")

        self.assertNotEqual(sizeC, sizeN)
        # Data file is so small that Lossy and Compressed are dominated
        # by the extra compression tables
        self.assertEqual(sizeL, sizeC)

    def testExposureFormatterAmpParameter(self):
        """Test the FitsExposureFormatter implementation of the Exposure
        StorageClass's 'amp' and 'detector' parameters.
        """
        # Our example exposure file has a realistic detector (looks like HSC),
        # but the image itself doesn't match it.  So we just load the detector,
        # and use it to make our own more useful images and put them.
        detector = ExposureFitsReader(os.path.join(TESTDIR, "data", "calexp.fits")).readDetector()
        trimmed_full = make_ramp_exposure_trimmed(detector)
        untrimmed_full = make_ramp_exposure_untrimmed(detector)
        trimmed_ref = self.butler.put(trimmed_full, "int_exp_trimmed")
        untrimmed_ref = self.butler.put(untrimmed_full, "int_exp_untrimmed")
        for n, amp in enumerate(detector):
            # Try to read each amp as it is on disk, with a variety of
            # parameters that should all do the same thing.
            for amp_parameter in [amp, amp.getName(), n]:
                for parameters in [{"amp": amp_parameter}, {"amp": amp_parameter, "detector": detector}]:
                    with self.subTest(parameters=parameters):
                        test_trimmed = self.butler.getDirect(trimmed_ref, parameters=parameters)
                        test_untrimmed = self.butler.getDirect(untrimmed_ref, parameters=parameters)
                        self.assertImagesEqual(test_trimmed.image, trimmed_full[amp.getBBox()].image)
                        self.assertImagesEqual(test_untrimmed.image, untrimmed_full[amp.getRawBBox()].image)
                        self.assertEqual(len(test_trimmed.getDetector()), 1)
                        self.assertEqual(len(test_untrimmed.getDetector()), 1)
                        self.assertAmplifiersEqual(test_trimmed.getDetector()[0], amp)
                        self.assertAmplifiersEqual(test_untrimmed.getDetector()[0], amp)
            # Try to read various transformed versions of the original amp,
            # to make sure flips and offsets are applied correctly.
            # First flip X only.
            amp_t1 = amp.rebuild().transform(outFlipX=True).finish()
            test_t1_trimmed = self.butler.getDirect(trimmed_ref, parameters={"amp": amp_t1})
            self.assertImagesEqual(
                test_t1_trimmed.image, flipImage(trimmed_full[amp.getBBox()].image, flipLR=True, flipTB=False)
            )
            self.assertAmplifiersEqual(test_t1_trimmed.getDetector()[0], amp_t1)
            test_t1_untrimmed = self.butler.getDirect(untrimmed_ref, parameters={"amp": amp_t1})
            self.assertImagesEqual(
                test_t1_untrimmed.image,
                flipImage(untrimmed_full[amp.getRawBBox()].image, flipLR=True, flipTB=False),
            )
            self.assertAmplifiersEqual(test_t1_trimmed.getDetector()[0], amp_t1)
            # Flip Y only.
            amp_t2 = amp.rebuild().transform(outFlipY=True).finish()
            test_t2_trimmed = self.butler.getDirect(trimmed_ref, parameters={"amp": amp_t2})
            self.assertImagesEqual(
                test_t2_trimmed.image, flipImage(trimmed_full[amp.getBBox()].image, flipLR=False, flipTB=True)
            )
            self.assertAmplifiersEqual(test_t2_trimmed.getDetector()[0], amp_t2)
            test_t2_untrimmed = self.butler.getDirect(untrimmed_ref, parameters={"amp": amp_t2})
            self.assertImagesEqual(
                test_t2_untrimmed.image,
                flipImage(untrimmed_full[amp.getRawBBox()].image, flipLR=False, flipTB=True),
            )
            self.assertAmplifiersEqual(test_t2_trimmed.getDetector()[0], amp_t2)
            # Add an XY offset only.
            amp_t3 = amp.rebuild().transform(outOffset=Extent2I(5, 4)).finish()
            test_t3_trimmed = self.butler.getDirect(trimmed_ref, parameters={"amp": amp_t3})
            self.assertImagesEqual(test_t3_trimmed.image, trimmed_full[amp.getBBox()].image)
            self.assertAmplifiersEqual(test_t3_trimmed.getDetector()[0], amp_t3)
            test_t3_untrimmed = self.butler.getDirect(untrimmed_ref, parameters={"amp": amp_t3})
            self.assertImagesEqual(test_t3_untrimmed.image, untrimmed_full[amp.getRawBBox()].image)
            self.assertAmplifiersEqual(test_t3_trimmed.getDetector()[0], amp_t3)


if __name__ == "__main__":
    unittest.main()
