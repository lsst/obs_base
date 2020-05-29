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
import unittest
import tempfile
import shutil

from typing import TYPE_CHECKING

import lsst.utils.tests

import lsst.afw.image
from lsst.afw.image import LOCAL
from lsst.geom import Box2I, Point2I

from lsst.daf.butler import Config
from lsst.daf.butler import StorageClassFactory
from lsst.daf.butler import DatasetType
from lsst.daf.butler.tests import DatasetTestHelper, makeTestRepo, addDatasetType, makeTestCollection

from lsst.obs.base.exposureAssembler import ExposureAssembler

if TYPE_CHECKING:
    from lsst.daf.butler import DatasetRef

TESTDIR = os.path.dirname(__file__)

BUTLER_CONFIG = """
storageClasses:
  ExposureCompositeF:
    inheritsFrom: ExposureF
datastore:
  # Want to check disassembly so can't use InMemory
  cls: lsst.daf.butler.datastores.posixDatastore.PosixDatastore
  formatters:
    ExposureCompositeF: lsst.obs.base.fitsExposureFormatter.FitsExposureFormatter
  composites:
    disassembled:
      ExposureCompositeF: True
"""

# Components present in the test file
COMPONENTS = {"wcs", "image", "mask", "coaddInputs", "psf", "visitInfo", "variance", "metadata", "photoCalib"}


class ButlerFitsTests(DatasetTestHelper, lsst.utils.tests.TestCase):

    @classmethod
    def setUpClass(cls):
        """Create a new butler once only."""

        cls.storageClassFactory = StorageClassFactory()

        cls.root = tempfile.mkdtemp(dir=TESTDIR)

        dataIds = {
            "instrument": ["DummyCam"],
            "physical_filter": ["d-r"],
            "visit": [42],
        }

        cls.creatorButler = makeTestRepo(cls.root, dataIds, config=Config.fromYaml(BUTLER_CONFIG))

        # Create dataset types used by the tests
        for datasetTypeName, storageClassName in (("calexp", "ExposureF"),
                                                  ("unknown", "ExposureCompositeF"),
                                                  ("testCatalog", "SourceCatalog"),
                                                  ):
            storageClass = cls.storageClassFactory.getStorageClass(storageClassName)
            addDatasetType(cls.creatorButler, datasetTypeName, set(dataIds), storageClass)

    @classmethod
    def tearDownClass(cls):
        if cls.root is not None:
            shutil.rmtree(cls.root, ignore_errors=True)

    def setUp(self):
        self.butler = makeTestCollection(self.creatorButler)

    def makeExampleCatalog(self) -> lsst.afw.table.SourceCatalog:
        catalogPath = os.path.join(TESTDIR, "data", "source_catalog.fits")
        return lsst.afw.table.SourceCatalog.readFits(catalogPath)

    def assertCatalogEqual(self, inputCatalog: lsst.afw.table.SourceCatalog,
                           outputCatalog: lsst.afw.table.SourceCatalog) -> None:
        self.assertIsInstance(outputCatalog, lsst.afw.table.SourceCatalog)
        inputTable = inputCatalog.getTable()
        inputRecord = inputCatalog[0]
        outputTable = outputCatalog.getTable()
        outputRecord = outputCatalog[0]
        self.assertEqual(inputRecord.getPsfInstFlux(), outputRecord.getPsfInstFlux())
        self.assertEqual(inputRecord.getPsfFluxFlag(), outputRecord.getPsfFluxFlag())
        self.assertEqual(inputTable.getSchema().getAliasMap().get("slot_Centroid"),
                         outputTable.getSchema().getAliasMap().get("slot_Centroid"))
        self.assertEqual(inputRecord.getCentroid(), outputRecord.getCentroid())
        self.assertFloatsAlmostEqual(
            inputRecord.getCentroidErr()[0, 0],
            outputRecord.getCentroidErr()[0, 0], rtol=1e-6)
        self.assertFloatsAlmostEqual(
            inputRecord.getCentroidErr()[1, 1],
            outputRecord.getCentroidErr()[1, 1], rtol=1e-6)
        self.assertEqual(inputTable.getSchema().getAliasMap().get("slot_Shape"),
                         outputTable.getSchema().getAliasMap().get("slot_Shape"))
        self.assertFloatsAlmostEqual(
            inputRecord.getShapeErr()[0, 0],
            outputRecord.getShapeErr()[0, 0], rtol=1e-6)
        self.assertFloatsAlmostEqual(
            inputRecord.getShapeErr()[1, 1],
            outputRecord.getShapeErr()[1, 1], rtol=1e-6)
        self.assertFloatsAlmostEqual(
            inputRecord.getShapeErr()[2, 2],
            outputRecord.getShapeErr()[2, 2], rtol=1e-6)

    def testFitsCatalog(self) -> None:
        catalog = self.makeExampleCatalog()
        dataId = {"visit": 42, "instrument": "DummyCam", "physical_filter": "d-r"}
        ref = self.butler.put(catalog, "testCatalog", dataId)
        stored = self.butler.get(ref)
        self.assertCatalogEqual(catalog, stored)

    def testExposureCompositePutGetConcrete(self) -> None:
        ref = self.runExposureCompositePutGetTest("calexp")

        uri = self.butler.getURI(ref)
        self.assertTrue(os.path.exists(uri.path), f"Checking URI {uri} existence")

    def testExposureCompositePutGetVirtual(self) -> None:
        ref = self.runExposureCompositePutGetTest("unknown")

        primary, components = self.butler.getURIs(ref)
        self.assertIsNone(primary)
        self.assertEqual(set(components), COMPONENTS)
        for compName, uri in components.items():
            self.assertTrue(os.path.exists(uri.path),
                            f"Checking URI {uri} existence for component {compName}")

    def runExposureCompositePutGetTest(self, datasetTypeName: str) -> DatasetRef:
        example = os.path.join(TESTDIR, "data", "small.fits")
        exposure = lsst.afw.image.ExposureF(example)

        dataId = {"visit": 42, "instrument": "DummyCam", "physical_filter": "d-r"}
        ref = self.butler.put(exposure, datasetTypeName, dataId)

        # Get the full thing
        composite = self.butler.get(datasetTypeName, dataId)

        # There is no assert for Exposure so just look at maskedImage
        self.assertMaskedImagesEqual(composite.maskedImage, exposure.maskedImage)

        # Helper for extracting components
        assembler = ExposureAssembler(ref.datasetType.storageClass)

        # Get each component from butler independently
        for compName in COMPONENTS:
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
                self.assertEqual(len(component.visits), len(reference.visits),
                                 f"cf visits {component.visits}")
                self.assertEqual(len(component.ccds), len(reference.ccds),
                                 f"cf CCDs {component.ccds}")
            elif compName == "psf":
                # Equality for PSF does not work
                pass
            elif compName == "visitInfo":
                self.assertEqual(component.getExposureId(), reference.getExposureId(),
                                 f"VisitInfo comparison")
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
                self.assertIn("spatially constant with mean: inf", str(component), "Checking photoCalib")
            else:
                raise RuntimeError(f"Unexpected component '{compName}' encountered in test")

        # With parameters
        inBBox = Box2I(minimum=Point2I(0, 0), maximum=Point2I(3, 3))
        parameters = dict(bbox=inBBox, origin=LOCAL)
        subset = self.butler.get(datasetTypeName, dataId, parameters=parameters)
        outBBox = subset.getBBox()
        self.assertEqual(inBBox, outBBox)
        self.assertImagesEqual(subset.getImage(), exposure.subset(inBBox, origin=LOCAL).getImage())

        return ref


if __name__ == "__main__":
    unittest.main()