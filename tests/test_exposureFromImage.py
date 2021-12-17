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

import unittest

import lsst.afw.image as afwImage
import lsst.utils.tests
from lsst.daf.base import PropertyList
from lsst.obs.base import exposureFromImage


class ExposureFromImageTestCase(lsst.utils.tests.TestCase):
    """A test case for exposureFromImage."""

    def setUp(self):
        self.maskedImage = makeRampMaskedImage(10, 11)

    def tearDown(self):
        del self.maskedImage

    def testDecoratedImage(self):
        image = self.maskedImage.getImage()
        decoImage = afwImage.DecoratedImageF(image)
        metadata = PropertyList()
        metadata.set("FOO", "BAR")
        decoImage.setMetadata(metadata)
        exposure = exposureFromImage(decoImage)
        self.assertImagesEqual(exposure.getMaskedImage().getImage(), image)
        md = exposure.getMetadata()
        self.assertEqual(md.getScalar("FOO"), "BAR")

    def testExposure(self):
        inExposure = afwImage.ExposureF(self.maskedImage)
        outExposure = exposureFromImage(inExposure)
        self.assertIs(inExposure, outExposure)

    def testImage(self):
        image = self.maskedImage.getImage()
        exposure = exposureFromImage(image)
        self.assertImagesEqual(image, exposure.getMaskedImage().getImage())

    def testMaskedImage(self):
        exposure = exposureFromImage(self.maskedImage)
        self.assertMaskedImagesEqual(self.maskedImage, exposure.getMaskedImage())

    def testDecoratedImageBadWcs(self):
        """Test that exposureFromImage() attaches a None wcs to the exposure
        when the WCS cannot be constructed
        """
        image = self.maskedImage.getImage()
        decoImage = afwImage.DecoratedImageF(image)
        metadata = PropertyList()
        metadata.set("CTYPE1", "RA---TPV")
        metadata.set("CTYPE2", "DEC--TPV")
        decoImage.setMetadata(metadata)
        exposure = exposureFromImage(decoImage)
        self.assertIs(exposure.getWcs(), None)


def makeRampMaskedImage(width, height, imgClass=afwImage.MaskedImageF):
    """Make a ramp image of the specified size and image class

    Image values start from 0 at the lower left corner and increase by 1 along
    rows. Variance values equal image values + 100.
    Mask values equal image values modulo 8 bits (leaving plenty of unused
    values).
    """
    mi = imgClass(width, height)
    image = mi.getImage()
    mask = mi.getMask()
    variance = mi.getVariance()
    val = 0
    for yInd in range(height):
        for xInd in range(width):
            image[xInd, yInd, afwImage.LOCAL] = val
            variance[xInd, yInd, afwImage.LOCAL] = val + 100
            mask[xInd, yInd, afwImage.LOCAL] = val % 0x100
            val += 1
    return mi


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
