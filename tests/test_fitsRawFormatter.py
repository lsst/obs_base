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

import astropy.units as u
from astro_metadata_translator import FitsTranslator, StubTranslator
from astro_metadata_translator.translators.helpers import tracking_from_degree_headers
from astropy.coordinates import Angle

import lsst.afw.geom
import lsst.afw.math
import lsst.daf.base
import lsst.daf.butler
import lsst.geom
import lsst.resources
import lsst.utils.tests
from lsst.afw.cameraGeom import makeUpdatedDetector
from lsst.afw.cameraGeom.testUtils import CameraWrapper, DetectorWrapper
from lsst.obs.base import (
    FilterDefinition,
    FilterDefinitionCollection,
    FitsRawFormatterBase,
    MakeRawVisitInfoViaObsInfo,
)
from lsst.obs.base.tests import make_ramp_exposure_untrimmed
from lsst.obs.base.utils import InitialSkyWcsError, createInitialSkyWcs


class SimpleTestingTranslator(FitsTranslator, StubTranslator):
    """Simple `~astro_metadata_translator.MetadataTranslator` used for
    testing.
    """

    _const_map = {
        "boresight_rotation_angle": Angle(90 * u.deg),
        "boresight_rotation_coord": "sky",
        "detector_exposure_id": 12345,
        # The following are defined to prevent warnings about
        # undefined translators
        "dark_time": 0.0 * u.s,
        "exposure_time": 0.0 * u.s,
        "physical_filter": "u",
        "detector_num": 0,
        "detector_name": "0",
        "detector_group": "",
        "detector_unique_name": "0",
        "detector_serial": "",
        "observation_id": "--",
        "science_program": "unknown",
        "object": "unknown",
        "exposure_id": 0,
        "visit_id": 0,
        "relative_humidity": 30.0,
        "pressure": 0.0 * u.MPa,
        "temperature": 273 * u.K,
        "altaz_begin": None,
    }
    _trivial_map = {"boresight_airmass": "AIRMASS", "observation_type": "OBSTYPE"}

    def to_tracking_radec(self):
        radecsys = ("RADESYS",)
        radecpairs = (("RA", "DEC"),)
        return tracking_from_degree_headers(self, radecsys, radecpairs, unit=(u.deg, u.deg))


class MakeTestingRawVisitInfo(MakeRawVisitInfoViaObsInfo):
    """Test class for VisitInfo creation."""

    metadataTranslator = SimpleTestingTranslator


class SimpleFitsRawFormatter(FitsRawFormatterBase):
    """Simple test formatter for datastore interaction."""

    filterDefinitions = FilterDefinitionCollection(FilterDefinition(physical_filter="u", band="u"))

    @property
    def translatorClass(self):
        return SimpleTestingTranslator

    def getDetector(self, id):
        """Use CameraWrapper to create a fake detector that can map from
        PIXELS to FIELD_ANGLE.

        Always return Detector #10, so all the tests are self-consistent, and
        make sure it is in "assembled" form, since that's what the base
        formatter implementations assume.
        """
        return makeUpdatedDetector(CameraWrapper().camera.get(10))


class FitsRawFormatterTestCase(lsst.utils.tests.TestCase):
    """Test that we can read and write FITS files with butler."""

    def setUp(self):
        # The FITS WCS and VisitInfo coordinates in this header are
        # intentionally different, to make comparisons between them more
        # obvious.
        self.boresight = lsst.geom.SpherePoint(10.0, 20.0, lsst.geom.degrees)
        self.header = {
            "TELESCOP": "TEST",
            "INSTRUME": "UNKNOWN",
            "AIRMASS": 1.2,
            "RADESYS": "ICRS",
            "OBSTYPE": "science",
            "EQUINOX": 2000,
            "OBSGEO-X": "-5464588.84421314",
            "OBSGEO-Y": "-2493000.19137644",
            "OBSGEO-Z": "2150653.35350771",
            "RA": self.boresight.getLatitude().asDegrees(),
            "DEC": self.boresight.getLongitude().asDegrees(),
            "CTYPE1": "RA---SIN",
            "CTYPE2": "DEC--SIN",
            "CRPIX1": 5,
            "CRPIX2": 6,
            "CRVAL1": self.boresight.getLatitude().asDegrees() + 1,
            "CRVAL2": self.boresight.getLongitude().asDegrees() + 1,
            "CD1_1": 1e-5,
            "CD1_2": 0,
            "CD2_2": 1e-5,
            "CD2_1": 0,
        }
        # make a property list of the above, for use by the formatter.
        self.metadata = lsst.daf.base.PropertyList()
        self.metadata.update(self.header)

        maker = MakeTestingRawVisitInfo()
        self.visitInfo = maker(self.header)

        self.metadataSkyWcs = lsst.afw.geom.makeSkyWcs(self.metadata, strip=False)
        self.boresightSkyWcs = createInitialSkyWcs(self.visitInfo, CameraWrapper().camera.get(10))

        # set this to `contextlib.nullcontext()` to print the log warnings
        self.warnContext = self.assertLogs(level="WARNING")

        # Make a ref to pass to the formatter.
        universe = lsst.daf.butler.DimensionUniverse()
        dataId = lsst.daf.butler.DataCoordinate.standardize(
            instrument="Cam1", exposure=2, detector=10, physical_filter="u", band="u", universe=universe
        )
        datasetType = lsst.daf.butler.DatasetType(
            "dummy",
            dimensions=("instrument", "exposure", "detector"),
            storageClass=lsst.daf.butler.StorageClass(),
            universe=universe,
        )
        ref = lsst.daf.butler.DatasetRef(dataId=dataId, datasetType=datasetType, run="test")

        # We have no file in these tests, so make an empty descriptor.
        fileDescriptor = lsst.daf.butler.FileDescriptor(None, None)
        self.formatter = SimpleFitsRawFormatter(fileDescriptor, ref=ref)
        # Force the formatter's metadata to be what we've created above.
        self.formatter._metadata = self.metadata

    def test_makeWcs(self):
        detector = self.formatter.getDetector(1)
        wcs = self.formatter.makeWcs(self.visitInfo, detector)
        self.assertNotEqual(wcs, self.metadataSkyWcs)
        self.assertEqual(wcs, self.boresightSkyWcs)

    def test_makeWcs_if_metadata_is_bad(self):
        """Always use the VisitInfo WCS if available."""
        detector = self.formatter.getDetector(1)
        self.metadata.remove("CTYPE1")
        wcs = self.formatter.makeWcs(self.visitInfo, detector)
        self.assertNotEqual(wcs, self.metadataSkyWcs)
        self.assertEqual(wcs, self.boresightSkyWcs)

    def test_makeWcs_warn_if_visitInfo_is_None(self):
        """If VisitInfo is None, log a warning and use the metadata WCS."""
        detector = self.formatter.getDetector(1)
        with self.warnContext:
            wcs = self.formatter.makeWcs(None, detector)
        self.assertEqual(wcs, self.metadataSkyWcs)
        self.assertNotEqual(wcs, self.boresightSkyWcs)

    def test_makeWcs_fail_if_visitInfo_is_None(self):
        """If VisitInfo is None and metadata failed, raise an exception."""
        detector = self.formatter.getDetector(1)
        self.metadata.remove("CTYPE1")
        with self.warnContext, self.assertRaises(InitialSkyWcsError):
            self.formatter.makeWcs(None, detector)

    def test_makeWcs_fail_if_detector_is_bad(self):
        """If Detector is broken, raise an exception."""
        # This detector doesn't know about FIELD_ANGLE, so can't be used to
        # make a SkyWcs.
        detector = DetectorWrapper().detector
        with self.assertRaises(InitialSkyWcsError):
            self.formatter.makeWcs(self.visitInfo, detector)

    def test_amp_parameter(self):
        """Test loading subimages with the 'amp' parameter."""
        with lsst.utils.tests.getTempFilePath(".fits") as tmpFile:
            # Get a detector; this must be the same one that's baked into the
            # simple formatter at the top of this file, so that's how we get
            # it.
            detector = self.formatter.getDetector(1)
            # Make full exposure with ramp values and save just the image to
            # the temp file (with metadata), so it looks like a raw.
            full = make_ramp_exposure_untrimmed(detector)
            full.image.writeFits(tmpFile, metadata=self.metadata)
            # Loop over amps and try to read them via the formatter.
            for n, amp in enumerate(detector):
                for amp_parameter in [amp, amp.getName(), n]:
                    for parameters in [{"amp": amp_parameter}, {"amp": amp_parameter, "detector": detector}]:
                        with self.subTest(parameters=repr(parameters)):
                            # Make a new formatter that points at the new file
                            # and has the right parameters.
                            formatter = SimpleFitsRawFormatter(
                                lsst.daf.butler.FileDescriptor(
                                    lsst.daf.butler.Location(None, path=lsst.resources.ResourcePath(tmpFile)),
                                    lsst.daf.butler.StorageClassFactory().getStorageClass("ExposureI"),
                                    parameters=parameters,
                                ),
                                ref=self.formatter.dataset_ref,
                            )
                            subexp = formatter.read()
                            self.assertImagesEqual(subexp.image, full[amp.getRawBBox()].image)
                            self.assertEqual(len(subexp.getDetector()), 1)
                            self.assertAmplifiersEqual(subexp.getDetector()[0], amp)
                            self.assertEqual(subexp.visitInfo.id, 2)
                # We could try transformed amplifiers here that involve flips
                # and offsets, but:
                # - we already test the low-level code that does that in afw;
                # - we test very similar high-level code (which calls that
                #   same afw code) in the non-raw Exposure formatter, in
                #   test_butlerFits.py;
                # - the only instruments that actually have those kinds of
                #   amplifiers are those in obs_lsst, and that has a different
                #   raw formatter implementation that we need to test there
                #   anyway;
                # - these are kind of expensive tests.


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    """Check for file leaks."""


def setup_module(module):
    """Initialize pytest."""
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
