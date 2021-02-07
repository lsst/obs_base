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

from astropy.coordinates import Angle
import astropy.units as u

import lsst.utils.tests

from astro_metadata_translator import FitsTranslator, StubTranslator
from astro_metadata_translator.translators.helpers import tracking_from_degree_headers
from lsst.afw.cameraGeom.testUtils import CameraWrapper, DetectorWrapper
import lsst.afw.geom
import lsst.daf.base
import lsst.daf.butler
import lsst.geom
from lsst.obs.base import FitsRawFormatterBase, MakeRawVisitInfoViaObsInfo, FilterDefinitionCollection
from lsst.obs.base.utils import createInitialSkyWcs, InitialSkyWcsError


class SimpleTestingTranslator(FitsTranslator, StubTranslator):
    _const_map = {"boresight_rotation_angle": Angle(90*u.deg),
                  "boresight_rotation_coord": "sky",
                  "detector_exposure_id": 12345,
                  # The following are defined to prevent warnings about
                  # undefined translators
                  "dark_time": 0.0*u.s,
                  "exposure_time": 0.0*u.s,
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
                  "pressure": 0.0*u.MPa,
                  "temperature": 273*u.K,
                  "altaz_begin": None,
                  }
    _trivial_map = {"boresight_airmass": "AIRMASS",
                    "observation_type": "OBSTYPE"}

    def to_tracking_radec(self):
        radecsys = ("RADESYS", )
        radecpairs = (("RA", "DEC"),)
        return tracking_from_degree_headers(self, radecsys, radecpairs, unit=(u.deg, u.deg))


class MakeTestingRawVisitInfo(MakeRawVisitInfoViaObsInfo):
    metadataTranslator = SimpleTestingTranslator


class SimpleFitsRawFormatter(FitsRawFormatterBase):
    filterDefinitions = FilterDefinitionCollection()

    @property
    def translatorClass(self):
        return SimpleTestingTranslator

    def getDetector(self, id):
        """Use CameraWrapper to create a fake detector that can map from
        PIXELS to FIELD_ANGLE.

        Always return Detector #10, so all the tests are self-consistent.
        """
        return CameraWrapper().camera.get(10)


class FitsRawFormatterTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        # reset the filters before we test anything
        FilterDefinitionCollection.reset()

        # The FITS WCS and VisitInfo coordinates in this header are
        # intentionally different, to make comparisons between them more
        # obvious.
        self.boresight = lsst.geom.SpherePoint(10., 20., lsst.geom.degrees)
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
            "CD2_1": 0
        }
        # make a property list of the above, for use by the formatter.
        self.metadata = lsst.daf.base.PropertyList()
        self.metadata.update(self.header)

        maker = MakeTestingRawVisitInfo()
        self.visitInfo = maker(self.header)

        self.metadataSkyWcs = lsst.afw.geom.makeSkyWcs(self.metadata, strip=False)
        self.boresightSkyWcs = createInitialSkyWcs(self.visitInfo, CameraWrapper().camera.get(10))

        # set these to `contextlib.nullcontext()` to print the log warnings
        self.warnContext = self.assertLogs(level="WARNING")
        self.logContext = lsst.log.UsePythonLogging()

        # Make a data ID to pass to the formatter.
        universe = lsst.daf.butler.DimensionUniverse()
        dataId = lsst.daf.butler.DataCoordinate.standardize(instrument="Cam1", exposure=2, detector=10,
                                                            physical_filter="u", band="u", universe=universe)

        # We have no file in these tests, so make an empty descriptor.
        fileDescriptor = lsst.daf.butler.FileDescriptor(None, None)
        self.formatter = SimpleFitsRawFormatter(fileDescriptor, dataId)
        # Force the formatter's metadata to be what we've created above.
        self.formatter._metadata = self.metadata

    def test_makeWcs(self):
        detector = self.formatter.getDetector(1)
        wcs = self.formatter.makeWcs(self.visitInfo, detector)
        self.assertNotEqual(wcs, self.metadataSkyWcs)
        self.assertEqual(wcs, self.boresightSkyWcs)

    def test_makeWcs_warn_if_metadata_is_bad(self):
        """If the metadata is bad, log a warning and use the VisitInfo WCS.
        """
        detector = self.formatter.getDetector(1)
        self.metadata.remove("CTYPE1")
        with self.warnContext, self.logContext:
            wcs = self.formatter.makeWcs(self.visitInfo, detector)
        self.assertNotEqual(wcs, self.metadataSkyWcs)
        self.assertEqual(wcs, self.boresightSkyWcs)

    def test_makeWcs_warn_if_visitInfo_is_None(self):
        """If VisitInfo is None, log a warning and use the metadata WCS.
        """
        detector = self.formatter.getDetector(1)
        with self.warnContext, self.logContext:
            wcs = self.formatter.makeWcs(None, detector)
        self.assertEqual(wcs, self.metadataSkyWcs)
        self.assertNotEqual(wcs, self.boresightSkyWcs)

    def test_makeWcs_fail_if_visitInfo_is_None(self):
        """If VisitInfo is None and metadata failed, raise an exception.
        """
        detector = self.formatter.getDetector(1)
        self.metadata.remove("CTYPE1")
        with self.warnContext, self.logContext, self.assertRaises(InitialSkyWcsError):
            self.formatter.makeWcs(None, detector)

    def test_makeWcs_fail_if_detector_is_bad(self):
        """If Detector is broken, raise an exception.
        """
        # This detector doesn't know about FIELD_ANGLE, so can't be used to
        # make a SkyWcs.
        detector = DetectorWrapper().detector
        with self.assertRaises(InitialSkyWcsError):
            self.formatter.makeWcs(self.visitInfo, detector)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == '__main__':
    lsst.utils.tests.init()
    unittest.main()
