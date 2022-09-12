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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import unittest

import astropy.units as u
import lsst.afw.image
import numpy as np
from astro_metadata_translator import FitsTranslator, ObservationInfo, StubTranslator
from astropy.time import Time
from lsst.daf.base import DateTime
from lsst.obs.base import MakeRawVisitInfoViaObsInfo


class NewTranslator(FitsTranslator, StubTranslator):
    _trivial_map = {
        "exposure_time": "EXPTIME",
        "exposure_id": "EXP-ID",
        "observation_type": "OBSTYPE",
        "science_program": "PROGRAM",
        "observation_reason": "REASON",
        "object": "OBJECT",
        "has_simulated_content": "SIMULATE",
    }

    def to_location(self):
        return None

    def to_detector_exposure_id(self):
        return self.to_exposure_id()

    def to_focus_z(self):
        return self._header["FOCUSZ"]


class MakeTestableVisitInfo(MakeRawVisitInfoViaObsInfo):
    metadataTranslator = NewTranslator


class TestMakeRawVisitInfoViaObsInfo(unittest.TestCase):
    def setUp(self):
        # Reference values
        self.exposure_time = 6.2 * u.s
        self.exposure_id = 54321
        self.datetime_begin = Time("2001-01-02T03:04:05.123456789", format="isot", scale="utc")
        self.datetime_begin.precision = 9
        self.datetime_end = Time("2001-01-02T03:04:07.123456789", format="isot", scale="utc")
        self.datetime_end.precision = 9
        self.focus_z = 1.5 * u.mm

        self.header = {
            "DATE-OBS": self.datetime_begin.isot,
            "DATE-END": self.datetime_end.isot,
            "INSTRUME": "SomeCamera",
            "TELESCOP": "LSST",
            "TIMESYS": "UTC",
            "EXPTIME": self.exposure_time,
            "EXP-ID": self.exposure_id,
            "FOCUSZ": self.focus_z,
            "EXTRA1": "an abitrary key and value",
            "EXTRA2": 5,
            "OBSTYPE": "test type",
            "PROGRAM": "test program",
            "REASON": "test reason",
            "OBJECT": "test object",
            "SIMULATE": True,
        }

    def testMakeRawVisitInfoViaObsInfo(self):
        maker = MakeTestableVisitInfo()
        beforeLength = len(self.header)

        # Capture the warnings from StubTranslator since they are
        # confusing to people but irrelevant for the test.
        with self.assertWarns(UserWarning):
            with self.assertLogs(level="WARNING"):
                visitInfo = maker(self.header)

        self.assertAlmostEqual(visitInfo.getExposureTime(), self.exposure_time.to_value("s"))
        with self.assertWarns(FutureWarning):
            # TODO: tested for backward-compatibility; remove on DM-32138
            self.assertEqual(visitInfo.getExposureId(), self.exposure_id)
        self.assertEqual(visitInfo.id, self.exposure_id)
        self.assertEqual(visitInfo.getDate(), DateTime("2001-01-02T03:04:06.123456789Z", DateTime.UTC))
        # The header can possibly grow with header fix up provenance.
        self.assertGreaterEqual(len(self.header), beforeLength)
        self.assertEqual(visitInfo.getInstrumentLabel(), "SomeCamera")
        # Check focusZ with default value from astro_metadata_translator
        self.assertEqual(visitInfo.getFocusZ(), self.focus_z.to_value("mm"))
        self.assertEqual(visitInfo.observationType, "test type")
        self.assertEqual(visitInfo.scienceProgram, "test program")
        self.assertEqual(visitInfo.observationReason, "test reason")
        self.assertEqual(visitInfo.object, "test object")
        self.assertEqual(visitInfo.hasSimulatedContent, True)

    def testMakeRawVisitInfoViaObsInfo_empty(self):
        """Test that empty metadata fields are set to appropriate defaults."""
        maker = MakeTestableVisitInfo()
        # Capture the warnings from StubTranslator since they are
        # confusing to people but irrelevant for the test.
        with self.assertWarns(UserWarning):
            with self.assertLogs(level="WARNING"):
                visitInfo = maker({})
        self.assertTrue(np.isnan(visitInfo.focusZ))
        self.assertEqual(visitInfo.observationType, "")
        self.assertEqual(visitInfo.scienceProgram, "")
        self.assertEqual(visitInfo.observationReason, "")
        self.assertEqual(visitInfo.object, "")
        self.assertEqual(visitInfo.hasSimulatedContent, False)

    def testObservationInfo2VisitInfo(self):

        with self.assertWarns(UserWarning):
            obsInfo = ObservationInfo(self.header, translator_class=NewTranslator)

        # No log specified so no log message should appear
        visitInfo = MakeRawVisitInfoViaObsInfo.observationInfo2visitInfo(obsInfo)
        self.assertIsInstance(visitInfo, lsst.afw.image.VisitInfo)
        self.assertAlmostEqual(visitInfo.getExposureTime(), self.exposure_time.to_value("s"))
        with self.assertWarns(FutureWarning):
            # TODO: tested for backward-compatibility; remove on DM-32138
            self.assertEqual(visitInfo.getExposureId(), self.exposure_id)
        self.assertEqual(visitInfo.id, self.exposure_id)
        self.assertEqual(visitInfo.getDate(), DateTime("2001-01-02T03:04:06.123456789Z", DateTime.UTC))
        self.assertEqual(visitInfo.getInstrumentLabel(), "SomeCamera")
        # Check focusZ with default value from astro_metadata_translator
        self.assertEqual(visitInfo.getFocusZ(), self.focus_z.to_value("mm"))


if __name__ == "__main__":
    unittest.main()
