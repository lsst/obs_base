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

import astropy.coordinates
import astropy.units as u
import numpy as np
from astro_metadata_translator import FitsTranslator, ObservationInfo, StubTranslator
from astro_metadata_translator.translators.helpers import (
    altaz_from_degree_headers,
    tracking_from_degree_headers,
)
from astropy.coordinates import EarthLocation
from astropy.time import Time

import lsst.afw.image
from lsst.daf.base import DateTime
from lsst.obs.base import MakeRawVisitInfoViaObsInfo


def _q_to_float(q: u.Quantity, unit: u.UnitBase | str) -> float:
    values = q.to_value(unit=unit)
    if isinstance(values, np.ndarray):
        raise ValueError(
            f"Converting quantity to a float failed because unexpectedly got more than one float: {values}"
        )
    return float(values)


class NewTranslator(FitsTranslator, StubTranslator):
    """Metadata translator to use for tests."""

    _trivial_map = {
        "exposure_time": ("EXPTIME", {"unit": u.s}),
        "dark_time": ("DARKTM", {"unit": u.s}),
        "exposure_id": "EXP-ID",
        "boresight_airmass": "AMASS",
        "boresight_rotation_angle": ("SKYROT", {"unit": u.deg}),
        "boresight_rotation_coord": "SKYC",
        "observation_type": "OBSTYPE",
        "science_program": "PROGRAM",
        "observation_reason": "REASON",
        "object": "OBJECT",
        "has_simulated_content": "SIMULATE",
        "relative_humidity": "RELHUM",
        "temperature": ("AIR_TEMP", {"unit": u.deg_C}),
        "pressure": ("PRESSURE", {"unit": u.kPa}),
        "focus_z": ("FOCUSZ", {"unit": u.mm}),
    }

    def to_location(self):
        return EarthLocation.from_geodetic(-70.747698, -30.244728, 2663.0)

    def to_detector_exposure_id(self):
        return self.to_exposure_id()

    def to_tracking_radec(self) -> astropy.coordinates.SkyCoord | None:
        return tracking_from_degree_headers(self, ["RADESYS"], (("RA_DEG", "DEC_DEG"),))

    def to_altaz_begin(self) -> astropy.coordinates.AltAz | None:
        return altaz_from_degree_headers(self, (("TELALT", "TELAZ"),), self.to_datetime_begin())


class MakeTestableVisitInfo(MakeRawVisitInfoViaObsInfo):
    """Test class for VisitInfo construction."""

    metadataTranslator = NewTranslator


class TestMakeRawVisitInfoViaObsInfo(unittest.TestCase):
    """Test VisitInfo construction."""

    def setUp(self):
        # Reference values
        self.exposure_time = 6.2 * u.s
        self.dark_time = 7.0 * u.s
        self.boresight_airmass = 1.5
        self.exposure_id = 54321
        self.datetime_begin = Time("2001-01-02T03:04:05.123456789", format="isot", scale="utc")
        self.datetime_begin.precision = 9
        self.datetime_end = self.datetime_begin + self.exposure_time
        self.datetime_end.precision = 9
        self.focus_z = 1.5 * u.mm
        self.pressure = 101.0 * u.kPa

        self.header = {
            "DATE-OBS": self.datetime_begin.isot,
            "DATE-END": self.datetime_end.isot,
            "INSTRUME": "SomeCamera",
            "TELESCOP": "LSST",
            "TIMESYS": "UTC",
            "SKYROT": 45.0,
            "SKYC": "sky",
            "EXPTIME": _q_to_float(self.exposure_time, u.s),
            "DARKTM": _q_to_float(self.dark_time, u.s),
            "EXP-ID": self.exposure_id,
            "FOCUSZ": _q_to_float(self.focus_z, u.mm),
            "EXTRA1": "an abitrary key and value",
            "EXTRA2": 5,
            "OBSTYPE": "test type",
            "PROGRAM": "test program",
            "REASON": "test reason",
            "OBJECT": "test object",
            "SIMULATE": True,
            "RELHUM": 42.0,
            "AIR_TEMP": 1.0,
            "PRESSURE": _q_to_float(self.pressure, u.kPa),
            "AMASS": self.boresight_airmass,
            "RADESYS": "FK5",
            "RA_DEG": 45.0,
            "DEC_DEG": -30.0,
            "TELALT": 60.0,
            "TELAZ": 180.0,
        }

    def testMakeRawVisitInfoViaObsInfo(self):
        maker = MakeTestableVisitInfo()
        beforeLength = len(self.header)

        # Capture the warnings from StubTranslator since they are
        # confusing to people but irrelevant for the test.
        with self.assertWarns(UserWarning):
            visitInfo = maker(self.header)

        self.assertAlmostEqual(visitInfo.getExposureTime(), _q_to_float(self.exposure_time, u.s))
        self.assertAlmostEqual(visitInfo.getDarkTime(), _q_to_float(self.dark_time, u.s))
        self.assertEqual(visitInfo.id, self.exposure_id)
        self.assertEqual(visitInfo.getDate(), DateTime("2001-01-02T03:04:08.223456789Z", DateTime.UTC))
        # The header can possibly grow with header fix up provenance.
        self.assertGreaterEqual(len(self.header), beforeLength)
        self.assertEqual(visitInfo.getInstrumentLabel(), "SomeCamera")
        self.assertAlmostEqual(visitInfo.getBoresightRaDec().getRa().asDegrees(), 45.0, places=5)
        self.assertAlmostEqual(visitInfo.getBoresightRaDec().getDec().asDegrees(), -30.0, places=5)
        self.assertAlmostEqual(visitInfo.getBoresightAzAlt().getLongitude().asDegrees(), 180.0, places=7)
        self.assertAlmostEqual(visitInfo.getBoresightAzAlt().getLatitude().asDegrees(), 60.0, places=7)
        self.assertEqual(visitInfo.getBoresightAirmass(), self.boresight_airmass)
        self.assertAlmostEqual(visitInfo.getBoresightRotAngle().asDegrees(), 45.0)
        self.assertEqual(visitInfo.getRotType(), lsst.afw.image.RotType.SKY)
        self.assertAlmostEqual(visitInfo.getObservatory().getLongitude().asDegrees(), -70.747698)
        self.assertAlmostEqual(visitInfo.getObservatory().getLatitude().asDegrees(), -30.244728)
        self.assertAlmostEqual(visitInfo.getObservatory().getElevation(), 2663.0)
        self.assertEqual(visitInfo.getWeather().getHumidity(), 42.0)
        self.assertEqual(visitInfo.getWeather().getAirTemperature(), 1.0)
        self.assertEqual(visitInfo.getWeather().getAirPressure(), _q_to_float(self.pressure, u.Pa))
        # Check focusZ with default value from astro_metadata_translator
        self.assertEqual(visitInfo.getFocusZ(), _q_to_float(self.focus_z, u.mm))
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

        # Check that a couple of values look correct.
        visitInfo = MakeRawVisitInfoViaObsInfo.observationInfo2visitInfo(obsInfo)
        self.assertIsInstance(visitInfo, lsst.afw.image.VisitInfo)
        self.assertAlmostEqual(visitInfo.getExposureTime(), _q_to_float(self.exposure_time, u.s))
        self.assertEqual(visitInfo.id, self.exposure_id)

        # Convert it back to an ObservationInfo.
        obsInfo2 = MakeRawVisitInfoViaObsInfo.visitInfo2observationInfo(visitInfo)

        # We can not compare all fields because some fields are lost in
        # conversion to VisitInfo.
        to_compare = (
            "datetime_begin",
            "altaz_begin",
            "boresight_airmass",
            "boresight_rotation_angle",
            "boresight_rotation_coord",
            "dark_time",
            "datetime_end",
            "exposure_group",
            "exposure_id",
            "exposure_time",
            "exposure_time_requested",
            "focus_z",
            "has_simulated_content",
            "location",
            "object",
            "observation_reason",
            "observation_type",
            "observing_day",
            "pressure",
            "relative_humidity",
            "science_program",
            "temperature",
            "tracking_radec",
        )

        for prop in to_compare:
            previous = getattr(obsInfo, prop)
            new = getattr(obsInfo2, prop)
            with self.subTest(prop=prop):
                match new:
                    case astropy.coordinates.AltAz():
                        self.assertAlmostEqual(new.az, previous.az)
                        self.assertAlmostEqual(new.alt, previous.alt)
                    case astropy.coordinates.SkyCoord():
                        self.assertAlmostEqual(new.icrs.ra, previous.icrs.ra)
                        self.assertAlmostEqual(new.icrs.dec, previous.icrs.dec)
                    case astropy.coordinates.EarthLocation():
                        for new_pos, previous_pos in zip(
                            new.to_geocentric(), previous.to_geocentric(), strict=True
                        ):
                            self.assertAlmostEqual(new_pos, previous_pos)
                    case astropy.time.Time():
                        self.assertEqual(new, previous)
                    case float() | u.Quantity():
                        self.assertEqual(new, previous, f"Comparing float-like property {prop}")
                    case int() | str():
                        self.assertEqual(new, previous, f"Comparing int or string property {prop}")
                    case _:
                        raise RuntimeError(f"Encountered unexpected type for property {prop}")


if __name__ == "__main__":
    unittest.main()
