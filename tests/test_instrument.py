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

"""Tests of the Instrument class."""

import datetime
import unittest

import astropy.units as u
from astro_metadata_translator import ObservationInfo
from astropy.coordinates import SkyCoord
from astropy.time import Time

from lsst.daf.butler import DimensionUniverse
from lsst.obs.base import Instrument, makeExposureRecordFromObsInfo
from lsst.obs.base.instrument_tests import DummyCam, InstrumentTestData, InstrumentTests


class InstrumentTestCase(InstrumentTests, unittest.TestCase):
    """Test for Instrument."""

    instrument = DummyCam()

    data = InstrumentTestData(
        name="DummyCam", nDetectors=2, firstDetectorName="RXX_S00", physical_filters={"dummy_g", "dummy_u"}
    )

    def test_getCamera(self):
        """No camera defined in DummyCam."""
        return

    def test_collectionTimestamps(self):
        self.assertEqual(
            Instrument.formatCollectionTimestamp("2018-05-03"),
            "20180503T000000Z",
        )
        self.assertEqual(
            Instrument.formatCollectionTimestamp("2018-05-03T14:32:16"),
            "20180503T143216Z",
        )
        self.assertEqual(
            Instrument.formatCollectionTimestamp("20180503T143216Z"),
            "20180503T143216Z",
        )
        self.assertEqual(
            Instrument.formatCollectionTimestamp(datetime.datetime(2018, 5, 3, 14, 32, 16)),
            "20180503T143216Z",
        )
        formattedNow = Instrument.makeCollectionTimestamp()
        self.assertIsInstance(formattedNow, str)
        datetimeThen1 = datetime.datetime.strptime(formattedNow, "%Y%m%dT%H%M%S%z")
        self.assertEqual(datetimeThen1.tzinfo, datetime.UTC)

    def test_group_name(self):
        """Test group name to ID conversion."""
        self.assertEqual(self.instrument.group_name_to_group_id("1:234-5.6"), 123456)
        with self.assertRaises(ValueError):
            self.instrument.group_name_to_group_id("no_int")

    def test_sky_angle_case_insensitive(self):
        """Test that makeExposureRecordFromObsInfo transfers sky_angle
        regardless of the boresight_rotation_coord capitalization.
        """
        obs_info = ObservationInfo(
            instrument="DummyCam",
            exposure_id=1,
            exposure_time=10.0 * u.s,
            exposure_time_requested=10.0 * u.s,
            datetime_begin=Time("2025-01-01T00:00:00", scale="utc"),
            datetime_end=Time("2025-01-01T00:00:10", scale="utc"),
            observing_day=20250101,
            observation_type="science",
            physical_filter="g",
            tracking_radec=SkyCoord("00:00:00.0 +00:00:00.0", unit="deg"),
            boresight_rotation_angle=45.0 * u.deg,
            boresight_rotation_coord="SKY",
        )
        record = makeExposureRecordFromObsInfo(obs_info, DimensionUniverse())
        self.assertEqual(record.sky_angle, 45.0)


if __name__ == "__main__":
    unittest.main()
