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

"""Tests of the Instrument class.
"""

import datetime
import unittest

from lsst.obs.base import Instrument
from lsst.obs.base.instrument_tests import DummyCam, InstrumentTestData, InstrumentTests


class InstrumentTestCase(InstrumentTests, unittest.TestCase):
    """Test for Instrument."""

    instrument = DummyCam()

    data = InstrumentTestData(
        name="DummyCam", nDetectors=2, firstDetectorName="RXX_S00", physical_filters={"dummy_g", "dummy_u"}
    )

    def test_getCamera(self):
        """No camera defined in DummyCam"""
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
        self.assertEqual(datetimeThen1.tzinfo, datetime.timezone.utc)


if __name__ == "__main__":
    unittest.main()
