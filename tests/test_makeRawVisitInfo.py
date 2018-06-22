#
# LSST Data Management System
# Copyright 2016-2017 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import math
import unittest

import astropy.coordinates
import astropy.units

import lsst.utils.tests
import lsst.pex.exceptions
from lsst.daf.base import DateTime, PropertySet
from lsst.obs.base import MakeRawVisitInfo
from lsst.afw.geom import degrees


class SimpleMakeRawVisitInfo(MakeRawVisitInfo):

    def getDateAvg(self, md, exposureTime):
        """Return date at the middle of the exposure"""
        dateObs = self.popIsoDate(md, "DATE-OBS")
        return self.offsetDate(dateObs, 0.5*exposureTime)


def getMetadata(keyValDict):
    md = PropertySet()
    for key, val in keyValDict.items():
        md.set(key, val)
    return md


class VisitInfoTestCase(lsst.utils.tests.TestCase):
    """Test lsst.afw.image.VisitInfo, a simple struct-like class"""

    def setUp(self):
        self.makeRawVisitInfo = SimpleMakeRawVisitInfo()

    def testMakeRawVisitInfo(self):
        """Test base class functor MakeRawVisitInfo

        The base class only sets date and exposureTime fields,
        reads DATE-OBS, TIMESYS and EXPTIME,
        and deletes DATE-OBS and EXPTIME, but not TIMESYS.
        """
        exposureTime = 6.2  # arbitrary value in seconds
        date = DateTime("2001-01-02T03:04:05.123456789", DateTime.TAI)
        dateNsec = date.nsecs(DateTime.TAI)
        startDate = DateTime(dateNsec - int(exposureTime*1e9/2), DateTime.TAI)

        exposureId = 54321  # an arbitrary value
        md = getMetadata({
            "DATE-OBS": startDate.toString(DateTime.UTC),
            "TIMESYS": "UTC",
            "EXPTIME": exposureTime,
            "EXTRA1": "an abitrary key and value",
            "EXTRA2": 5,
        })
        visitInfo = self.makeRawVisitInfo(md=md, exposureId=exposureId)
        self.assertEqual(visitInfo.getExposureId(), exposureId)
        self.assertEqual(md.nameCount(), 3)  # TIMESYS and two EXTRAn keywords
        self.assertEqual(visitInfo.getExposureTime(), exposureTime)
        self.assertEqual(visitInfo.getDate(), date)

        # try TIMESYS=TAI
        md = getMetadata({
            "DATE-OBS": startDate.toString(DateTime.TAI),
            "TIMESYS": "TAI",
            "EXPTIME": exposureTime,
        })
        visitInfo = self.makeRawVisitInfo(md=md, exposureId=exposureId)
        self.assertEqual(md.nameCount(), 1)  # TIMESYS
        self.assertEqual(visitInfo.getExposureTime(), exposureTime)
        self.assertEqual(visitInfo.getDate(), date)

        # try omitting TIMESYS, which defaults to UTC
        md = getMetadata({
            "DATE-OBS": startDate.toString(DateTime.UTC),
            "EXPTIME": exposureTime,
        })
        visitInfo = self.makeRawVisitInfo(md=md, exposureId=exposureId)
        self.assertEqual(md.nameCount(), 0)
        self.assertEqual(visitInfo.getExposureTime(), exposureTime)
        self.assertEqual(visitInfo.getDate(), date)

        # omit DATE-OBS; date should be default-constructed
        md = getMetadata({
            "EXPTIME": exposureTime,
        })
        visitInfo = self.makeRawVisitInfo(md=md, exposureId=exposureId)
        self.assertEqual(md.nameCount(), 0)
        self.assertEqual(visitInfo.getExposureTime(), exposureTime)
        self.assertEqual(visitInfo.getDate(), DateTime())

        # omit EXPTIME; date should be start date, not avg date, and exposureTime should be nan
        md = getMetadata({
            "DATE-OBS": startDate.toString(DateTime.UTC),
        })
        visitInfo = self.makeRawVisitInfo(md=md, exposureId=exposureId)
        self.assertEqual(md.nameCount(), 0)
        self.assertTrue(math.isnan(visitInfo.getExposureTime()))
        self.assertEqual(visitInfo.getDate(), startDate)

    def testPopItem(self):
        md = getMetadata({
            "TIMESYS": "UTC",
            "OTHER": 5,
        })

        timesys = self.makeRawVisitInfo.popItem(md, "TIMESYS")
        self.assertEqual(timesys, "UTC")
        self.assertEqual(set(md.names()), set(["OTHER"]))

        defVal = self.makeRawVisitInfo.popItem(md, "BADKEY", default=7)
        self.assertEqual(set(md.names()), set(["OTHER"]))
        self.assertEqual(defVal, 7)
        missingItem = self.makeRawVisitInfo.popItem(md, "BADKEY")
        self.assertIsNone(missingItem)

    def testPopFloat(self):
        dataDict = {
            "FLOAT": 5.5,
            "INT": 5,
            "FLOATSTR": "6.1",
            "INTSTR": "6",
            "STR": "FOO",
        }

        for key, desValue in dataDict.items():
            if key == "STR":
                continue
            md = getMetadata(dataDict)
            value = self.makeRawVisitInfo.popFloat(md, key)
            self.assertAlmostEqual(value, float(desValue))
            self.assertEqual(len(md.names()), 4)

        badFloat = self.makeRawVisitInfo.popFloat(md, "STR")
        self.assertTrue(math.isnan(badFloat))

        missingValue = self.makeRawVisitInfo.popFloat(md, "BADKEY")
        self.assertTrue(math.isnan(missingValue))

    def testPopAngle(self):
        dataDict = {
            "DEG1": 270.5,
            "DEG2": "45:30",
            "DEG3": "-310:12:32",
            "HR_1": -23.5,
            "HR_2": "23:30",
            "HR_3": "-13:15:16.7",
            "STR": "FOO",
        }

        for key, desValue in dataDict.items():
            if key == "STR":
                continue
            elif key.startswith("DEG"):
                units = astropy.units.deg
            else:
                units = astropy.units.h

            desAngleDeg = astropy.coordinates.Angle(desValue, unit=units).deg
            md = getMetadata(dataDict)
            angle = self.makeRawVisitInfo.popAngle(md, key, units=units)
            self.assertAnglesAlmostEqual(angle, desAngleDeg*degrees)

        badAngle = self.makeRawVisitInfo.popAngle(md, "STR")
        self.assertTrue(math.isnan(badAngle.asDegrees()))

        missingAngle = self.makeRawVisitInfo.popAngle(md, "BADKEY")
        self.assertTrue(math.isnan(missingAngle.asDegrees()))

    def testPopIsoDate(self):
        for timesys in (None, "UTC", "TAI"):
            dataDict = {
                "DATE1": "2001-02-03T04:05:06.123456789",
                "DATE2": "2001-02-03T04:05:06.123456789Z",  # Z should be ignored
                "DATE3": "1980-03-04T01:02:03.999999999",
                "BADISODATE": "51234.354",
            }
            if timesys is not None:
                dataDict["TIMESYS"] = timesys

            for key, dateStr in dataDict.items():
                if not key.startswith("DATE"):
                    continue

                lsstSys = dict(
                    UTC=DateTime.UTC,
                    TAI=DateTime.TAI,
                ).get(timesys, DateTime.UTC)

                # lsstDateStr = dateStr with trailing Z if UTC, else no trailing Z,
                # because lsst.daf.base.DateTime is very picky
                lsstDateStr = dateStr
                if lsstSys == DateTime.UTC:
                    if not lsstDateStr.endswith("Z"):
                        lsstDateStr = lsstDateStr + "Z"
                elif dateStr.endswith("Z"):
                    lsstDateStr = lsstDateStr[0:-1]
                desDate = DateTime(lsstDateStr, lsstSys)

                md = getMetadata(dataDict)
                numKeys = len(md.names())
                date = self.makeRawVisitInfo.popIsoDate(md, key, timesys=timesys)
                self.assertEqual(len(md.names()), numKeys - 1)
                self.assertEqual(date, desDate)

            badDate = self.makeRawVisitInfo.popIsoDate(md, "BADISODATE")
            self.assertEqual(badDate, DateTime())

            missingDate = self.makeRawVisitInfo.popIsoDate(md, "BADKEY")
            self.assertEqual(missingDate, DateTime())

    def testPopMjdDate(self):
        for timesys in (None, "UTC", "TAI"):

            dataDict = {
                "DATE1": 51943.1705801,
                "DATE2": 44302.0433218,
                "BADMJDDATE": "2001-02-03T04:05:06.123456789",
            }
            if timesys is not None:
                dataDict["TIMESYS"] = timesys

            for key, mjdDate in dataDict.items():
                if not key.startswith("DATE"):
                    continue

                lsstSys = dict(
                    UTC=DateTime.UTC,
                    TAI=DateTime.TAI,
                ).get(timesys, DateTime.UTC)

                desDate = DateTime(mjdDate, DateTime.MJD, lsstSys)

                md = getMetadata(dataDict)
                numKeys = len(md.names())
                date = self.makeRawVisitInfo.popMjdDate(md, key, timesys=timesys)
                self.assertEqual(len(md.names()), numKeys - 1)
                self.assertAlmostEqual(date.get(), desDate.get())

            badDate = self.makeRawVisitInfo.popMjdDate(md, "BADMJDDATE")
            self.assertEqual(badDate, DateTime())

            missingDate = self.makeRawVisitInfo.popMjdDate(md, "BADKEY")
            self.assertEqual(missingDate, DateTime())

    def testEraFromLstAndLongitude(self):
        LST = 90*degrees
        Longitude = 50*degrees
        era = self.makeRawVisitInfo.eraFromLstAndLongitude(LST, Longitude)
        self.assertAnglesAlmostEqual(era, LST-Longitude)

    def testEraFromLstAndLongitude_float_vs_Angle_fails(self):
        val1 = 90*degrees
        val2 = 50.0
        with self.assertRaises(TypeError):
            self.makeRawVisitInfo.eraFromLstAndLongitude(val1, val2)
        with self.assertRaises(TypeError):
            self.makeRawVisitInfo.eraFromLstAndLongitude(val2, val1)

    def testAltitudeFromZenithDistance(self):
        for zdDeg in (0, 35.6, 89.999, 90.0):
            desAltDeg = 90-zdDeg
            self.assertAnglesAlmostEqual(
                desAltDeg*degrees,
                self.makeRawVisitInfo.altitudeFromZenithDistance(zdDeg*degrees),
            )

    def testCentigradeFromKelvin(self):
        for tempK, desTempC in (  # a few values from http://www.convertunits.com/from/kelvin/to/centigrade
            (0, -273.15),
            (301.5, 28.35),
        ):
            self.assertAlmostEqual(desTempC, self.makeRawVisitInfo.centigradeFromKelvin(tempK))

    def testPascalFromMmHg(self):
        for mmHg, desPascal in (  # a few values from http://www.convertunits.com/from/mm+Hg/to/pascal
            (1, 133.32239),
            (0.062, 8.26598818),
        ):
            self.assertAlmostEqual(desPascal, self.makeRawVisitInfo.pascalFromMmHg(mmHg), places=5)


def setup_module(module):
    lsst.utils.tests.init()


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
