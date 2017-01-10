#
# LSST Data Management System
# Copyright 2016 LSST Corporation.
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
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import math

import astropy.coordinates
import astropy.time
import astropy.units

from lsst.log import Log
from lsst.daf.base import DateTime
from lsst.afw.geom import degrees
from lsst.afw.image import makeVisitInfo

__all__ = ["MakeRawVisitInfo"]


PascalPerMillibar = 100.0
PascalPerMmHg = 133.322387415  # from Wikipedia; exact
PascalPerTorr = 101325.0/760.0  # from Wikipedia; exact
KelvinMinusCentigrade = 273.15  # from Wikipedia; exact

# have these read at need, to avoid unexpected errors later
NaN = float("nan")
BadDate = DateTime()


class MakeRawVisitInfo(object):
    """Base class functor to make a VisitInfo from the FITS header of a raw image

    A subclass will be wanted for each camera. Subclasses should override
    - setArgDict: the override can call the base implementation,
                    which simply sets exposure time and date of observation
    - setDateAvg

    The design philosophy is to make a best effort and log warnings of problems,
    rather than raising exceptions, in order to extract as much VisitInfo information as possible
    from a messy FITS header without the user needing to add a lot of error handling.

    However, the methods that transform units are less forgiving; they assume
    the user provides proper data types, since type errors in arguments to those
    are almost certainly due to coding mistakes.
    """

    def __init__(self, log=None):
        """Construct a MakeRawVisitInfo
        """
        if log is None:
            log = Log.getLogger("MakeRawVisitInfo")
        self.log = log

    def __call__(self, md, exposureId):
        """Construct a VisitInfo and strip associated data from the metadata

        @param[in,out] md  metadata, as an lsst.daf.base.PropertyList or PropertySet;
            items that are used are stripped from the metadata
            (except TIMESYS, because it may apply to more than one other keyword).
        @param[in] exposureId  exposure ID

        The basic implementation sets date and exposureTime using typical values
        found in FITS files and logs a warning if neither can be set.
        """
        argDict = dict(exposureId=exposureId)
        self.setArgDict(md, argDict)
        for key in list(argDict.keys()):  # use a copy because we may delete items
            if argDict[key] is None:
                self.log.warn("argDict[{}] is None; stripping".format(key, argDict[key]))
                del argDict[key]
        return makeVisitInfo(**argDict)

    def setArgDict(self, md, argDict):
        """Fill an argument dict with arguments for makeVisitInfo and pop associated metadata

        Subclasses are expected to override this method, though the override
        may wish to call this default implementation, which:
        - sets exposureTime from "EXPTIME"
        - sets date by calling getDateAvg

        @param[in,out] md  metadata, as an lsst.daf.base.PropertyList or PropertySet;
            items that are used are stripped from the metadata
            (except TIMESYS, because it may apply to more than one other keyword).
        @param[in,out] argdict  a dict of arguments

        Subclasses should expand this or replace it.
        """
        argDict["exposureTime"] = self.popFloat(md, "EXPTIME")
        darkTime = self.popFloat(md, "DARKTIME") if md.exists("DARKTIME") else NaN
        argDict["darkTime"] = darkTime if not math.isnan(darkTime) else argDict["exposureTime"]
        argDict["date"] = self.getDateAvg(md=md, exposureTime=argDict["exposureTime"])

    def getDateAvg(self, md, exposureTime):
        """Return date at the middle of the exposure

        @param[in,out] md  metadata, as an lsst.daf.base.PropertyList or PropertySet;
            items that are used are stripped from the metadata
            (except TIMESYS, because it may apply to more than one other keyword).
        @param[in] exposureTime  exposure time (sec)

        Subclasses must override. Here is a typical implementation:
        dateObs = self.popIsoDate(md, "DATE-OBS")
        return self.offsetDate(dateObs, 0.5*exposureTime)
        """
        raise NotImplementedError()

    def offsetDate(self, date, offsetSec):
        """Return a date offset by a specified number of seconds

        @param[in] date  date (an lsst.daf.base.DateTime)
        @param[in] offsetSec  offset, in seconds (float)
        @return the offset date (an lsst.daf.base.DateTime)
        """
        if not date.isValid():
            self.log.warn("date is invalid; cannot offset it")
            return date
        if math.isnan(offsetSec):
            self.log.warn("offsetSec is invalid; cannot offset date")
            return date
        dateNSec = date.nsecs(DateTime.TAI)
        return DateTime(dateNSec + int(offsetSec*1.0e9), DateTime.TAI)

    def popItem(self, md, key, default=None):
        """Remove an item of metadata and return the value

        @param[in,out] md  metadata, as an lsst.daf.base.PropertyList or PropertySet;
            the popped key is removed
        @param[in] key  metadata key
        @param[in] default  default value to return if key not found; ignored if doRaise true
        @return the value of the specified key, using whatever type md.get(key) returns

        Log a warning if the key is not found
        """
        try:
            if not md.exists(key):
                self.log.warn("Key=\"{}\" not in metadata".format(key))
                return default
            val = md.get(key)
            md.remove(key)
            return val
        except Exception as e:
            # this should never happen, but is a last ditch attempt to avoid exceptions
            self.log.warn("Could not read key=\"{}\" in metadata: {}" % (key, e))
        return default

    def popFloat(self, md, key):
        """Pop a float with a default of Nan

        @param[in,out] md  metadata, as an lsst.daf.base.PropertyList or PropertySet
        @param[in] key  date key to read and remove from md
        @return the value of the specified key as a float; float("nan") if the key is not found
        """
        val = self.popItem(md, key, default=NaN)
        try:
            return float(val)
        except Exception as e:
            self.log.warn("Could not interpret {} value {} as a float: {}".format(key, repr(val), e))
        return NaN

    def popAngle(self, md, key, units=astropy.units.deg):
        """Pop an lsst.afw.geom.Angle, whose metadata is in the specified units, with a default of Nan

        The angle may be specified as a float or sexagesimal string with 1-3 fields.

        @param[in,out] md  metadata, as an lsst.daf.base.PropertyList or PropertySet
        @param[in] key  date key to read and remove from md
        @return angle, as an lsst.afw.geom.Angle; Angle(NaN) if the key is not found
        """
        angleStr = self.popItem(md, key, default=None)
        if angleStr is not None:
            try:
                return (astropy.coordinates.Angle(angleStr, unit=units).deg)*degrees
            except Exception as e:
                self.log.warn("Could not intepret {} value {} as an angle: {}".format(key, repr(angleStr), e))
        return NaN*degrees

    def popIsoDate(self, md, key, timesys=None):
        """Pop a FITS ISO date as an lsst.daf.base.DateTime

        @param[in,out] md  metadata, as an lsst.daf.base.PropertyList or PropertySet
        @param[in] key  date key to read and remove from md
        @param[in] timesys  time system as a string (not case sensitive), e.g. "UTC" or None;
                    if None then look for TIMESYS (but do NOT pop it, since it may be used
                    for more than one date) and if not found, use UTC
        @return date as an lsst.daf.base.DateTime; DateTime() if the key is not found
        """
        isoDateStr = self.popItem(md=md, key=key)
        if isoDateStr is not None:
            try:
                if timesys is None:
                    timesys = md.get("TIMESYS") if md.exists("TIMESYS") else "UTC"
                if isoDateStr.endswith("Z"):  # illegal in FITS
                    isoDateStr = isoDateStr[0:-1]
                astropyTime = astropy.time.Time(isoDateStr, scale=timesys.lower(), format="fits")
                # DateTime uses nanosecond resolution, regardless of the resolution of the original date
                astropyTime.precision = 9
                # isot is ISO8601 format with "T" separating date and time and no time zone
                return DateTime(astropyTime.tai.isot, DateTime.TAI)
            except Exception as e:
                self.log.warn("Could not parse {} = {} as an ISO date: {}".format(key, isoDateStr, e))
        return BadDate

    def popMjdDate(self, md, key, timesys=None):
        """Get a FITS MJD date as an lsst.daf.base.DateTime

        @param[in,out] md  metadata, as an lsst.daf.base.PropertyList or PropertySet
        @param[in] dateKey  date key to read and remove from md
        @param[in] timesys  time system as a string, e.g. "UTC" or None;
                    if None then look for TIMESYS (but do NOT pop it, since it may be used
                    for more than one date) and if not found, use UTC
        @return date as an lsst.daf.base.DateTime; DateTime() if the key is not found
        """
        mjdDate = self.popFloat(md, key)
        try:
            if timesys is None:
                timesys = md.get("TIMESYS") if md.exists("TIMESYS") else "UTC"
            astropyTime = astropy.time.Time(mjdDate, format="mjd", scale=timesys.lower())
            # DateTime uses nanosecond resolution, regardless of the resolution of the original date
            astropyTime.precision = 9
            # isot is ISO8601 format with "T" separating date and time and no time zone
            return DateTime(astropyTime.tai.isot, DateTime.TAI)
        except Exception as e:
            self.log.warn("Could not parse {} = {} as an MJD date: {}".format(key, mjdDate, e))
        return BadDate

    @staticmethod
    def eraFromLstAndLongitude(lst, longitude):
        """
        Return an approximate Earth Rotation Angle (afw:Angle) computed from
        local sidereal time and longitude (both as afw:Angle; Longitude shares
        the afw:Observatory covention: positive values are E of Greenwich).

        NOTE: if we properly compute ERA via UT1 a la DM-8053, we should remove
        this method.
        """
        return lst - longitude

    @staticmethod
    def altitudeFromZenithDistance(zd):
        """Convert zenith distance to altitude (lsst.afw.geom.Angle)"""
        return 90*degrees - zd

    @staticmethod
    def centigradeFromKelvin(tempK):
        """Convert temperature from Kelvin to Centigrade"""
        return tempK - KelvinMinusCentigrade

    @staticmethod
    def pascalFromMBar(mbar):
        """Convert pressure from millibars to Pascals
        """
        return mbar*PascalPerMillibar

    @staticmethod
    def pascalFromMmHg(mmHg):
        """Convert pressure from mm Hg to Pascals

        @note could use the following, but astropy.units.cds is not fully compatible with Python 2
        as of astropy 1.2.1 (see https://github.com/astropy/astropy/issues/5350#issuecomment-248612824):
        astropy.units.cds.mmHg.to(astropy.units.pascal, mmHg)
        """
        return mmHg*PascalPerMmHg

    @staticmethod
    def pascalFromTorr(torr):
        """Convert pressure from torr to Pascals
        """
        return torr*PascalPerTorr
