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

import math
import numpy as np

import astropy.coordinates
import astropy.time
import astropy.units

from lsst.log import Log
from lsst.daf.base import DateTime
from lsst.geom import degrees
from lsst.afw.image import VisitInfo

__all__ = ["MakeRawVisitInfo"]


PascalPerMillibar = 100.0
PascalPerMmHg = 133.322387415  # from Wikipedia; exact
PascalPerTorr = 101325.0/760.0  # from Wikipedia; exact
KelvinMinusCentigrade = 273.15  # from Wikipedia; exact

# have these read at need, to avoid unexpected errors later
NaN = float("nan")
BadDate = DateTime()


class MakeRawVisitInfo(object):
    """Base class functor to make a VisitInfo from the FITS header of a raw
    image.

    A subclass will be wanted for each camera. Subclasses should override:

    - `setArgDict`, The override can call the base implementation,
        which simply sets exposure time and date of observation
    - `getDateAvg`

    The design philosophy is to make a best effort and log warnings of
    problems, rather than raising exceptions, in order to extract as much
    VisitInfo information as possible from a messy FITS header without the
    user needing to add a lot of error handling.

    However, the methods that transform units are less forgiving; they assume
    the user provides proper data types, since type errors in arguments to
    those are almost certainly due to coding mistakes.

    Parameters
    ----------
    log : `lsst.log.Log` or None
        Logger to use for messages.
        (None to use ``Log.getLogger("MakeRawVisitInfo")``).
    doStripHeader : `bool`, optional
        Strip header keywords from the metadata as they are used?
    """

    def __init__(self, log=None, doStripHeader=False):
        if log is None:
            log = Log.getLogger("MakeRawVisitInfo")
        self.log = log
        self.doStripHeader = doStripHeader

    def __call__(self, md, exposureId):
        """Construct a VisitInfo and strip associated data from the metadata.

        Parameters
        ----------
        md : `lsst.daf.base.PropertyList` or `lsst.daf.base.PropertySet`
            Metadata to pull from.
            Items that are used are stripped from the metadata (except TIMESYS,
            because it may apply to other keywords) if ``doStripHeader``.
        exposureId : `int`
            exposure ID

        Notes
        -----
        The basic implementation sets `date` and `exposureTime` using typical
        values found in FITS files and logs a warning if neither can be set.
        """
        argDict = dict(exposureId=exposureId)
        self.setArgDict(md, argDict)
        for key in list(argDict.keys()):  # use a copy because we may delete items
            if argDict[key] is None:
                self.log.warn("argDict[%s] is None; stripping", key)
                del argDict[key]
        return VisitInfo(**argDict)

    def setArgDict(self, md, argDict):
        """Fill an argument dict with arguments for VisitInfo and pop
        associated metadata

        Subclasses are expected to override this method, though the override
        may wish to call this default implementation, which:

        - sets exposureTime from "EXPTIME"
        - sets date by calling getDateAvg

        Parameters
        ----------
        md : `lsst.daf.base.PropertyList` or `PropertySet`
            Metadata to pull from.
            Items that are used are stripped from the metadata (except TIMESYS,
            because it may apply to other keywords).
        argdict : `dict`
            dict of arguments

        Notes
        -----
        Subclasses should expand this or replace it.
        """
        argDict["exposureTime"] = self.popFloat(md, "EXPTIME")
        argDict["date"] = self.getDateAvg(md=md, exposureTime=argDict["exposureTime"])

    def getDateAvg(self, md, exposureTime):
        """Return date at the middle of the exposure.

        Parameters
        ----------
        md : `lsst.daf.base.PropertyList` or `PropertySet`
            Metadata to pull from.
            Items that are used are stripped from the metadata (except TIMESYS,
            because it may apply to other keywords).
        exposureTime : `float`
            Exposure time (sec)

        Notes
        -----
        Subclasses must override. Here is a typical implementation::

            dateObs = self.popIsoDate(md, "DATE-OBS")
            return self.offsetDate(dateObs, 0.5*exposureTime)
        """
        raise NotImplementedError()

    def getDarkTime(self, argDict):
        """Get the darkTime from the DARKTIME keyword, else expTime, else NaN,

        If dark time is available then subclasses should call this method by
        putting the following in their `__init__` method::

            argDict['darkTime'] = self.getDarkTime(argDict)

        Parameters
        ----------
        argdict : `dict`
            Dict of arguments.

        Returns
        -------
        `float`
            Dark time, as inferred from the metadata.
        """
        darkTime = argDict.get("darkTime", NaN)
        if np.isfinite(darkTime):
            return darkTime

        self.log.info("darkTime is NaN/Inf; using exposureTime")
        exposureTime = argDict.get("exposureTime", NaN)
        if not np.isfinite(exposureTime):
            raise RuntimeError("Tried to substitute exposureTime for darkTime but it is not available")

        return exposureTime

    def offsetDate(self, date, offsetSec):
        """Return a date offset by a specified number of seconds.

        date : `lsst.daf.base.DateTime`
            Date baseline to offset from.
        offsetSec : `float`
            Offset, in seconds.

        Returns
        -------
        `lsst.daf.base.DateTime`
            The offset date.
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
        """Return an item of metadata.

        The item is removed if ``doStripHeader`` is ``True``.

        Log a warning if the key is not found.

        Parameters
        ----------
        md : `lsst.daf.base.PropertyList` or `PropertySet`
            Metadata to pull `key` from and (optionally) remove.
        key : `str`
            Metadata key to extract.
        default : `object`
            Value to return if key not found.

        Returns
        -------
        `object`
            The value of the specified key, using whatever type
            md.getScalar(key) returns.
        """
        try:
            if not md.exists(key):
                self.log.warn("Key=\"{}\" not in metadata".format(key))
                return default
            val = md.getScalar(key)
            if self.doStripHeader:
                md.remove(key)
            return val
        except Exception as e:
            # this should never happen, but is a last ditch attempt to avoid
            # exceptions
            self.log.warn('Could not read key="{}" in metadata: {}'.format(key, e))
        return default

    def popFloat(self, md, key):
        """Pop a float with a default of NaN.

        Parameters
        ----------
        md : `lsst.daf.base.PropertyList` or `PropertySet`
            Metadata to pull `key` from.
        key : `str`
            Key to read.

        Returns
        -------
        `float`
            Value of the requested key as a float; float("nan") if the key is
            not found.
        """
        val = self.popItem(md, key, default=NaN)
        try:
            return float(val)
        except Exception as e:
            self.log.warn("Could not interpret {} value {} as a float: {}".format(key, repr(val), e))
        return NaN

    def popAngle(self, md, key, units=astropy.units.deg):
        """Pop an lsst.afw.geom.Angle, whose metadata is in the specified
        units, with a default of Nan

        The angle may be specified as a float or sexagesimal string with 1-3
        fields.

        Parameters
        ----------
        md : `lsst.daf.base.PropertyList` or `PropertySet`
            Metadata to pull `key` from.
        key : `str`
            Key to read from md.

        Returns
        -------
        `lsst.afw.geom.Angle`
            Value of the requested key as an angle; Angle(NaN) if the key is
            not found.
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

        Parameters
        ----------
        md : `lsst.daf.base.PropertyList` or `PropertySet`
            Metadata to pull `key` from.
        key : `str`
            Date key to read from md.
        timesys : `str`
            Time system as a string (not case sensitive), e.g. "UTC" or None;
            if None then look for TIMESYS (but do NOT pop it, since it may be
            used for more than one date) and if not found, use UTC.

        Returns
        -------
        `lsst.daf.base.DateTime`
            Value of the requested date; `DateTime()` if the key is not found.
        """
        isoDateStr = self.popItem(md=md, key=key)
        if isoDateStr is not None:
            try:
                if timesys is None:
                    timesys = md.getScalar("TIMESYS") if md.exists("TIMESYS") else "UTC"
                if isoDateStr.endswith("Z"):  # illegal in FITS
                    isoDateStr = isoDateStr[0:-1]
                astropyTime = astropy.time.Time(isoDateStr, scale=timesys.lower(), format="fits")
                # DateTime uses nanosecond resolution, regardless of the
                # resolution of the original date
                astropyTime.precision = 9
                # isot is ISO8601 format with "T" separating date and time and
                # no time zone
                return DateTime(astropyTime.tai.isot, DateTime.TAI)
            except Exception as e:
                self.log.warn("Could not parse {} = {} as an ISO date: {}".format(key, isoDateStr, e))
        return BadDate

    def popMjdDate(self, md, key, timesys=None):
        """Get a FITS MJD date as an ``lsst.daf.base.DateTime``.

        Parameters
        ----------
        md : `lsst.daf.base.PropertyList` or `PropertySet`
            Metadata to pull `key` from.
        key : `str`
            Date key to read from md.
        timesys : `str`
            Time system as a string (not case sensitive), e.g. "UTC" or None;
            if None then look for TIMESYS (but do NOT pop it, since it may be
            used for more than one date) and if not found, use UTC.

        Returns
        -------
        `lsst.daf.base.DateTime`
            Value of the requested date; `DateTime()` if the key is not found.
        """
        mjdDate = self.popFloat(md, key)
        try:
            if timesys is None:
                timesys = md.getScalar("TIMESYS") if md.exists("TIMESYS") else "UTC"
            astropyTime = astropy.time.Time(mjdDate, format="mjd", scale=timesys.lower())
            # DateTime uses nanosecond resolution, regardless of the resolution
            # of the original date
            astropyTime.precision = 9
            # isot is ISO8601 format with "T" separating date and time and no
            # time zone
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

        Notes
        -----
        Could use the following, but astropy.units.cds is not fully compatible
        with Python 2 as of astropy 1.2.1 (see
        https://github.com/astropy/astropy/issues/5350#issuecomment-248612824):
        astropy.units.cds.mmHg.to(astropy.units.pascal, mmHg)
        """
        return mmHg*PascalPerMmHg

    @staticmethod
    def pascalFromTorr(torr):
        """Convert pressure from torr to Pascals
        """
        return torr*PascalPerTorr

    @staticmethod
    def defaultMetadata(value, defaultValue, minimum=None, maximum=None):
        """Return the value if it is not NaN and within min/max, otherwise
        return defaultValue.

        Parameters
        ----------
        value : `float`
            metadata value returned by popItem, popFloat, or popAngle
        defaultValue : `float``
            default value to use if the metadata value is invalid
        minimum : `float`
            Minimum possible valid value, optional
        maximum : `float`
            Maximum possible valid value, optional

        Returns
        -------
        `float`
            The "validated" value.
        """
        if np.isnan(value):
            retVal = defaultValue
        else:
            if minimum is not None and value < minimum:
                retVal = defaultValue
            elif maximum is not None and value > maximum:
                retVal = defaultValue
            else:
                retVal = value
        return retVal
