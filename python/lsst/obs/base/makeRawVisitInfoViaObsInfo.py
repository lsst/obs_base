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

__all__ = ["MakeRawVisitInfoViaObsInfo"]

import logging
import warnings
from typing import ClassVar, Optional

import astropy.units
import astropy.utils.exceptions
from astropy.utils import iers

# Prefer the standard pyerfa over the Astropy version
try:
    import erfa

    ErfaWarning = erfa.ErfaWarning
except ImportError:
    import astropy._erfa as erfa

    ErfaWarning = None

from astro_metadata_translator import MetadataTranslator, ObservationInfo
from lsst.afw.coord import Observatory, Weather
from lsst.afw.image import RotType, VisitInfo
from lsst.daf.base import DateTime
from lsst.geom import SpherePoint, degrees, radians


class MakeRawVisitInfoViaObsInfo:
    """Base class functor to make a VisitInfo from the FITS header of a
    raw image using `~astro_metadata_translator.ObservationInfo` translators.

    Subclasses can be used if a specific
    `~astro_metadata_translator.MetadataTranslator` translator should be used.

    The design philosophy is to make a best effort and log warnings of
    problems, rather than raising exceptions, in order to extract as much
    VisitInfo information as possible from a messy FITS header without the
    user needing to add a lot of error handling.

    Parameters
    ----------
    log : `logging.Logger` or None
        Logger to use for messages.
        (None to use
        ``Log.getLogger("lsst.obs.base.makeRawVisitInfoViaObsInfo")``).
    doStripHeader : `bool`, optional
        Strip header keywords from the metadata as they are used?
    """

    metadataTranslator: ClassVar[Optional[MetadataTranslator]] = None
    """Header translator to use to construct VisitInfo, defaulting to
    automatic determination."""

    def __init__(self, log=None, doStripHeader=False):
        if log is None:
            log = logging.getLogger(__name__)
        self.log = log
        self.doStripHeader = doStripHeader

    def __call__(self, md, exposureId=None):
        """Construct a VisitInfo and strip associated data from the metadata.

        Parameters
        ----------
        md : `lsst.daf.base.PropertyList` or `lsst.daf.base.PropertySet`
            Metadata to pull from.
            May be modified if ``stripHeader`` is ``True``.
        exposureId : `int`, optional
            Ignored.  Here for compatibility with `MakeRawVisitInfo`.

        Returns
        -------
        visitInfo : `lsst.afw.image.VisitInfo`
            `~lsst.afw.image.VisitInfo` derived from the header using
            a `~astro_metadata_translator.MetadataTranslator`.
        """

        obsInfo = ObservationInfo(md, translator_class=self.metadataTranslator)

        if self.doStripHeader:
            # Strip all the cards out that were used
            for c in obsInfo.cards_used:
                del md[c]

        return self.observationInfo2visitInfo(obsInfo, log=self.log)

    @staticmethod
    def observationInfo2visitInfo(obsInfo, log=None):
        """Construct a `~lsst.afw.image.VisitInfo` from an
        `~astro_metadata_translator.ObservationInfo`

        Parameters
        ----------
        obsInfo : `astro_metadata_translator.ObservationInfo`
            Information gathered from the observation metadata.
        log : `logging.Logger` or `lsst.log.Log`, optional
            Logger to use for logging informational messages.
            If `None` logging will be disabled.

        Returns
        -------
        visitInfo : `lsst.afw.image.VisitInfo`
            `~lsst.afw.image.VisitInfo` derived from the supplied
            `~astro_metadata_translator.ObservationInfo`.
        """
        argDict = dict()

        # Map the translated information into a form suitable for VisitInfo
        if obsInfo.exposure_time is not None:
            argDict["exposureTime"] = obsInfo.exposure_time.to_value("s")
        if obsInfo.dark_time is not None:
            argDict["darkTime"] = obsInfo.dark_time.to_value("s")
        argDict["exposureId"] = obsInfo.detector_exposure_id
        argDict["id"] = obsInfo.exposure_id
        argDict["instrumentLabel"] = obsInfo.instrument

        # VisitInfo uses the middle of the observation for the date
        if obsInfo.datetime_begin is not None and obsInfo.datetime_end is not None:
            tdelta = obsInfo.datetime_end - obsInfo.datetime_begin
            middle = obsInfo.datetime_begin + 0.5 * tdelta

            # DateTime uses nanosecond resolution, regardless of the resolution
            # of the original date
            middle.precision = 9
            # isot is ISO8601 format with "T" separating date and time and no
            # time zone
            argDict["date"] = DateTime(middle.tai.isot, DateTime.TAI)

            # Derive earth rotation angle from UT1 (being out by a second is
            # not a big deal given the uncertainty over exactly what part of
            # the observation we are needing it for).
            # ERFA needs a UT1 time split into two floats
            # We ignore any problems with DUT1 not being defined for now.
            try:
                # Catch any warnings about the time being in the future
                # since there is nothing we can do about that for simulated
                # data and it tells us nothing for data from the past.
                with warnings.catch_warnings():
                    # If we are using the real erfa it is not an AstropyWarning
                    # During transition period filter both
                    warnings.simplefilter("ignore", category=astropy.utils.exceptions.AstropyWarning)
                    if ErfaWarning is not None:
                        warnings.simplefilter("ignore", category=ErfaWarning)
                    ut1time = middle.ut1
            except iers.IERSRangeError:
                ut1time = middle

            era = erfa.era00(ut1time.jd1, ut1time.jd2)
            argDict["era"] = era * radians
        else:
            argDict["date"] = DateTime()

        # Coordinates
        if obsInfo.tracking_radec is not None:
            icrs = obsInfo.tracking_radec.transform_to("icrs")
            argDict["boresightRaDec"] = SpherePoint(icrs.ra.degree, icrs.dec.degree, units=degrees)

        altaz = obsInfo.altaz_begin
        if altaz is not None:
            argDict["boresightAzAlt"] = SpherePoint(altaz.az.degree, altaz.alt.degree, units=degrees)

        argDict["boresightAirmass"] = obsInfo.boresight_airmass

        if obsInfo.boresight_rotation_angle is not None:
            argDict["boresightRotAngle"] = obsInfo.boresight_rotation_angle.degree * degrees

        if obsInfo.boresight_rotation_coord is not None:
            rotType = RotType.UNKNOWN
            if obsInfo.boresight_rotation_coord == "sky":
                rotType = RotType.SKY
            argDict["rotType"] = rotType

        # Weather and Observatory Location
        temperature = float("nan")
        if obsInfo.temperature is not None:
            temperature = obsInfo.temperature.to_value("deg_C", astropy.units.temperature())
        pressure = float("nan")
        if obsInfo.pressure is not None:
            pressure = obsInfo.pressure.to_value("Pa")
        relative_humidity = float("nan")
        if obsInfo.relative_humidity is not None:
            relative_humidity = obsInfo.relative_humidity
        argDict["weather"] = Weather(temperature, pressure, relative_humidity)

        if obsInfo.location is not None:
            geolocation = obsInfo.location.to_geodetic()
            argDict["observatory"] = Observatory(
                geolocation.lon.degree * degrees,
                geolocation.lat.degree * degrees,
                geolocation.height.to_value("m"),
            )

        for key in list(argDict.keys()):  # use a copy because we may delete items
            if argDict[key] is None:
                if log is not None:
                    log.warning("argDict[%s] is None; stripping", key)
                del argDict[key]

        return VisitInfo(**argDict)
