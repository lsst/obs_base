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

__all__ = ("InitialSkyWcsError", "createInitialSkyWcs", "createInitialSkyWcsFromBoresight", "bboxFromIraf")

import re

import lsst.geom as geom
import lsst.pex.exceptions
from lsst.afw.cameraGeom import FIELD_ANGLE, PIXELS
from lsst.afw.geom.skyWcs import makeSkyWcs
from lsst.afw.image import RotType


class InitialSkyWcsError(Exception):
    """For handling failures when creating a SkyWcs from a camera geometry and
    boresight.

    Typically used as a chained exception from a lower level exception.
    """

    pass


def createInitialSkyWcs(visitInfo, detector, flipX=False):
    """Create a SkyWcs from the visit information and detector geometry.

    A typical use case for this is to create the initial WCS for a newly-read
    raw exposure.


    Parameters
    ----------
    visitInfo : `lsst.afw.image.VisitInfo`
        Where to get the telescope boresight and rotator angle from.
    detector : `lsst.afw.cameraGeom.Detector`
        Where to get the camera geometry from.
    flipX : `bool`, optional
        If False, +X is along W, if True +X is along E.

    Returns
    -------
    skyWcs : `lsst.afw.geom.SkyWcs`
        The new composed WCS.

    Raises
    ------
    InitialSkyWcsError
        Raised if there is an error generating the SkyWcs, chained from the
        lower-level exception if available.
    """
    if visitInfo.getRotType() != RotType.SKY:
        msg = (
            "Cannot create SkyWcs from camera geometry: rotator angle defined using "
            f"RotType={visitInfo.getRotType()} instead of SKY."
        )
        raise InitialSkyWcsError(msg)
    orientation = visitInfo.getBoresightRotAngle()
    boresight = visitInfo.getBoresightRaDec()
    return createInitialSkyWcsFromBoresight(boresight, orientation, detector, flipX)


def createInitialSkyWcsFromBoresight(boresight, orientation, detector, flipX=False):
    """Create a SkyWcs from the telescope boresight and detector geometry.

    A typical usecase for this is to create the initial WCS for a newly-read
    raw exposure.

    Parameters
    ----------
    boresight : `lsst.geom.SpherePoint`
        The ICRS boresight RA/Dec
    orientation : `lsst.geom.Angle`
        The rotation angle of the focal plane on the sky.
    detector : `lsst.afw.cameraGeom.Detector`
        Where to get the camera geometry from.
    flipX : `bool`, optional
        If False, +X is along W, if True +X is along E.

    Returns
    -------
    skyWcs : `lsst.afw.geom.SkyWcs`
        The new composed WCS.

    Raises
    ------
    InitialSkyWcsError
        Raised if there is an error generating the SkyWcs, chained from the
        lower-level exception if available.
    """
    try:
        pixelsToFieldAngle = detector.getTransform(
            detector.makeCameraSys(PIXELS), detector.makeCameraSys(FIELD_ANGLE)
        )
    except lsst.pex.exceptions.InvalidParameterError as e:
        raise InitialSkyWcsError("Cannot compute PIXELS to FIELD_ANGLE Transform.") from e
    return makeSkyWcs(pixelsToFieldAngle, orientation, flipX, boresight)


def bboxFromIraf(irafBBoxStr):
    """Return a Box2I corresponding to an IRAF-style BBOX.

    [x0:x1,y0:y1] where x0 and x1 are the one-indexed start and end columns,
    and correspondingly y0 and y1 are the start and end rows.
    """
    mat = re.search(r"^\[([-\d]+):([-\d]+),([-\d]+):([-\d]+)\]$", irafBBoxStr)
    if not mat:
        raise RuntimeError(f'Unable to parse IRAF-style bbox "{irafBBoxStr}"')
    x0, x1, y0, y1 = (int(_) for _ in mat.groups())

    return geom.BoxI(geom.PointI(x0 - 1, y0 - 1), geom.PointI(x1 - 1, y1 - 1))
