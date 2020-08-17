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

__all__ = ('InitialSkyWcsError', 'createInitialSkyWcs', 'createInitialSkyWcsFromBoresight', 'bboxFromIraf')

import re
import lsst.geom as geom

from . import Instrument
from lsst.afw.cameraGeom import PIXELS, FIELD_ANGLE
from lsst.afw.image import RotType
from lsst.afw.geom.skyWcs import makeSkyWcs
import lsst.pex.exceptions
from lsst.utils import doImport


class InitialSkyWcsError(Exception):
    """For handling failures when creating a SkyWcs from a camera geometry and
    boresight.

    Typically used as a chained exception from a lower level exception.
    """
    pass


def createInitialSkyWcs(visitInfo, detector, flipX=False):
    """Create a SkyWcs from the visit information and detector geometry.

    A typical usecase for this is to create the initial WCS for a newly-read
    raw exposure.


    Parameters
    ----------
    visitInfo : `lsst.afw.image.VisitInfo`
        Where to get the telescope boresight and rotator angle from.
    detector : `lsst.afw.cameraGeom.Detector`
        Where to get the camera geomtry from.
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
        msg = (f"Cannot create SkyWcs from camera geometry: rotator angle defined using "
               f"RotType={visitInfo.getRotType()} instead of SKY.")
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
        Where to get the camera geomtry from.
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
        pixelsToFieldAngle = detector.getTransform(detector.makeCameraSys(PIXELS),
                                                   detector.makeCameraSys(FIELD_ANGLE))
    except lsst.pex.exceptions.InvalidParameterError as e:
        raise InitialSkyWcsError("Cannot compute PIXELS to FIELD_ANGLE Transform.") from e
    return makeSkyWcs(pixelsToFieldAngle, orientation, flipX, boresight)


def bboxFromIraf(irafBBoxStr):
    """Return a Box2I corresponding to an IRAF-style BBOX

    [x0:x1,y0:y1] where x0 and x1 are the one-indexed start and end columns, and correspondingly
    y0 and y1 are the start and end rows.
    """

    mat = re.search(r"^\[([-\d]+):([-\d]+),([-\d]+):([-\d]+)\]$", irafBBoxStr)
    if not mat:
        raise RuntimeError("Unable to parse IRAF-style bbox \"%s\"" % irafBBoxStr)
    x0, x1, y0, y1 = [int(_) for _ in mat.groups()]

    return geom.BoxI(geom.PointI(x0 - 1, y0 - 1), geom.PointI(x1 - 1, y1 - 1))


def getInstrument(instrumentName, registry=None):
    """Return an instance of a named instrument.

    If the instrument name not is qualified (does not contain a '.') and a
    butler registry is provided, this will attempt to load the instrument using
    Instrument.fromName. Otherwise the instrument will be imported and
    instantiated.

    Parameters
    ----------
    instrumentName : string
        The name or fully-qualified class name of an instrument.
    registry : `lsst.daf.butler.Registry`, optional
        Butler registry to query to find information about the instrument, by
        default None

    Returns
    -------
    Instrument subclass instance
        The instantiated instrument.

    Raises
    ------
    RuntimeError
        If the instrument can not be imported, instantiated, or obtained from
        the registry.
    TypeError
        If the instrument is not a subclass of lsst.obs.base.Instrument.
    """
    if "." not in instrumentName and registry is not None:
        try:
            instr = Instrument.fromName(instrumentName, registry)
        except Exception as err:
            raise RuntimeError(
                f"Could not get instrument from name: {instrumentName}. Failed with exception: {err}")
    else:
        try:
            instr = doImport(instrumentName)
        except Exception as err:
            raise RuntimeError(f"Could not import instrument: {instrumentName}. Failed with exception: {err}")
        instr = instr()
    if not isinstance(instr, Instrument):
        raise TypeError(f"{instrumentName} is not an Instrument subclass.")
    return instr


#  TODO remove the impl in pipe_base? (NB this combines setDottedAtr AND the
# handling in ConfigValueAction.__call__)
def setDottedAttr(item, name, value):
    """Set an instance attribute (like `setattr` but accepting hierarchical
    names such as ``foo.bar.baz``) If the attribute can not be set as a string,
    will attempt to set the attribute with the result of eval'ing the value.

    Parameters
    ----------
    item : obj
        Object whose attribute is to be set.
    name : `str`
        Name of attribute to set.
    value : obj
        New value for the attribute.

    Notes
    -----
    For example if name is ``foo.bar.baz`` then ``item.foo.bar.baz``
    is set to the specified value.

    Raises
    ------
    AttributeError
        If the item does not have a field specified by name that can be set.
    RuntimeError
        If the value can not be set as a string or rendered by eval, or if
        there is an error setting the attribute with the rendered value.
    """
    subitem = item
    subnameList = name.split(".")
    for subname in subnameList[:-1]:
        subitem = getattr(subitem, subname)
    try:
        setattr(subitem, subnameList[-1], value)
    except AttributeError:
        raise AttributeError(f"No field: {name!r}")
    except Exception:
        try:
            v = eval(value, {})
        except Exception:
            raise RuntimeError(f"Cannot render {value!r} as a value for {name!r}")
        try:
            setattr(subitem, subnameList[-1], v)
        except Exception as e:
            raise RuntimeError(f"Cannot set config. {name}={value!r}: {e}")


def setDottedAttrs(item, attrs):
    for name, value in attrs:
        setDottedAttr(item, name, value)
