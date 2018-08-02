#
# LSST Data Management System
# Copyright 2017 AURA/LSST
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import re
import lsst.afw.geom as afwGeom


def bboxFromIraf(irafBBoxStr):
    """Return a Box2I corresponding to an IRAF-style BBOX

    [x0:x1,y0:y1] where x0 and x1 are the one-indexed start and end columns, and correspondingly
    y0 and y1 are the start and end rows.
    """

    mat = re.search(r"^\[([-\d]+):([-\d]+),([-\d]+):([-\d]+)\]$", irafBBoxStr)
    if not mat:
        raise RuntimeError("Unable to parse IRAF-style bbox \"%s\"" % irafBBoxStr)
    x0, x1, y0, y1 = [int(_) for _ in mat.groups()]

    return afwGeom.Box2I(afwGeom.PointI(x0 - 1, y0 - 1), afwGeom.PointI(x1 - 1, y1 - 1), invert=False)
