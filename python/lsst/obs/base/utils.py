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

__all__ = (
    "InitialSkyWcsError",
    "add_provenance_to_fits_header",
    "bboxFromIraf",
    "createInitialSkyWcs",
    "createInitialSkyWcsFromBoresight",
    "strip_provenance_from_fits_header",
)

import re
from collections.abc import MutableMapping

import lsst.geom as geom
import lsst.pex.exceptions
from lsst.afw.cameraGeom import FIELD_ANGLE, PIXELS
from lsst.afw.geom.skyWcs import makeSkyWcs
from lsst.afw.image import RotType
from lsst.daf.base import PropertyList
from lsst.daf.butler import DatasetProvenance, DatasetRef

_PROVENANCE_PREFIX = "LSST BUTLER"
"""The prefix used in FITS headers for butler provenance."""


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


def _store_str_header(
    hdr: PropertyList, key: str, value: str, comment: str | None = None, allow_long_headers: bool = True
) -> None:
    """Examine string header and return value that can be used.

    Parameters
    ----------
    key : `str`
        The key to use in the FITS header (without HIERARCH).
    value : `str`
        The value that is to be stored in the header.
    allow_long_headers : `bool`, optional
        If `True` the value will be used unchanged. If `False` the value
        could have some content elided and a modified version used that will
        fit in a FITS header. If `False` the maximum length of ``key`` is
        58 characters to prevent inconsistency where a very short string will
        fit in 80 characters with ``HIERARCH`` but a very long string will
        not.
    """
    if not allow_long_headers:
        # Declare that we do not allow keys longer than a fixed number of
        # characters in this mode. This is enough to allow
        # HIERARCH {key} = 'abc     '
        # to fit into 80 characters.
        max_card_length = 80
        if len(key) > 58:
            raise ValueError(
                f"Given keyword {key} is too long when requiring {max_card_length}-character header cards"
            )

        # FITS pads strings to 8 characters.
        min_len = 8
        formatted = f"HIERARCH {key} = '{value:{min_len}s}'"
        if (n_over := len(formatted) - max_card_length) > 0:
            # Elide some content.
            # Remove the over run + 3 for the ... and 1 for rounding.
            n_remove = n_over + 1 + 3
            middle = len(value) // 2
            half_remove = n_remove // 2
            value = f"{value[: middle - half_remove]}...{value[middle + half_remove :]}"

            # Do not forward comment if we have elided.
            comment = None

    hdr.set(key, value, comment)


def add_provenance_to_fits_header(
    hdr: PropertyList | MutableMapping | None,
    ref: DatasetRef,
    provenance: DatasetProvenance | None = None,
    allow_long_headers: bool = True,
) -> None:
    """Modify the given header to include provenance headers.

    Parameters
    ----------
    hdr : `lsst.daf.base.PropertyList` or `collections.abc.MutableMapping`
        The FITS header to modify. Assumes ``HIERARCH`` will be handled
        implicitly by the writer.
    ref : `lsst.daf.butler.DatasetRef`
        The butler dataset associated with this FITS file.
    provenance : `lsst.daf.butler.DatasetProvenance` or `None`, optional
        Provenance for this dataset.
    allow_long_headers : `bool`, optional
        If `True` it is assumed that there is no limit to the length of the
        values being stored. If `False`, assumes that FITS header cards must be
        kept within 80 characters including quoting, ``HIERARCH``, and `=`.
    """
    # Some datasets do not define a header.
    if hdr is None:
        return

    if provenance is None:
        provenance = DatasetProvenance()

    # Get the flat dict in form suitable for FITS.
    prov_dict = provenance.to_flat_dict(ref, prefix="LSST BUTLER", sep=" ", simple_types=True)

    # Copy keys into a PropertyList so we have the option of including
    # comments.
    extras = PropertyList()

    for k, v in prov_dict.items():
        if not allow_long_headers and isinstance(v, str):
            _store_str_header(extras, k, v, allow_long_headers=allow_long_headers)
        else:
            extras.set(k, v)
        # Only add comments if we think it's safe.
        if allow_long_headers:
            comment = ""
            if re.search(r"\bRUN\b", k):
                comment = "Run collection"
            elif re.search(r"\bID\b", k):
                comment = "Dataset ID"
            elif re.search(r"\bDATAID\b", k):
                comment = "Data identifier"
            elif re.search(r"\bDATASETTYPE\b", k):
                comment = "Dataset type"
            if comment:
                extras.setComment(k, comment)

    # Purge old headers from metadata (important for data ID and input headers
    # and to prevent headers accumulating in a PropertyList).
    strip_provenance_from_fits_header(hdr)

    # Update the header.
    hdr.update(extras)


def strip_provenance_from_fits_header(hdr: MutableMapping | PropertyList) -> None:
    """Remove all FITS headers relating to butler provenance.

    Parameters
    ----------
    hdr : `lsst.daf.base.PropertyList` or `collections.abc.MutableMapping`
        The FITS header to modify. Assumes ``HIERARCH`` will be handled
        implicitly by the writer.

    Notes
    -----
    These headers will have been added by, for example, the butler formatter
    via a call to `add_provenance_to_fits_header`.
    """
    for k in list(hdr):
        if k.startswith(_PROVENANCE_PREFIX):
            del hdr[k]
