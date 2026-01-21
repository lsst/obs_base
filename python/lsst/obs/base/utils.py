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
    "TableVStack",
    "add_provenance_to_fits_header",
    "bboxFromIraf",
    "createInitialSkyWcs",
    "createInitialSkyWcsFromBoresight",
    "iso_date_to_curated_calib_file_root",
    "strip_provenance_from_fits_header",
)

import datetime
import re
from collections.abc import Iterable, MutableMapping
from typing import Any, Self

import astropy.table
import numpy as np
from astro_metadata_translator.headers import merge_headers
from numpy.typing import ArrayLike

import lsst.geom as geom
import lsst.pex.exceptions
from lsst.afw.cameraGeom import FIELD_ANGLE, PIXELS
from lsst.afw.geom.skyWcs import makeSkyWcs
from lsst.afw.image import RotType
from lsst.daf.base import PropertyList
from lsst.daf.butler import DatasetProvenance, DatasetRef

_PROVENANCE_PREFIX = "LSST BUTLER"
"""The prefix used in FITS headers for butler provenance."""

DeferredHandles = Iterable[lsst.daf.butler.DeferredDatasetHandle]


class InitialSkyWcsError(Exception):
    """For handling failures when creating a SkyWcs from a camera geometry and
    boresight.

    Typically used as a chained exception from a lower level exception.
    """

    pass


def createInitialSkyWcs(visitInfo, detector, flipX=None):
    """Create a SkyWcs from the visit information and detector geometry.

    A typical use case for this is to create the initial WCS for a newly-read
    raw exposure.

    Parameters
    ----------
    visitInfo : `lsst.afw.image.VisitInfo`
        Where to get the telescope boresight and rotator angle from.
    detector : `lsst.afw.cameraGeom.Detector`
        Where to get the camera geometry from.
    flipX : `bool` or `None`, optional
        If `False`, +X is along W, if `True` +X is along E.  If `None`, the
        focal plane parity in the camera geometry is assumed to be correct,
        which may not be true for old camera datasets.

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


def createInitialSkyWcsFromBoresight(boresight, orientation, detector, flipX=None):
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
    flipX : `bool` or `None`, optional
        If `False`, +X is along W, if `True` +X is along E.  If `None`, the
        focal plane parity in the camera geometry is assumed to be correct,
        which may not be true for old camera datasets.

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
    camera_parity = detector.getTransformMap().getFocalPlaneParity()
    actual_flip_x = False
    if flipX and not camera_parity:
        # This is probably an old camera definition from before those could
        # hold parity flips.  We need to work around this since data
        # repositories aren't necessarily updated in sync with code.
        actual_flip_x = True
    elif flipX is not None and not flipX and camera_parity:
        raise InitialSkyWcsError(
            "Camera geometry reports a parity flip between FOCAL_PLANE and FIELD_ANGLE, but this "
            "was not expected by caller for this instrument."
        )
    try:
        pixelsToFieldAngle = detector.getTransform(
            detector.makeCameraSys(PIXELS), detector.makeCameraSys(FIELD_ANGLE)
        )
    except lsst.pex.exceptions.InvalidParameterError as e:
        raise InitialSkyWcsError("Cannot compute PIXELS to FIELD_ANGLE Transform.") from e
    return makeSkyWcs(pixelsToFieldAngle, orientation, actual_flip_x, boresight)


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
    hdr : `lsst.daf.base.PropertyList`
        Header to examine.
    key : `str`
        The key to use in the FITS header (without HIERARCH).
    value : `str`
        The value that is to be stored in the header.
    comment : `str` or `None`, optional
        Possible comment to add.
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
        kept within 80 characters including quoting, ``HIERARCH``, and ``=``.
    """
    # Some datasets do not define a header.
    if hdr is None:
        return

    if provenance is None:
        provenance = DatasetProvenance()

    # Get the flat dict in form suitable for FITS. Protect against abnormally
    # large number of inputs.
    prov_dict = provenance.to_flat_dict(
        ref, prefix="LSST BUTLER", sep=" ", simple_types=True, max_inputs=3_000
    )

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


def iso_date_to_curated_calib_file_root(valid_start: str) -> str:
    """Parse an ISO time, potentially from a calibration valid start time,
    and convert to form suitable for use in a curated calibration filename.

    Parameters
    ----------
    valid_start : `str`
        Validity start time in ISOT format.

    Returns
    -------
    date_str : `str`
        Form suitable for use in a curated calibration file name in a
        data package: YYYYMMDDTHHMMSS.
    """
    # Parse the ISO string to ensure conformity if there is an edit to
    # valid_start.
    valid_date = datetime.datetime.fromisoformat(valid_start)

    # Drop any fractional seconds.
    return re.sub(r"\W", "", valid_date.isoformat(timespec="seconds"))


class TableVStack:
    """A helper class for stacking astropy tables without having them all in
    memory at once.

    Parameters
    ----------
    capacity : `int`
        Full size of the final table.

    Notes
    -----
    Unlike `astropy.table.vstack`, this class requires all tables to have the
    exact same columns (it's slightly more strict than even the
    ``join_type="exact"`` argument to `astropy.table.vstack`).
    """

    def __init__(self, capacity: int) -> None:
        self.index: int = 0
        self.capacity = capacity
        self.result: astropy.table.Table | None = None

    @classmethod
    def set_extra_values(
        cls,
        table: astropy.table.Table,
        key: str,
        values: Any,
        capacity: int,
        slicer: slice | None = None,
        validate: bool = True,
    ) -> None:
        """Set extra column values in a slice of a table.

        Parameters
        ----------
        table : `astropy.table.Table`
            The table to set values for.
        key : `str`
            The column key.
        values : `Any`
            The value(s) to set. Can be a scalar.
        capacity : `int`
            The size to initialize the column with, if it doesn't exist yet.
        slicer : `slice` or `None`, optional
            A slice to select values to update.
        validate : `bool`, optional
            If True and the column already exists, will raise if the new values
            do not match the existing ones.
        """
        if key in table.colnames:
            column = table[key]
            if validate and not np.all((table[key] if (slicer is None) else table[key][slicer]) == values):
                raise RuntimeError(
                    f"table already contains {column=} with {key=} but values differ from {values=}"
                )
            table[key][slicer] = values
        else:
            if slicer is None:
                table[key] = values
            else:
                try:
                    dtype = values.dtype
                except AttributeError:
                    dtype = np.dtype(type(values))
                table[key] = np.empty(capacity, dtype=dtype)
                table[key][slicer] = values

    @classmethod
    def from_handles(cls, handles: DeferredHandles) -> Self:
        """Construct from an iterable of
          `lsst.daf.butler.DeferredDatasetHandle`.

        Parameters
        ----------
        handles : `~collections.abc.Iterable` [ \
                `lsst.daf.butler.DeferredDatasetHandle` ]
            Iterable of handles.   Must have a storage class that supports the
            "rowcount" component, which is all that will be fetched.

        Returns
        -------
        vstack : `TableVStack`
            An instance of this class, initialized with capacity equal to the
            sum of the rowcounts of all the given table handles.
        """
        capacity = sum(handle.get(component="rowcount") for handle in handles)
        return cls(capacity=capacity)

    def extend(self, table: astropy.table.Table, extra_values: dict[str, ArrayLike] | None = None) -> None:
        """Add a single table to the stack.

        Parameters
        ----------
        table : `astropy.table.Table`
            An astropy table instance.
        extra_values : `dict`
            Dict keyed by column name with an array-like of values to set
            for this table only.
        """
        if extra_values is None:
            extra_values = {}
        if self.result is None:
            self.result = astropy.table.Table()
            slicer = slice(None, len(table))
            for name in table.colnames:
                column = table[name]
                column_cls = type(column)
                self.result[name] = column_cls.info.new_like([column], self.capacity, name=name)
                self.result[name][: len(table)] = column
            for name, values in extra_values.items():
                self.set_extra_values(
                    table=self.result,
                    key=name,
                    values=values,
                    capacity=self.capacity,
                    slicer=slicer,
                )
            self.index = len(table)
            self.result.meta = table.meta.copy()
        else:
            next_index = self.index + len(table)
            slicer = slice(self.index, next_index)
            for name in table.colnames:
                out_col = self.result[name]
                in_col = table[name]
                if out_col.dtype != in_col.dtype:
                    raise TypeError(f"Type mismatch on column {name!r}: {out_col.dtype} != {in_col.dtype}.")
                self.result[name][slicer] = table[name]
            for name, values in extra_values.items():
                self.set_extra_values(
                    table=self.result,
                    key=name,
                    values=values,
                    capacity=self.capacity,
                    slicer=slicer,
                    validate=False,
                )
            self.index = next_index
            # Butler provenance should be stripped on merge. It will be
            # added by butler on write. No attempt is made here to combine
            # provenance from multiple input tables.
            self.result.meta = merge_headers([self.result.meta, table.meta], mode="drop")
            strip_provenance_from_fits_header(self.result.meta)

    @classmethod
    def vstack_handles(
        cls,
        handles: DeferredHandles,
        extra_values: dict[int, dict[str, ArrayLike]] | None = None,
        kwargs_get: dict[str, Any] | None = None,
    ) -> astropy.table.Table:
        """Vertically stack tables represented by deferred dataset handles.

        Parameters
        ----------
        handles : `~collections.abc.Iterable` [ \
                `lsst.daf.butler.DeferredDatasetHandle` ]
            Iterable of handles.   Must have the "ArrowAstropy" storage class
            and identical columns.
        extra_values : `dict`
            Dictionary keyed by index of handle of additional values
            to pass to extend.
        kwargs_get : `dict` [`str`, `Any`]
            Keyword argument-value pairs to pass to handle.get().

        Returns
        -------
        table : `astropy.table.Table`
            Concatenated table with the same columns as each input table and
            the rows of all of them, or an empty table if there are no handles.
        """
        if not handles:
            return astropy.table.Table()
        if extra_values is None:
            extra_values = {}
        if kwargs_get is None:
            kwargs_get = {}
        handles = tuple(handles)  # guard against single-pass iterators
        # Ensure that zero length catalogs are not included
        rowcounts = tuple(handle.get(component="rowcount") for handle in handles)
        handles = tuple(handle for handle, count in zip(handles, rowcounts) if count > 0)

        vstack = cls(capacity=np.sum(rowcounts))
        for idx, handle in enumerate(handles):
            vstack.extend(handle.get(**kwargs_get), extra_values=extra_values.get(idx))
        return astropy.table.Table() if vstack.result is None else vstack.result
