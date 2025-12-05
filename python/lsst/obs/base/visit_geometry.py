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

from __future__ import annotations

__all__ = ("VisitGeometry",)


import pydantic

from lsst.daf.butler.pydantic_utils import SerializableRegion


class VisitGeometry(pydantic.BaseModel):
    """A serializable struct that holds geometry information for a visit.

    This is intended to be used as a butler output dataset intermediary between
    a tasks that fit coordinate mappings from stars to reference catalogs and
    tool that update butler dimension record regions.
    """

    boresight_ra: float
    """Re-fit boresight right ascension in degrees."""

    boresight_dec: float
    """Re-fit boresight declination in degrees."""

    orientation: float
    """Re-fit rotation angle in degress."""

    visit_region: SerializableRegion
    """Updated region for the visit."""

    detector_regions: dict[int, SerializableRegion] = pydantic.Field(default_factory=dict)
    """Updated region for each detector in this visit."""
