# This file is part of pipe_tasks.
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

__all__ = ["ColortermModel", "InstrumentRefCatData", "InstrumentRefCatLibrary"]

import fnmatch
from collections.abc import Sequence
from typing import TypeVar

import pydantic

from .filters import FilterDefinition


class ColortermModel(pydantic.BaseModel):
    """Colorterm correction for one pair of filters.

    Notes
    -----
    The transformed magnitude p' is given by::

        p' = primary + c0 + c1*(primary - secondary)
             + c2*(primary - secondary)**2

    """

    primary: str
    """Name of primary filter."""

    secondary: str
    """Name of the secondary filter."""

    c0: float = 0.0
    """Constant parameter."""

    c1: float = 0.0
    """First-order parameter."""

    c2: float = 0.0
    """Second-order parameter."""


class InstrumentRefCatData(pydantic.BaseModel):
    """A data model for a single Instrument's relationship with a reference
    catalog.
    """

    filter_map: str | dict[str, str] = pydantic.Field(default_factory=dict)
    """Mapping from band name to the reference catalog filter nanme, or a
    single reference catalog filter name that corresponds to all bands.
    """

    colorterms: dict[str, ColortermModel] = pydantic.Field(default_factory=dict)
    """Mapping from physical filter to a colorterm on reference catalog filter
    names.
    """


class InstrumentRefCatLibrary(pydantic.RootModel):
    """A mapping from a shell-style glob expression on reference catalog name
    to information about how that reference catalog's filters interact with
    a particular instrument.
    """

    root: dict[str, InstrumentRefCatData]

    def add_gaia(self) -> None:
        """Add an entry for the Gaia catalog, in which every filter maps to
        the broad 'g' band, and there are no colorterms.
        """
        self.root["gaia*"] = InstrumentRefCatData(filter_map="phot_g_mean", colorterms={})

    def add_physical_filter_mappings(self, filter_definitions: Sequence[FilterDefinition]) -> None:
        """Add mappings to each `filter_map` for all physical filters for which
        a band is defined and is already present in the filter map.

        This can be used to support tasks that do not expect filter maps to
        have bands rather than physical filter namess.

        Parameters
        ----------
        filter_definitions : `~collections.abc.Sequence` [ `FilterDefinition` ]
            Iterable of filter definitions.
        """
        for ref_cat_data in self.root.values():
            if isinstance(ref_cat_data.filter_map, str):
                continue
            for filter_definition in filter_definitions:
                if filter_definition.band is None:
                    continue
                ref_cat_filter = ref_cat_data.filter_map.get(filter_definition.band, filter_definition.band)
                ref_cat_data.filter_map.setdefault(filter_definition.physical_filter, ref_cat_filter)

    def find(self, ref_cat_name: str) -> InstrumentRefCatData:
        """Return the data for the given reference catalog.

        Parameters
        ----------
        ref_cat_name : `str`
            Full reference catalog name.

        Returns
        -------
        data : `InstrumentRefCatData`
            Information about an instrument's relationship to the given
            reference catalog.
        """
        return _find_ref_cat_in_library(ref_cat_name, self.root)


_T = TypeVar("_T")


# This is a free function so lsst.pipe.tasks.colorterms can use it until we get
# a chance to more fully unify that code with what's here.
def _find_ref_cat_in_library(ref_cat_name: str, library: dict[str, _T]) -> _T:
    """Find the value in a dictionary whose key is a shell-style glob match to
    a reference catalog name.

    Parameters
    ----------
    ref_cat_name : `str`
        Full reference catalog name.
    library : `dict`
        Dictionary whose keys are reference catalog name glob patterns.

    Returns
    -------
    value
        Matching value.

    Raises
    ------
    LookupError
        Raised if no unambiguous match was found.
    """
    if (data := library.get(ref_cat_name, None)) is not None:
        return data
    # try glob expression
    matched_keys = [key for key in library if fnmatch.fnmatch(ref_cat_name, key)]
    if len(matched_keys) == 1:
        return library[matched_keys[0]]
    elif len(matched_keys) > 1:
        raise LookupError(f"Multiple globs match refcat name {ref_cat_name}: {matched_keys}")
    else:
        raise LookupError(f"No globs match refcat name {ref_cat_name}.")
