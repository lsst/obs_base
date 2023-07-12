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

"""Classes to allow obs packages to define the filters used by an Instrument
and for use by `lsst.afw.image.Filter`, gen2 dataIds, and gen3 Dimensions.
"""

from __future__ import annotations

__all__ = ("FilterDefinition", "FilterDefinitionCollection")

import dataclasses
from collections.abc import Sequence, Set
from typing import Any, overload

import lsst.afw.image.utils


@dataclasses.dataclass(frozen=True)
class FilterDefinition:
    """The definition of an instrument's filter bandpass.

    This class is used to declare ``physical_filter`` and ``band``
    information for an instrument.

    This class is likely temporary, until we have a better versioned filter
    definition system that includes complete transmission information.
    """

    physical_filter: str
    """The name of a filter associated with a particular instrument: unique for
    each piece of glass. This should match the exact filter name used in the
    observatory's metadata.

    This name is used to define the ``physical_filter`` gen3 Butler Dimension.

    If neither ``band`` or ``afw_name`` is defined, this is used
    as the `~lsst.afw.image.Filter` ``name``, otherwise it is added to the
    list of `~lsst.afw.image.Filter` aliases.
    """

    band: str | None = None
    """The generic name of a filter not associated with a particular instrument
    (e.g. `r` for the SDSS Gunn r-band, which could be on SDSS, LSST, or HSC).

    Not all filters have an abstract filter: engineering or test filters may
    not have a genericly-termed filter name.

    If specified and if `afw_name` is None, this is used as the
    `~lsst.afw.image.Filter` ``name`` field, otherwise it is added to the list
    of `~lsst.afw.image.Filter` aliases.
    """

    doc: str | None = None
    """A short description of this filter, possibly with a link to more
    information.
    """

    afw_name: str | None = None
    """If not None, the name of the `~lsst.afw.image.Filter` object.

    This is distinct from physical_filter and band to maintain
    backwards compatibility in some obs packages.
    For example, for HSC there are two distinct ``r`` and ``i`` filters, named
    ``r/r2`` and ``i/i2``.
    """

    alias: Set[str] = frozenset()
    """Alternate names for this filter. These are added to the
    `~lsst.afw.image.Filter` alias list.
    """

    def __post_init__(self) -> None:
        # force alias to be immutable, so that hashing works
        if not isinstance(self.alias, frozenset):
            object.__setattr__(self, "alias", frozenset(self.alias))

    def __str__(self) -> str:
        txt = f"FilterDefinition(physical_filter='{self.physical_filter}'"
        if self.band is not None:
            txt += f", band='{self.band}'"
        if self.afw_name is not None:
            txt += f", afw_name='{self.afw_name}'"
        if len(self.alias) != 0:
            txt += f", alias='{self.alias}'"
        return txt + ")"

    def makeFilterLabel(self) -> lsst.afw.image.FilterLabel:
        """Create a complete FilterLabel for this filter."""
        return lsst.afw.image.FilterLabel(band=self.band, physical=self.physical_filter)


class FilterDefinitionCollection(Sequence[FilterDefinition]):
    """An order-preserving collection of multiple `FilterDefinition`.

    Parameters
    ----------
    filters : `Sequence`
        The filters in this collection.
    """

    physical_to_band: dict[str, str | None]
    """A mapping from physical filter name to band name.
    This is a convenience feature to allow file readers to create a FilterLabel
    when reading a raw file that only has a physical filter name, without
    iterating over the entire collection.
    """

    def __init__(self, *filters: FilterDefinition):
        self._filters = list(filters)
        self.physical_to_band = {filter.physical_filter: filter.band for filter in self._filters}

    @overload
    def __getitem__(self, i: int) -> FilterDefinition:
        pass

    @overload
    def __getitem__(self, s: slice) -> Sequence[FilterDefinition]:
        pass

    def __getitem__(self, index: Any) -> Any:
        return self._filters[index]

    def __len__(self) -> int:
        return len(self._filters)

    def __str__(self) -> str:
        return "FilterDefinitions(" + ", ".join(str(f) for f in self._filters) + ")"

    def findAll(self, name: str) -> set[FilterDefinition]:
        """Return the FilterDefinitions that match a particular name.

        This method makes no attempt to prioritize, e.g., band names over
        physical filter names; any definition that makes *any* reference
        to the name is returned.

        Parameters
        ----------
        name : `str`
            The name to search for. May be any band, physical, or alias name.

        Returns
        -------
        matches : `set` [`FilterDefinition`]
            All FilterDefinitions containing ``name`` as one of their
            filter names.
        """
        matches = set()
        for filter in self._filters:
            if (
                name == filter.physical_filter
                or name == filter.band
                or name == filter.afw_name
                or name in filter.alias
            ):
                matches.add(filter)
        return matches
