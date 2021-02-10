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

__all__ = ("FilterDefinition", "FilterDefinitionCollection")

import dataclasses
import collections.abc
import re
import warnings

import numpy as np

import lsst.afw.image.utils


class FilterDefinitionCollection(collections.abc.Sequence):
    """An order-preserving collection of multiple `FilterDefinition`.

    Parameters
    ----------
    filters : `~collections.abc.Sequence`
        The filters in this collection.
    """

    _defined = None
    """Whether these filters have been defined via
    `~lsst.afw.image.utils.defineFilter`. If so, set to ``self`` to identify
    the filter collection that defined them.
    """

    physical_to_band = {}
    """A mapping from physical filter name to band name.
    This is a convenience feature to allow file readers to create a FilterLabel
    when reading a raw file that only has a physical filter name, without
    iterating over the entire collection.
    """

    def __init__(self, *filters):
        self._filters = list(filters)
        self.physical_to_band = {filter.physical_filter: filter.band for filter in self._filters}

    def __getitem__(self, key):
        return self._filters[key]

    def __len__(self):
        return len(self._filters)

    def __str__(self):
        return "FilterDefinitions(" + ', '.join(str(f) for f in self._filters) + ')'

    def defineFilters(self):
        """Define all the filters to `lsst.afw.image.Filter`.

        `~lsst.afw.image.Filter` objects are singletons, so we protect against
        filters being defined multiple times.

        Raises
        ------
        RuntimeError
            Raised if any other `FilterDefinitionCollection` has already called
            ``defineFilters``.
        """
        if self._defined is None:
            with warnings.catch_warnings():
                # surpress Filter warnings; we already know this is deprecated
                warnings.simplefilter('ignore', category=FutureWarning)
                self.reset()
                for filter in self._filters:
                    filter.defineFilter()
            FilterDefinitionCollection._defined = self
        elif self._defined is self:
            # noop: we've already defined these filters, so do nothing
            pass
        else:
            msg = f"afw Filters were already defined on: {self._defined}"
            raise RuntimeError(msg)

    @classmethod
    def reset(cls):
        """Reset the afw Filter definitions and clear the `defined` singleton.
        Use this in unittests that define different filters.
        """
        with warnings.catch_warnings():
            # surpress Filter warnings; we already know this is deprecated
            warnings.simplefilter('ignore', category=FutureWarning)
            lsst.afw.image.utils.resetFilters()
        cls._defined = None

    def findAll(self, name):
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
            if name == filter.physical_filter or name == filter.band or name == filter.afw_name \
                    or name in filter.alias:
                matches.add(filter)
        return matches


@dataclasses.dataclass(frozen=True)
class FilterDefinition:
    """The definition of an instrument's filter bandpass.

    This class is used to interface between the `~lsst.afw.image.Filter` class
    and the Gen2 `~lsst.daf.persistence.CameraMapper` and Gen3
    `~lsst.obs.base.Instruments` and ``physical_filter``/``band``
    `~lsst.daf.butler.Dimension`.

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

    lambdaEff: float
    """The effective wavelength of this filter (nm)."""

    band: str = None
    """The generic name of a filter not associated with a particular instrument
    (e.g. `r` for the SDSS Gunn r-band, which could be on SDSS, LSST, or HSC).

    Not all filters have an abstract filter: engineering or test filters may
    not have a genericly-termed filter name.

    If specified and if `afw_name` is None, this is used as the
    `~lsst.afw.image.Filter` ``name`` field, otherwise it is added to the list
    of `~lsst.afw.image.Filter` aliases.
    """

    doc: str = None
    """A short description of this filter, possibly with a link to more
    information.
    """

    afw_name: str = None
    """If not None, the name of the `~lsst.afw.image.Filter` object.

    This is distinct from physical_filter and band to maintain
    backwards compatibility in some obs packages.
    For example, for HSC there are two distinct ``r`` and ``i`` filters, named
    ``r/r2`` and ``i/i2``.
    """

    lambdaMin: float = np.nan
    """The minimum wavelength of this filter (nm; defined as 1% throughput)"""
    lambdaMax: float = np.nan
    """The maximum wavelength of this filter (nm; defined as 1% throughput)"""

    alias: set = frozenset()
    """Alternate names for this filter. These are added to the
    `~lsst.afw.image.Filter` alias list.
    """

    def __post_init__(self):
        # force alias to be immutable, so that hashing works
        if not isinstance(self.alias, frozenset):
            object.__setattr__(self, 'alias', frozenset(self.alias))

    def __str__(self):
        txt = f"FilterDefinition(physical_filter='{self.physical_filter}', lambdaEff='{self.lambdaEff}'"
        if self.band is not None:
            txt += f", band='{self.band}'"
        if self.afw_name is not None:
            txt += f", afw_name='{self.afw_name}'"
        if not np.isnan(self.lambdaMin):
            txt += f", lambdaMin='{self.lambdaMin}'"
        if not np.isnan(self.lambdaMax):
            txt += f", lambdaMax='{self.lambdaMax}'"
        if len(self.alias) != 0:
            txt += f", alias='{self.alias}'"
        return txt + ")"

    def defineFilter(self):
        """Declare the filters via afw.image.Filter.
        """
        aliases = set(self.alias)
        name = self.physical_filter

        current_names = set(lsst.afw.image.Filter.getNames())
        # band can be defined multiple times -- only use the first
        # occurrence
        band = self.band
        if band is not None:
            # Special case code that uses afw_name to override the band.
            # This was generally used as a workaround.
            if self.afw_name is not None and (mat := re.match(fr"{band}(\d)+$", self.afw_name)):
                i = int(mat.group(1))
            else:
                i = 0
            while i < 50:  # in some instruments gratings or ND filters are combined
                if i == 0:
                    nband = band
                else:
                    nband = f"{band}{i}"
                if nband not in current_names:
                    band = nband
                    name = band
                    aliases.add(self.physical_filter)
                    break
                i += 1
            else:
                warnings.warn(f"Too many band aliases found for physical_filter {self.physical_filter}"
                              f" with band {band}")

        if self.physical_filter == self.band and self.physical_filter in current_names:
            # We have already defined a filter like this
            return

        # Do not add an alias for band if the given afw_name matches
        # the dynamically calculated band.
        if self.afw_name is not None and self.afw_name != band and self.afw_name not in current_names:
            # This will override the band setting above but it
            # is still used as an alias below
            name = self.afw_name
            aliases.add(self.physical_filter)

            # Only add physical_filter/band as an alias if afw_name is defined.
            if band is not None:
                aliases.add(band)

        # Aliases are a serious issue so as a last attempt to clean up
        # remove any registered names from the new aliases
        # This usually means some variant filter name is being used
        aliases.difference_update(current_names)

        with warnings.catch_warnings():
            # surpress Filter warnings; we already know this is deprecated
            warnings.simplefilter('ignore', category=FutureWarning)
            lsst.afw.image.utils.defineFilter(name,
                                              lambdaEff=self.lambdaEff,
                                              lambdaMin=self.lambdaMin,
                                              lambdaMax=self.lambdaMax,
                                              alias=sorted(aliases))

    def makeFilterLabel(self):
        """Create a complete FilterLabel for this filter.
        """
        return lsst.afw.image.FilterLabel(band=self.band, physical=self.physical_filter)
