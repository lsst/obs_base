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

__all__ = ["FilePathParser"]

import re
from typing import Dict

from ..mapping import Mapping


class FilePathParser:
    """A callable object that extracts Gen2 data IDs from filenames
    corresponding to a particular Gen2 DatasetType.

    External code should use the `fromMapping` method to construct instances.

    Parameters
    ----------
    keys : `dict`
        Dictionary mapping Gen2 data ID key to the type of its associated
        value.
    regex : regular expression object
        Regular expression pattern with named groups for all data ID keys.
    """
    def __init__(self, keys: Dict[str, type], regex: re.Pattern):
        self.keys = keys
        self.regex = regex

    # Regular expression that matches a single substitution in
    # Gen2 CameraMapper template, such as "%(tract)04d".
    TEMPLATE_RE = re.compile(r"\%\((?P<name>\w+)\)[^\%]*?(?P<type>[idrs])")

    @classmethod
    def fromMapping(cls, mapping: Mapping) -> FilePathParser:
        """Construct a FilePathParser instance from a Gen2
        `lsst.obs.base.Mapping` instance.

        Parameters
        ----------
        mapping : `lsst.obs.base.Mapping`
            Mapping instance from a Gen2 `CameraMapper`.

        Returns
        -------
        parser : `FilePathParser`
            A new `FilePathParser` instance that can extract Gen2 data IDs
            from filenames with the given mapping's template.

        Raises
        ------
        RuntimeError
            Raised if the mapping has no template.
        """
        template = mapping.template
        keys = {}
        # The template string is something like
        # "deepCoadd/%(tract)04d-%(patch)s/%(filter)s"; each step of this
        # iterator corresponds to a %-tagged substitution string.
        # Our goal in all of this parsing is to turn the template into a regex
        # we can use to extract the associated values when matching strings
        # generated with the template.
        last = 0
        terms = []
        allKeys = mapping.keys()
        for match in cls.TEMPLATE_RE.finditer(template):
            # Copy the (escaped) regular string between the last substitution
            # and this one to the terms that will form the regex.
            terms.append(re.escape(template[last:match.start()]))
            # Pull out the data ID key from the name used in the
            # substitution string.  Use that and the substition
            # type to come up with the pattern to use in the regex.
            name = match.group("name")
            if name == "patch":
                pattern = r"\d+,\d+"
            elif match.group("type") in "id":  # integers
                pattern = r"0*\d+"
            else:
                pattern = ".+"
            # only use named groups for the first occurence of a key
            if name not in keys:
                terms.append(r"(?P<%s>%s)" % (name, pattern))
                keys[name] = allKeys[name]
            else:
                terms.append(r"(%s)" % pattern)
            # Remember the end of this match
            last = match.end()
        # Append anything remaining after the last substitution string
        # to the regex.
        terms.append(re.escape(template[last:]))
        return cls(keys, regex=re.compile("".join(terms)))

    def __call__(self, filePath: str) -> dict:
        """Extract a Gen2 data ID dictionary from the given path.

        Parameters
        ----------
        filePath : `str`
            Path and filename relative to the repository root.

        Returns
        -------
        dataId : `dict`
            Dictionary used to identify the dataset in the Gen2 butler, or
            None if the file was not recognized.
        """
        m = self.regex.fullmatch(filePath)
        if m is None:
            return None
        return {k: v(m.group(k)) for k, v in self.keys.items()}
