# This file is part of obs_base.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
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
"""Classes that transform (part of) a Gen2 filename template into a regular
expression that we can use to extract Gen2 data IDs from files.
"""
from __future__ import annotations

__all__ = ["PathElementParser"]


import logging
from abc import ABC, abstractmethod
import re
from typing import ClassVar, Dict, Optional


class FormattableRegEx(ABC):
    """An interface that generates a regular expression from a template and
    a data ID.

    This is used by `PathElementParser` to abstract over whether a path
    element's regex needs to include values from a data ID extracted from
    parent path elements or not.
    """

    @abstractmethod
    def format(self, dataId: dict) -> re.Pattern:
        """Substitute values from the given data ID and return a regular
        expression.

        Parameters
        ----------
        dataId : `dict`
            A dictionary whose entries may be used to format the regular
            expression.  May include unused entries.
        """
        raise NotImplementedError()


class FixedRegEx(FormattableRegEx):
    """A trivial implementation of `FormattableRegEx` that does no formatting.

    Parameters
    ----------
    regex : `re.Pattern`
        The fixed regular expression to return.
    """
    def __init__(self, regex: re.Pattern):
        self.regex = regex

    __slots__ = ("regex",)

    def format(self, dataId: dict) -> re.Pattern:
        # Docstring inherited from FormattableRegEx.
        return self.regex

    def __str__(self):
        return f"{type(self).__name__}({self.regex})"


class SubstitutableRegEx:
    """An implementation of `FormattableRegEx` formed from a concatenation of
    actual regular terms and %-style format strings.
    """
    def __init__(self):
        self._terms = []

    __slots__ = ("_terms",)

    def addRegexTerm(self, regex: str):
        """Add a regular expression term.
        """
        self._terms.append((regex, False))

    def addSubstitutionTerm(self, template: str):
        """Add a %-style format template term.
        """
        self._terms.append((template, True))

    def format(self, dataId: dict) -> re.Pattern:
        # Docstring inherited from FormattableRegEx.
        return re.compile("".join(re.escape(s % dataId) if isSub else s
                                  for s, isSub in self._terms))

    def simplify(self) -> FormattableRegEx:
        """Return a possibly-simplified version of this object.

        If `addSubstitionTerm` was never called, this returns a simple
        `FixedRegEx`.
        """
        if not any(isSub for _, isSub in self._terms):
            return FixedRegEx(re.compile("".join(s for s, _ in self._terms)))
        else:
            return self


class PathElementParser:
    """An object that matches Gen2 file names and extracts Gen2 data IDs.

    Parameters
    ----------
    target : `str`
        Either a full Gen2 path template or the part of one the corresponds to
        a single path element (a subdirectory or file name).
    allKeys : `dict` [`str`, `type`]
        A dictionary that provides types for all Gen2 data ID keys that are
        substituted into the given template.  Additional key-value pairs may
        be present and will be ignored.
    previousKeys : `dict` [`str`, `type`], optional
        A dictionary containing key strings and types for Gen2 data ID keys
        that have been extracted from previous path elements of the same
        template.  Values for these keys must be provided via the
        ``lastDataId`` argument when calling `parse`.
    """
    def __init__(self, template: str, allKeys: Dict[str, type], *,
                 previousKeys: Optional[Dict[str, type]] = None):
        self.template = template
        self.keys = {}
        # For each template path element, we iterate over each %-tagged
        # substitution string.
        last = 0
        self.regex = SubstitutableRegEx()
        for match in self.TEMPLATE_RE.finditer(self.template):
            # Copy the (escaped) regular string between the last substitution
            # and this one, escaping it appropriately.
            self.regex.addRegexTerm(re.escape(self.template[last:match.start()]))
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
            # Create a new named groups for the first occurence of a key
            # within an element.
            if name not in self.keys:
                if previousKeys and name in previousKeys:
                    # Key is new to this part of the template, but it appeared
                    # in some previous part of the template.  We'll format the
                    # original template with the data ID from that previous
                    # step later.
                    start, stop = match.span()
                    self.regex.addSubstitutionTerm(self.template[start:stop])
                else:
                    # Key is new; expect to extract a data ID value from it.
                    self.regex.addRegexTerm(r"(?P<%s>%s)" % (name, pattern))
                    self.keys[name] = allKeys[name]
            else:
                # Require a match with the last group for a second
                # occurrence.
                self.regex.addRegexTerm(r"(?P=<%s>)" % name)
            # Remember the end of this match
            last = match.end()
        # Append anything remaining after the last substitution string.
        self.regex.addRegexTerm(re.escape(self.template[last:]))
        # If there are no substitutions, join and compile into a single regex
        # now.
        self.regex = self.regex.simplify()

    __slots__ = ("keys", "template", "regex")

    TEMPLATE_RE: ClassVar[re.Pattern] = re.compile(r"\%\((?P<name>\w+)\)[^\%]*?(?P<type>[idrs])")
    """Regular expression that matches a single substitution in
    Gen2 CameraMapper template, such as "%(tract)04d".
    """

    def __str__(self):
        return f"{type(self).__name__}({self.regex})"

    def parse(self, name: str, lastDataId: dict, *, log: Optional[logging.Logger] = None) -> Optional[dict]:
        """Parse the path element.

        Parameters
        ----------
        name : `str`
            The path name to parse.
        lastDataId : `dict`
            The cumulative Gen2 data ID obtaining by calling `parse` on parsers
            for parent directories of the same path.
        log : `logging.Logger`, optional
            Log to use to report warnings and debug information.

        Returns
        -------
        dataId : `dict` or `None`
            Gen2 data ID that combines key-value pairs obtained from this path
            with those from ``lastDataId``.  `None` if ``name`` is not matched
            by this parser.  If the keys extracted are inconsistent with those
            in ``lastDataID``, a warning is sent to ``log`` and `None` is
            returned.
        """
        m = self.regex.format(lastDataId).fullmatch(name)
        if m is None:
            return None
        newDataId = {k: v(m.group(k)) for k, v in self.keys.items()}
        for commonKey in lastDataId.keys() & newDataId.keys():
            if newDataId[commonKey] != lastDataId[commonKey]:
                if log is not None:
                    log.warning("Inconsistent value %s=%r when parsing %r with %r.",
                                commonKey, newDataId[commonKey], name, lastDataId)
                return None
        newDataId.update(lastDataId)
        return newDataId

    keys: Dict[str, type]
    """Dictionary mapping Gen2 data ID key to the type of its associated
    value, covering only those keys that can be extracted from this path
    element.
    """

    template: str
    """The portion of the original Gen2 filename template that this parser was
    constructed with.
    """

    regex: re.Pattern
    """A regular expression that can be used to match the path element and
    populate the Gen2 data ID items whose keys are in ``keys``.
    """
