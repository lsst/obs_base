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
"""Interfaces and common code for recursively scanning directories for Gen2
dataset files.

The `PathElementHandler` ABC is defined here instead of ``handlers.py`` for
dependency reasons: `DirectoryScanner` uses the ABC, while its concrete
implementations use `DirectorySCanner`.
"""
from __future__ import annotations

__all__ = ["PathElementHandler", "DirectoryScanner"]

from abc import ABC, abstractmethod
import bisect
import logging
import os
from typing import (
    Callable,
    Iterator,
    List,
    Mapping,
    Optional,
    Tuple,
)

from lsst.daf.butler import (
    DataCoordinate,
    DatasetType,
    FileDataset,
    Progress,
)


class PathElementHandler(ABC):
    """An interface for objects that handle a single path element (directory or
    file) in a Gen2 data repository.

    Handlers are added to a `DirectoryScanner` instance, which then calls them
    until one succeeds when it processes each element in a directory.
    """
    def __init__(self):
        self.lastDataId2 = {}

    __slots__ = ("lastDataId2", "log")

    @abstractmethod
    def isForFiles(self) -> bool:
        """Report what kind of path element this object handlers.

        Returns
        -------
        Return `True` if this handler is for file entries, or `False` if it
        is for directories.
        """
        raise NotImplementedError()

    @abstractmethod
    def __call__(self, path: str, name: str,
                 datasets: Mapping[DatasetType, Mapping[Optional[str], List[FileDataset]]], *,
                 predicate: Callable[[DataCoordinate], bool]) -> bool:
        """Apply the handler to a file path.

        Parameters
        ----------
        path : `str`
            Full path of the file or directory.
        name : `str`
            Local name of the file or directory within its parent directory.
        datasets : `dict` [`DatasetType`, `dict` ]
            Dictionary that found datasets should be added to.  Nested dicts
            are keyed by either `None` (for most datasets) or a `str`
            "CALIBDATE" for calibration datasets.
        predicate : `~collections.abc.Callable`
            A callable taking a single `DataCoordinate` argument and returning
            `bool`, indicating whether that (Gen3) data ID represents one
            that should be included in the scan.'

        Returns
        -------
        matched : `bool`
            `True` if this handler was a match for the given path and no other
            handlers need to be tried on it, `False` otherwise.
        """
        raise NotImplementedError()

    @property
    @abstractmethod
    def rank(self) -> int:
        """Return a rough indication of how flexible this handler is in terms
        of the path element names it can match.

        Handlers that match a constant path element should always return zero.
        """
        raise NotImplementedError()

    def translate(self, dataId2: dict, *, partial: bool = False
                  ) -> Tuple[Optional[DataCoordinate], Optional[str]]:
        """Translate the given data ID from Gen2 to Gen3.

        The default implementation returns `None`.  Subclasses that are able
        to translate data IDs should override this method.

        Parameters
        ----------
        dataId2 : `dict`
            Gen2 data ID.
        partial : `bool`, optional
            If `True` (`False` is default) this is a partial data ID for some
            dataset, and missing keys are expected.

        Returns
        -------
        dataId3 : `lsst.daf.butler.DataCoordinate` or `None`
            A Gen3 data ID, or `None` if this handler cannot translate data
            IDs.
        calibDate : `str` or `None`
            A Gen2 calibration "CALIBDATE" value, or `None` if there was no
            such value in the template.
        """
        return None, None

    def __lt__(self, other: PathElementHandler):
        """Handlers are sorted by rank to reduce the possibility that more
        flexible handlers will have a chance to match something they shouldn't.
        """
        return self.rank < other.rank

    lastDataId2: dict
    """The Gen2 data ID obtained by processing parent levels in the directory
    tree.

    This attribute should be reset by calling code whenever a new parent
    directory is entered, before invoking `__call__`.
    """

    log: logging.Logger
    """A logger to use for all diagnostic messages (`logging.Logger`).

    This attribute is set on a handler in `DirectoryScanner.add`; this avoids
    needing to forward one through all subclass constructors.
    """


class DirectoryScanner:
    """An object that uses `PathElementHandler` instances to process the files
    and subdirectories in a directory tree.

    Parameters
    ----------
    log : `logging.Logger`, optional
        Log to use to report warnings and debug information.
    progress : `Progress`, optional
        Object to use to report incremental progress.
    """
    def __init__(self, log: Optional[logging.Logger] = None, progress: Optional[Progress] = None):
        self._files = []
        self._subdirectories = []
        if log is None:
            log = logging.getLogger("lsst.obs.base.gen2to3.repoWalker")
        self.log = log
        self.progress = progress

    __slots__ = ("_files", "_subdirectories", "log", "progress")

    def add(self, handler: PathElementHandler):
        """Add a new handler to the scanner.

        Parameters
        ----------
        handler : `PathElementHandler`
            The handler to be added.
        """
        handler.log = self.log
        if handler.isForFiles():
            bisect.insort(self._files, handler)
        else:
            bisect.insort(self._subdirectories, handler)

    def __iter__(self) -> Iterator[PathElementHandler]:
        """Iterate over all handlers.
        """
        yield from self._files
        yield from self._subdirectories

    def scan(self, path: str, datasets: Mapping[DatasetType, Mapping[Optional[str], List[FileDataset]]], *,
             predicate: Callable[[DataCoordinate], bool]):
        """Process a directory.

        Parameters
        ----------
        path : `str`
            Full path to the directory to be processed.
        datasets : `dict` [`DatasetType`, `list` ]
            Dictionary that found datasets should be added to.  Nested lists
            elements are tuples of `FileDataset` and an optional "CALIBDATE"
            `str` value (for calibration datasets only).
        predicate : `~collections.abc.Callable`
            A callable taking a single `DataCoordinate` argument and returning
            `bool`, indicating whether that (Gen3) data ID represents one
            that should be included in the scan.
        """
        with os.scandir(path) as iterator:
            unrecognized = []
            recognized = []
            for entry in iterator:
                if entry.is_file():
                    handlers = self._files
                elif entry.is_dir():
                    handlers = self._subdirectories
                else:
                    continue
                if self.progress is None:
                    # No progress reporting; look for a matching handler
                    # with an immediate depth-first search.
                    for handler in handlers:
                        if handler(entry.path, entry.name, datasets, predicate=predicate):
                            break
                    else:
                        unrecognized.append(entry.name)
                else:
                    # Caller wants progress reporting, but we won't know how
                    # many entries we'll have until we're done scanning.  So we
                    # save them in a list and process them in together later
                    # (essentially breadth-first search at this level).
                    recognized.append((entry.path, entry.name, handlers))
            if self.progress is not None:
                # Loop through the previously-recognized entries and process
                # them.
                for filepath, filename, handlers in self.progress.wrap(recognized, desc=f"Scanning {path}"):
                    for handler in handlers:
                        if handler(filepath, filename, datasets, predicate=predicate):
                            break
                    else:
                        unrecognized.append(entry.name)
            if unrecognized:
                self.log.warning("Skipped unrecognized entries in %s: %s", path, unrecognized)
