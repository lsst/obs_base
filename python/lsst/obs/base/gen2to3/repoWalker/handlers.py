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
"""Concrete implementations of `PathElementHandler`.

The `PathElementHandler` ABC is defined in ``scanner.py`` instead of here to
avoid a circular dependency between modules.
"""
from __future__ import annotations

__all__ = ["IgnoreHandler", "SkipHandler", "SubdirectoryHandler", "TargetFileHandler"]

from abc import abstractmethod
import re
from typing import (
    Callable,
    List,
    Mapping,
    Optional,
)

from lsst.log import Log
from lsst.daf.butler import (
    DataCoordinate,
    DatasetRef,
    DatasetType,
    FileDataset,
)
from ..translators import Translator
from .parser import PathElementParser
from .scanner import PathElementHandler, DirectoryScanner


class IgnoreHandler(PathElementHandler):
    """A `PathElementHandler` that matches via a regular expression, and does
    nothing.

    An `IgnoreHandler` is used to ignore file or directory patterns that can
    occur at any level in the directory tree, and have no relation to any
    Gen2 filename template.

    Parameters
    ----------
    pattern : `re.Pattern`
        A regular expression pattern.
    isForFiles : `bool`
        Whether this handler should be applied to files (`True`) or
        directories (`False`).
    """
    def __init__(self, pattern: re.Pattern, isForFiles: bool):
        super().__init__()
        self._pattern = pattern
        self._isForFiles = isForFiles

    __slots__ = ("_pattern", "_isForFiles")

    def isForFiles(self) -> bool:
        # Docstring inherited from PathElementHandler.
        return self._isForFiles

    @property
    def rank(self) -> int:
        # Docstring inherited from PathElementHandler.
        return 0

    def __call__(self, path: str, name: str, datasets: Mapping[DatasetType, List[FileDataset]], *,
                 log: Log, predicate: Callable[[DataCoordinate], bool]) -> bool:
        # Docstring inherited from PathElementHandler.
        if self._pattern.fullmatch(name):
            return True
        else:
            return False


class ParsedPathElementHandler(PathElementHandler):
    """An intermediate base class for `PathElementHandler` classes that utilize
    a `PathElementParser` to match a Gen2 filename template.

    Parameters
    ----------
    parser : `PathElementParser`
        An object that matches the path element this handler is responsible for
        and extracts a (partial) Gen2 data ID from it.
    """
    def __init__(self, parser: PathElementParser):
        super().__init__()
        self._parser = parser

    __slots__ = ("_parser",)

    def __call__(self, path: str, name: str, datasets: Mapping[DatasetType, List[FileDataset]], *,
                 log: Log, predicate: Callable[[DataCoordinate], bool]) -> bool:
        # Docstring inherited from PathElementParser.
        nextDataId2 = self._parser.parse(name, self.lastDataId2, log=log)
        if nextDataId2 is None:
            return False
        self.handle(path, nextDataId2, datasets, log=log, predicate=predicate)
        return True

    @property
    def rank(self) -> int:
        # Docstring inherited from PathElementParser.
        return len(self._parser.keys)

    @abstractmethod
    def handle(self, path: str, nextDataId2: dict, datasets: Mapping[DatasetType, List[FileDataset]], *,
               log: Log, predicate: Callable[[DataCoordinate], bool]):
        """Customization hook for ``__call__``.

        Subclasses must override this method, while external callers (i.e.
        `DirectoryScanner` should instead invoke `__call__`.

        Parameters
        ----------
        path : `str`
            Full path of the file or directory.
        nextDataId2 : `dict`
            Gen2 data ID (usually partial) extracted from the path so far.
        datasets : `dict` [`DatasetType`, `list` [`FileDataset`] ]
            Dictionary that found datasets should be added to.
        log : `Log`, optional
            Log to use to report warnings and debug information.
        predicate : `~collections.abc.Callable`
            A callable taking a single `DataCoordinate` argument and returning
            `bool`, indicating whether that (Gen3) data ID represents one
            that should be included in the scan.
        """
        raise NotImplementedError()


class SkipHandler(ParsedPathElementHandler):
    """A `ParsedPathElementHandler` that does nothing with an entry other
    optionally logging a warning message.

    A `SkipHandler` is used for Gen2 datasets that we can recognize but do not
    want to (or cannot) extract Gen3 datasets from, or other files/directories
    that alway appears at a fixed level in the diectory tree.

    Parameters
    ----------
    parser : `PathElementParser`
        An object that matches the path element this handler is responsible for
        and extracts a (partial) Gen2 data ID from it.
    isForFiles : `bool`
        Whether this handler should be applied to files (`True`) or
        directories (`False`).
    message : `str`, optional
        A message to log at warning level when this handler matches a path
        entry.  If `None`, matched entrie will be silently skipped.
    """
    def __init__(self, parser: PathElementParser, isForFiles: bool, message: Optional[str]):
        super().__init__(parser=parser)
        self._isForFiles = isForFiles
        self._message = message

    __slots__ = ("_message", "_isForFiles")

    def isForFiles(self) -> bool:
        # Docstring inherited from PathElementHandler.
        return self._isForFiles

    def handle(self, path: str, nextDataId2: dict, datasets: Mapping[DatasetType, List[FileDataset]], *,
               log: Log, predicate: Callable[[DataCoordinate], bool]):
        # Docstring inherited from ParsedPathElementHandler.
        if self._message is not None:
            log.warn("Skipping %s: %s", path, self._message)


class SubdirectoryHandler(ParsedPathElementHandler):
    """A `PathElementHandler` that uses a `DirectoryScanner` to recurse.

    Parameters
    ----------
    parser : `PathElementParser`
        An object that matches the path element this handler is responsible for
        and extracts a (partial) Gen2 data ID from it.

    Notes
    -----
    The nested `DirectoryScanner` is default-constructed and should be
    populated with child handlers after the `SubdirectoryHandler` is created.
    """

    def __init__(self, parser: PathElementParser):
        super().__init__(parser=parser)
        self.scanner = DirectoryScanner()

    __slots__ = ("scanner",)

    def isForFiles(self) -> bool:
        # Docstring inherited from PathElementHandler.
        return False

    def handle(self, path: str, nextDataId2, datasets: Mapping[DatasetType, List[FileDataset]], *,
               log: Log, predicate: Callable[[DataCoordinate], bool]):
        # Docstring inherited from ParsedPathElementHandler.
        if not nextDataId2:
            # We matched, and there's no data ID at all yet.  That means the
            # full path so far is just a fixed string so we should descend
            # and the match is exclusive.
            scan = True
        else:
            dataId3 = self.translate(nextDataId2, partial=True, log=log)
            if dataId3 is not None:
                scan = predicate(dataId3)
            else:
                scan = True
        if scan:
            for handler in self.scanner:
                handler.lastDataId2 = nextDataId2
            self.scanner.scan(path, datasets, log=log, predicate=predicate)

    def translate(self, dataId2: dict, *, partial: bool = False, log: Log) -> Optional[DataCoordinate]:
        # Docstring inherited from PathElementHandler.
        for handler in self.scanner:
            # Since we're recursing, we're always asking for a partial match,
            # because the data ID we have corresponds to different level than
            # the one child handlers operate at.
            result = handler.translate(dataId2, partial=True, log=log)
            if result is not None:
                return result
        return None

    scanner: DirectoryScanner
    """Scanner object that holds handlers for the entries of the subdirectory
    matched by this handler (`DirectoryScanner`).
    """


class TargetFileHandler(ParsedPathElementHandler):
    """A `PathElementHandler` that matches files that correspond to target
    datasets and outputs `FileDataset` instances for them.

    Parameters
    ----------
    parser : `PathElementParser`
        An object that matches the path element this handler is responsible for
        and extracts a (partial) Gen2 data ID from it.
    translator : `Translator`
        Object that translates data IDs from Gen2 to Gen3.
    datasetType : `lsst.daf.butler.DatasetType`
        Gen3 dataset type for the datasets this handler matches.
    """
    def __init__(self, parser: PathElementParser, translator: Translator, datasetType: DatasetType):
        super().__init__(parser=parser)
        self._translator = translator
        self._datasetType = datasetType

    __slots__ = ("_translator", "_datasetType")

    def isForFiles(self) -> bool:
        # Docstring inherited from PathElementHandler.
        return True

    def handle(self, path: str, nextDataId2, datasets: Mapping[DatasetType, List[FileDataset]], *,
               log: Log, predicate: Callable[[DataCoordinate], bool]):
        # Docstring inherited from ParsedPathElementHandler.
        dataId3 = self.translate(nextDataId2, partial=False, log=log)
        if predicate(dataId3):
            datasets[self._datasetType].append(FileDataset(refs=[DatasetRef(self._datasetType, dataId3)],
                                                           path=path))

    def translate(self, dataId2: dict, *, partial: bool = False, log: Log) -> Optional[DataCoordinate]:
        # Docstring inherited from PathElementHandler.
        rawDataId3 = self._translator(dataId2, partial=partial, log=log)
        if partial:
            return DataCoordinate.standardize(rawDataId3, universe=self._datasetType.dimensions.universe)
        else:
            return DataCoordinate.standardize(rawDataId3, graph=self._datasetType.dimensions)
