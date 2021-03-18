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
    Tuple,
    TYPE_CHECKING
)

import lsst.afw.fits
from lsst.daf.butler import (
    DataCoordinate,
    DatasetRef,
    DatasetType,
    FileDataset,
    Progress,
)
from ..translators import Translator
from .parser import PathElementParser
from .scanner import PathElementHandler, DirectoryScanner

if TYPE_CHECKING:
    from lsst.daf.butler import FormatterParameter


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

    def __str__(self):
        return f"{type(self).__name__}({self._pattern}, isForFiles={self._isForFiles})"

    def isForFiles(self) -> bool:
        # Docstring inherited from PathElementHandler.
        return self._isForFiles

    @property
    def rank(self) -> int:
        # Docstring inherited from PathElementHandler.
        return 0

    def __call__(self, path: str, name: str,
                 datasets: Mapping[DatasetType, Mapping[Optional[str], List[FileDataset]]], *,
                 predicate: Callable[[DataCoordinate], bool]) -> bool:
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

    def __str__(self):
        return f"{type(self).__name__}(parser={self._parser})"

    def __call__(self, path: str, name: str,
                 datasets: Mapping[DatasetType, Mapping[Optional[str], List[FileDataset]]], *,
                 predicate: Callable[[DataCoordinate], bool]) -> bool:
        # Docstring inherited from PathElementParser.
        nextDataId2 = self._parser.parse(name, self.lastDataId2)
        if nextDataId2 is None:
            return False
        self.handle(path, nextDataId2, datasets, predicate=predicate)
        return True

    @property
    def rank(self) -> int:
        # Docstring inherited from PathElementParser.
        return len(self._parser.keys)

    @abstractmethod
    def handle(self, path: str, nextDataId2: dict,
               datasets: Mapping[DatasetType, Mapping[Optional[str], List[FileDataset]]], *,
               predicate: Callable[[DataCoordinate], bool]):
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
        predicate : `~collections.abc.Callable`
            A callable taking a single `DataCoordinate` argument and returning
            `bool`, indicating whether that (Gen3) data ID represents one
            that should be included in the scan.
        formatterMap : `dict`, optional
            Map dataset type to specialist formatter.
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

    def handle(self, path: str, nextDataId2: dict,
               datasets: Mapping[DatasetType, Mapping[Optional[str], List[FileDataset]]], *,
               predicate: Callable[[DataCoordinate], bool]):
        # Docstring inherited from ParsedPathElementHandler.
        if self._message is not None:
            self.log.warn("Skipping %s: %s", path, self._message)


class SubdirectoryHandler(ParsedPathElementHandler):
    """A `PathElementHandler` that uses a `DirectoryScanner` to recurse.

    Parameters
    ----------
    parser : `PathElementParser`
        An object that matches the path element this handler is responsible for
        and extracts a (partial) Gen2 data ID from it.
    progress : `Progress`, optional
        Object to use to report incremental progress.

    Notes
    -----
    The nested `DirectoryScanner` is default-constructed and should be
    populated with child handlers after the `SubdirectoryHandler` is created.
    """

    def __init__(self, parser: PathElementParser, progress: Optional[Progress] = None):
        super().__init__(parser=parser)
        self.scanner = DirectoryScanner(progress=progress)

    __slots__ = ("scanner",)

    def isForFiles(self) -> bool:
        # Docstring inherited from PathElementHandler.
        return False

    def handle(self, path: str, nextDataId2,
               datasets: Mapping[DatasetType, Mapping[Optional[str], List[FileDataset]]], *,
               predicate: Callable[[DataCoordinate], bool]):
        # Docstring inherited from ParsedPathElementHandler.
        if not nextDataId2:
            # We matched, and there's no data ID at all yet.  That means the
            # full path so far is just a fixed string so we should descend
            # and the match is exclusive.
            scan = True
        else:
            dataId3, _ = self.translate(nextDataId2, partial=True)
            if dataId3 is not None:
                scan = predicate(dataId3)
            else:
                scan = True
        if scan:
            for handler in self.scanner:
                handler.lastDataId2 = nextDataId2
            self.scanner.scan(path, datasets, predicate=predicate)

    def translate(self, dataId2: dict, *, partial: bool = False
                  ) -> Tuple[Optional[DataCoordinate], Optional[str]]:
        # Docstring inherited from PathElementHandler.
        for handler in self.scanner:
            # Since we're recursing, we're always asking for a partial match,
            # because the data ID we have corresponds to different level than
            # the one child handlers operate at.
            result, calibDate = handler.translate(dataId2, partial=True)
            if result is not None:
                return result, calibDate
        return None, None

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
    formatter : `lsst.daf.butler.Formatter` or `str`, optional
        A Gen 3 formatter class or fully-qualified name.
    """
    def __init__(self, parser: PathElementParser, translator: Translator, datasetType: DatasetType,
                 formatter: FormatterParameter = None):
        super().__init__(parser=parser)
        self._translator = translator
        self._datasetType = datasetType
        self._formatter = formatter

    __slots__ = ("_translator", "_datasetType", "_formatter")

    def __str__(self):
        return f"{type(self).__name__}({self._translator}, {self._datasetType})"

    def isForFiles(self) -> bool:
        # Docstring inherited from PathElementHandler.
        return True

    def handle(self, path: str, nextDataId2,
               datasets: Mapping[DatasetType, Mapping[Optional[str], List[FileDataset]]], *,
               predicate: Callable[[DataCoordinate], bool]):
        # Docstring inherited from ParsedPathElementHandler.
        dataId3, calibDate = self.translate(nextDataId2, partial=False)
        if predicate(dataId3):
            datasets[self._datasetType][calibDate].append(
                FileDataset(
                    refs=[DatasetRef(self._datasetType, dataId3)],
                    path=path, formatter=self._formatter
                )
            )

    def translate(self, dataId2: dict, *, partial: bool = False
                  ) -> Tuple[Optional[DataCoordinate], Optional[str]]:
        # Docstring inherited from PathElementHandler.
        rawDataId3, calibDate = self._translator(dataId2, partial=partial)
        if partial:
            return (
                DataCoordinate.standardize(rawDataId3, universe=self._datasetType.dimensions.universe),
                calibDate,
            )
        else:
            return (
                DataCoordinate.standardize(rawDataId3, graph=self._datasetType.dimensions),
                calibDate
            )


class MultiExtensionFileHandler(TargetFileHandler):
    """Handler for FITS files that store image and metadata in multiple HDUs
    per file, for example DECam raw and Community Pipeline calibrations.

    Notes
    -----
    For now, this is only used by DECam, and may need to be made more generic
    (e.g. making ``metadata['CCDNUM']`` use a configurable field) to be used
    with other obs packages.
    """
    def handle(self, path: str, nextDataId2,
               datasets: Mapping[DatasetType, Mapping[Optional[str], List[FileDataset]]], *,
               predicate: Callable[[DataCoordinate], bool]):
        dataId3, calibDate = self.translate(nextDataId2, partial=True)

        def get_detectors(filename):
            fitsData = lsst.afw.fits.Fits(filename, 'r')
            # NOTE: The primary header (HDU=0) does not contain detector data.
            detectors = []
            for i in range(1, fitsData.countHdus()):
                fitsData.setHdu(i)
                metadata = fitsData.readMetadata()
                detectors.append(metadata['CCDNUM'])
            return detectors

        if predicate(dataId3):
            detectors = get_detectors(path)
            refs = []
            for detector in detectors:
                newDataId3 = DataCoordinate.standardize(dataId3,
                                                        graph=self._datasetType.dimensions,
                                                        detector=detector)
                refs.append(DatasetRef(self._datasetType, newDataId3))

            datasets[self._datasetType][calibDate].append(
                FileDataset(refs=refs, path=path, formatter=self._formatter)
            )

    def translate(self, dataId2: dict, *, partial: bool = False
                  ) -> Tuple[Optional[DataCoordinate], Optional[str]]:
        assert partial is True, "We always require partial, to ignore 'ccdnum'"
        rawDataId3, calibDate = self._translator(dataId2, partial=partial)
        return (
            DataCoordinate.standardize(rawDataId3, universe=self._datasetType.dimensions.universe),
            calibDate,
        )
