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
"""High-level interface to the Gen2 repository-walking functionality defined
by this package.
"""
from __future__ import annotations

__all__ = ["RepoWalker"]

from collections import defaultdict
import re
from typing import (
    Callable,
    ClassVar,
    Dict,
    Iterable,
    List,
    Mapping,
    Optional,
    Union,
)

from lsst.log import Log
from lsst.daf.butler import (
    DataCoordinate,
    DatasetType,
    FileDataset,
)
from .builders import BuilderTargetInput, BuilderSkipInput, BuilderTree
from .scanner import DirectoryScanner


class RepoWalker:
    """An object that recursively walks a Gen2 data repository tree, extracting
    Gen3 `FileDataset` objects and warning about unrecognized or unconvertable
    Gen2 datasets.

    Parameters
    ----------
    inputs : `~collections.abc.Iterable` of `Target` or `Skip`
        Structs that indicate dataset types to be extracted (`Target`) or
        explicitly skipped (`Skip`).  Skips may include a warning message to
        log when matching entries are encountered.
    fileIgnoreRegEx : `re.Pattern`, optional
        A regular expression pattern that identifies non-dataset files that
        can be ignored, to be applied at all levels of the directory tree.
    dirIgnoreRegEx : `re.Pattern`, optional
        A regular expression pattern that identifies non-dataset subdirectories
        that can be ignored, to be applied at all levels of the directory tree.
    log : `Log`, optional
            Logger for warnings and diagnostic information.
    """
    def __init__(self, inputs: Iterable[Union[Target, Skip]], *,
                 fileIgnoreRegEx: Optional[re.Pattern] = None,
                 dirIgnoreRegEx: Optional[re.Pattern] = None,
                 log: Optional[Log] = None):
        super().__init__()
        if log is None:
            log = Log.getLogger("obs.base.gen2to3.TranslatorFactory")
        self.log = log
        tree = BuilderTree()
        allKeys: Dict[str, type] = {}
        for leaf in inputs:
            tree.insert(0, leaf)
            for key, dtype in leaf.keys.items():
                if allKeys.setdefault(key, dtype) != dtype:
                    raise ValueError(f"Multiple types for key '{key}': {dtype} "
                                     f"(from {leaf.template}) vs. {allKeys[key]}.")
        tree, messages, pruned = tree.prune()
        if not pruned:
            self._scanner = DirectoryScanner(log=self.log)
            tree.fill(self._scanner, allKeys, {}, fileIgnoreRegEx=fileIgnoreRegEx,
                      dirIgnoreRegEx=dirIgnoreRegEx)
        else:
            # Nothing to do; just remember this for later to avoid disturbing
            # higher-level code with the fact that walk() will be a no-op.
            self._scanner = None

    Target: ClassVar[type] = BuilderTargetInput
    """An input struct type whose instances represent a dataset type to be
    extracted (`type`).
    """

    Skip: ClassVar[type] = BuilderSkipInput
    """An input struct type whose instances represent a dataset type to be
    explicitly skipped.
    """

    def walk(self, root: str, *, predicate: Optional[Callable[[DataCoordinate], bool]]
             ) -> Mapping[DatasetType, Mapping[Optional[str], List[FileDataset]]]:
        """Walk a Gen2 repository root to extract Gen3 `FileDataset` instances
        from it.

        Parameters
        ----------
        root : `str`
            Absolute path to the repository root.
        predicate : `~collections.abc.Callable`, optional
            If not `None`, a callable that returns `True` if a `DataCoordinate`
            is consistent with what we want to extract.  If ``predicate``
            returns `False`, the file or directory that data ID was extracted
            from will not be processed, even if it includes target dataset
            types.

        Returns
        -------
        datasets : `defaultdict` [`DatasetType`, `defaultdict` ]
            Extracted datasets, grouped by Gen3 `DatasetType`.  Nested dict
            keys are "CALIBDATE" strings (for calibration datasets) or `None`
            (otherwise).  Nested dict values are lists of `FileDataset`.
        """
        if predicate is None:
            def predicate(dataId: DataCoordinate) -> bool:
                return True
        datasets = defaultdict(lambda: defaultdict(list))
        if self._scanner is not None:
            self._scanner.scan(root, datasets, predicate=predicate)
        return datasets
