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
"""Classes used in `RepoWalker` construction.

The objects here form a temporary tree that is pruned and then transformed
into a similar tree of `PathElementHandler` instances.  See `BuilderNode`
method documentation for more information.
"""
from __future__ import annotations

__all__ = ["BuilderSkipInput", "BuilderTargetInput", "BuilderTree"]

from abc import ABC, abstractmethod
import os
import re
from typing import (
    Any,
    Dict,
    List,
    Optional,
    Tuple,
)

from lsst.daf.butler import DatasetType, DimensionUniverse, StorageClass
from ..translators import Translator
from .parser import PathElementParser
from .scanner import PathElementHandler, DirectoryScanner
from .handlers import IgnoreHandler, SubdirectoryHandler, SkipHandler, TargetFileHandler


class BuilderNode(ABC):
    """Abstract interface for nodes in the temporary tree that is used to
    construct a `RepoWalker`.
    """

    @abstractmethod
    def prune(self) -> Tuple[BuilderNode, List[str], bool]:
        """Attempt to prune this node and its children from the tree.

        Returns
        -------
        replacement : `BuilderNode`
            The result of recursively pruning child nodes; often just ``self``.
        messages : `list` [`str`]
            Warning messages that should be logged by a parent node when a
            matching path element is encountered, if this node is pruned.
        prune : `bool`
            If `True`, this node may be pruned from the tree (but will not
            necessarily be - it may correspond to a path element that should
            be skipped with siblings that should not be).
        """
        raise NotImplementedError()

    @abstractmethod
    def build(self, parser: PathElementParser, allKeys: Dict[str, type], cumulativeKeys: Dict[str, type], *,
              fileIgnoreRegEx: Optional[re.Pattern], dirIgnoreRegEx: Optional[re.Pattern]
              ) -> PathElementHandler:
        """Transform this node in the build tree into a corresponding
        `PathElementHandler`, recursing to any children.

        Must be called after `prune`.

        Parameters
        ----------
        parser : `PathElementParser`
            An object that matches the path element the new handler is
            responsible for and extracts a (partial) Gen2 data ID from it.
        allKeys : `dict` [`str`, `type`]
            A mapping from Gen2 data ID key to the type of its value.  Will
            contain all keys that may be extracted by the given parser, and
            possibly others.
        cumulativeKeys : `dict` [`str`, `type`], optional
            A dictionary containing key strings and types for Gen2 data ID keys
            that have been extracted from previous path elements for this
            template, including those extracted by ``parser``.

        Returns
        -------
        handler : `PathElementHandler`
            A new handler object.
        """
        raise NotImplementedError()


class BuilderInput(BuilderNode):
    """An intermediate base for `BuilderNode` classes that are provided as
    direct inputs to a `RepoWalker`, and generally correspond to exactly one
    Gen2 dataset type.

    Parameters
    ----------
    template : `str`
        The complete Gen2 template to be matched (not just the template for
        one path element).
    keys : `dict` [`str`, `type`]
        A mapping from Gen2 data ID key to the type of its value.
    """
    def __init__(self, template: str, keys: Dict[str, type]):
        self.template = template
        self.keys = keys
        self.elements = self.template.split(os.path.sep)

    template: str
    """The complete Gen2 template to be matched (`str`).
    """

    keys: Dict[str, type]
    """A mapping from Gen2 data ID key to the type of its value
    (`dict` [`str`, `type`]).
    """

    elements: List[str]
    """The path elements (file or directory levels) of `template`
    (`list` of `str`).
    """


class BuilderSkipInput(BuilderInput):
    """An input to a `RepoWalker` that indicates that matched files should be
    skipped, possibly with a warning message.

    BuilderSkipInputs can be pruned.  When they are not pruned, they build
    `SkipHandler` instances.

    Parameters
    ----------
    template : `str`
        The complete Gen2 template to be matched (not just the template for
        one path element).
    keys : `dict` [`str`, `type`]
        A mapping from Gen2 data ID key to the type of its value.
    message : `str`, optional
        If not `None`, a warning message that should be printed either when a
        matching file is enountered or a directory that may contain such files
        is skipped.
    isForFiles : `bool`, optional
        If `True` (default), this handler should be run on files.  Otherwise it
        should be run on directories.
    """
    def __init__(self, template: str, keys: Dict[str, type], message: Optional[str] = None, *,
                 isForFiles: bool = True):
        super().__init__(template=template, keys=keys)
        self._message = message
        self._isForFiles = isForFiles

    def build(self, parser: PathElementParser, allKeys: Dict[str, type], cumulativeKeys: Dict[str, type], *,
              fileIgnoreRegEx: Optional[re.Pattern], dirIgnoreRegEx: Optional[re.Pattern]
              ) -> PathElementHandler:
        # Docstring inherited from BuilderNode.
        return SkipHandler(parser=parser, isForFiles=self._isForFiles, message=self._message)

    def prune(self) -> Tuple[BuilderNode, List[str], bool]:
        # Docstring inherited from BuilderNode.
        return self, [self._message] if self._message is not None else [], True


class BuilderTargetInput(BuilderInput):
    """An input to a `RepoWalker` that matches files that correspond to
    datasets that we want to extract.

    BuilderTargetInputs can never be pruned, and always build
    `TargetFileHandler` instances.

    Parameters
    ----------
    datasetTypeName : `str`
        Name of the dataset type.
    template : `str`
        Full Gen2 filename template.
    keys : `dict` [`str`, `type`]
        Dictionary that maps Gen2 data ID key to the type of its value.
    storageClass : `StorageClass`
        `StorageClass` for the Gen3 dataset type.
    universe : `DimensionUniverse`
        All candidate dimensions for the Gen3 dataset type.
    kwargs:
        Additional keyword argumetns are passed to `Translator.makeMatching`,
        in along with ``datasetTypeName`` and ``keys``.
    """
    def __init__(self, *, datasetTypeName: str, template: str, keys: Dict[str, type],
                 storageClass: StorageClass, universe: DimensionUniverse, **kwargs: Any):
        super().__init__(template=template, keys=keys)
        self._translator = Translator.makeMatching(datasetTypeName, keys, **kwargs)
        self.datasetType = DatasetType(datasetTypeName, dimensions=self._translator.dimensionNames,
                                       storageClass=storageClass, universe=universe)

    def build(self, parser: PathElementParser, allKeys: Dict[str, type], cumulativeKeys: Dict[str, type], *,
              fileIgnoreRegEx: Optional[re.Pattern], dirIgnoreRegEx: Optional[re.Pattern]
              ) -> PathElementHandler:
        # Docstring inherited from BuilderNode.
        return TargetFileHandler(parser=parser, translator=self._translator, datasetType=self.datasetType)

    def prune(self) -> Tuple[BuilderNode, List[str], bool]:
        # Docstring inherited from BuilderNode.
        return self, [], False

    datasetType: DatasetType
    """The Gen3 dataset type extracted by the hander this object builds
    (`lsst.daf.butler.DatasetType`).
    """


class BuilderPrunedTree(BuilderNode):
    """A `BuilderNode` that represents a subdirectory to be skipped,
    created by pruning `BuilderTree` that contained only `BuilderSkipInput`
    instances.

    BuilderPrunedTrees can be pruned.  When they are not pruned, they
    build `SkipHandler` instances.

    Parameters
    ----------
    messages : `list` [`str`]
        A list of warning messages to be printed when the handler produced by
        this builder matches a subdirectory.
    """

    def __init__(self, messages: List[str]):
        self._messages = messages

    def build(self, parser: PathElementParser, allKeys: Dict[str, type], cumulativeKeys: Dict[str, type], *,
              fileIgnoreRegEx: Optional[re.Pattern], dirIgnoreRegEx: Optional[re.Pattern]
              ) -> PathElementHandler:
        # Docstring inherited from BuilderNode.
        message = ": ".join(self._messages) if self._messages else None
        return SkipHandler(parser=parser, isForFiles=False, message=message)

    def prune(self) -> Tuple[BuilderNode, List[str], bool]:
        # Docstring inherited from BuilderNode.
        return self, self._messages, True


class BuilderTree(BuilderNode):
    """A `BuilderNode` that represents a directory.

    This is the only `BuilderNode` class that is not a leaf node.  If all
    of its children can be pruned, it is replaced by a `BuilderPrunedTree`
    (which can then be pruned itself).  It builds `SubdirectoryHandler`
    instances when not pruned.
    """
    def __init__(self):
        self._children = {}  # Maps template path element to BuilderNode

    def insert(self, level: int, leaf: BuilderInput):
        """Insert an input leaf node into the tree, recursively constructing
        intermediate parents in order to put it at the right level.

        Parameters
        ----------
        level : `int`
            The level ``self``is at in the larger tree, with zero the
            repository root.  The right level for the leaf is given by the
            length of ``leaf.elements``.
        leaf : `BuilderInput`
            The leaf node to insert.
        """
        nextLevel = level + 1
        if nextLevel == len(leaf.elements):
            self._children[leaf.elements[level]] = leaf
        else:
            child = self._children.setdefault(leaf.elements[level], BuilderTree())
            child.insert(nextLevel, leaf)

    def fill(self, scanner: DirectoryScanner, allKeys: Dict[str, type], previousKeys: Dict[str, type], *,
             fileIgnoreRegEx: Optional[re.Pattern], dirIgnoreRegEx: Optional[re.Pattern]):
        """Fill a `DirectoryScanner` instance by recursively building all
        child nodes.

        Parameters
        ----------
        scanner : `DirectoryScanner`
            Object to populate.
        allKeys : `dict` [`str`, `type`]
            Mapping from Gen2 data ID key to its value type, covering all keys
            that could be used in any child template.
        previousKeys : `dict` [`str`, `type`], optional
            A dictionary containing key strings and types for Gen2 data ID keys
            that have been extracted from previous path elements of the same
            template.
        """
        if fileIgnoreRegEx is not None:
            scanner.add(IgnoreHandler(fileIgnoreRegEx, isForFiles=True))
        if dirIgnoreRegEx is not None:
            scanner.add(IgnoreHandler(dirIgnoreRegEx, isForFiles=False))
        for template, child in self._children.items():
            parser = PathElementParser(template, allKeys, previousKeys=previousKeys)
            cumulativeKeys = previousKeys.copy()
            cumulativeKeys.update(parser.keys)
            scanner.add(child.build(parser, allKeys, cumulativeKeys, fileIgnoreRegEx=fileIgnoreRegEx,
                                    dirIgnoreRegEx=dirIgnoreRegEx))

    def prune(self) -> Tuple[BuilderNode, List[str], bool]:
        # Docstring inherited from BuilderNode.
        toPruneThis = True
        newChildren = {}
        messages = []
        # Recursively prune children.
        for template, child in list(self._children.items()):
            newChild, childMessages, toPruneChild = child.prune()
            newChildren[template] = newChild
            messages.extend(childMessages)
            if not toPruneChild:
                toPruneThis = False
        self._children = newChildren
        if toPruneThis:
            return BuilderPrunedTree(messages), messages, True
        else:
            return self, [], False

    def build(self, parser: PathElementParser, allKeys: Dict[str, type], cumulativeKeys: Dict[str, type], *,
              fileIgnoreRegEx: Optional[re.Pattern], dirIgnoreRegEx: Optional[re.Pattern]
              ) -> PathElementHandler:
        # Docstring inherited from BuilderNode.
        built = SubdirectoryHandler(parser)
        self.fill(built.scanner, allKeys, cumulativeKeys, fileIgnoreRegEx=fileIgnoreRegEx,
                  dirIgnoreRegEx=dirIgnoreRegEx)
        return built
