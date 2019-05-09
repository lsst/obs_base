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

__all__ = ("RepoConverter", "DataIdExtractor")

import os
import pickle
from collections import OrderedDict  # for move_to_end
from contextlib import contextmanager

import yaml

# register YAML loader for repositoryCfg.yaml files.
import lsst.daf.persistence.repositoryCfg   # noqa: F401

from lsst.daf.butler import DataId, DatasetType, DatasetRef
from lsst.daf.butler.gen2convert import FilePathParser, Translator
from lsst.log import Log
from lsst.utils import doImport


def findMapperClass(root):
    """Find the mapper class associated with a Gen2 data repository root.

    Parameters
    ----------
    root : `str`
        Path to a Gen2 repository root directory.

    Returns
    -------
    cls : `type`
        A subclass of `lsst.obs.base.CameraMapper`.

    Raises
    ------
    ValueError
        Raised if the directory does not appear to be the root of a
        Gen2 data repository.
    """
    cfgPath = os.path.join(root, "repositoryCfg.yaml")
    if os.path.exists(cfgPath):
        with open(cfgPath, "r") as f:
            repoCfg = yaml.load(f, Loader=yaml.UnsafeLoader)
            return repoCfg.mapper
    parentLinkPath = os.path.join(root, "_parent")
    if os.path.exists(parentLinkPath):
        return findMapperClass(os.readlink(parentLinkPath))
    mapperFilePath = os.path.join(root, "_mapper")
    if os.path.exists(mapperFilePath):
        with open(mapperFilePath, "r") as f:
            mapperClassPath = f.read().strip()
        return doImport(mapperClassPath)
    calibRegistryPath = os.path.join(root, "calibRegistry.sqlite3")
    if os.path.exists(calibRegistryPath):
        return findMapperClass(os.path.normpath(os.path.join(root, os.path.pardir)))
    raise ValueError(f"Could not determine (Gen2) mapper class for repo at '{root}'.")


@contextmanager
def loggedAt(name, level):
    """A context manager that temporarily sets the level of an `lsst.log.Log`.

    Parameters
    ----------
    name : `str`
        Name of the log to modify.
    level : `int`
        Integer enumeration constant indicating the temporary log level.
    """
    log = Log.getLogger(name)
    old = log.getLevel()
    log.setLevel(level)
    try:
        yield
    finally:
        log.setLevel(old)


class DataIdExtractor:
    """A class that extracts Gen3 data IDs from Gen2 filenames for a
    particular dataset type.

    Parameters
    ----------
    datasetTypeName : `str`
        Name of the dataset type the object will process.
    storageClass : `str` or `lsst.daf.butler.StorageClass`
        Gen3 storage class of the dataset type.
    universe : `lsst.daf.butler.DimensionUniverse`
        Object containing all dimension definitions.
    baseDataId : `dict`
        Key-value pairs that may need to appear in the Gen3 data ID, but can
        never be inferred from a Gen2 filename.  This should always include
        the instrument name (even Gen3 data IDs that don't involve the
        instrument dimension have instrument-dependent Gen2 filenames) and
        should also include the skymap name for any data ID that involves
        tracts or patches.
    filePathParser : `lsst.daf.butler.gen2convert.FilePathParser`, optional
        Object responsible for reading a Gen2 data ID from a filename.  Will
        be created from ``mapper`` if not provided.
    translator : `lsst.daf.butler.gen2convert.Translator`, optional
        Object responsible for converting a Gen2 data ID into a Gen3 data ID.
        Will be created if not provided.
    mapper : `lsst.obs.base.CameraMapper`, optional
        Object that defines Gen2 filename templates.  Must be provided if
        ``filePathParser`` is not.
    skyMap : `lsst.skymap.BaseSkyMap`, optional
        SkyMap that defines tracts and patches.  Must be provided for datasets
        with a ``patch`` key in their data IDs.
    """

    def __init__(self, datasetTypeName, storageClass, *, universe, baseDataId,
                 filePathParser=None, translator=None, mapper=None, skyMap=None):
        if filePathParser is None:
            filePathParser = FilePathParser.fromMapping(mapper.mappings[datasetTypeName])
        self.filePathParser = filePathParser
        if translator is None:
            translator = Translator.makeMatching(filePathParser.datasetType, baseDataId, skyMap=skyMap)
        self.translator = translator
        self.datasetType = DatasetType(datasetTypeName, dimensions=self.translator.dimensions,
                                       storageClass=storageClass)
        self.datasetType.normalize(universe=universe)

    def apply(self, fileNameInRoot):
        """Extract a Gen3 data ID from the given filename,

        Parameters
        ----------
        fileNameInRoot : `str`
            Filename relative to a Gen2 data repository root.

        Returns
        -------
        dataId : `lsst.daf.butler.DataId` or `None`
            The Gen3 data ID, or `None` if the file was not recognized as an
            instance of the extractor's dataset type.
        """
        gen2id = self.filePathParser(fileNameInRoot)
        if gen2id is None:
            return None
        return DataId(self.translator(gen2id), dimensions=self.datasetType.dimensions)


class RepoConverter:
    """A helper class that ingests (some of) the contents of a Gen2 data
    repository into a Gen3 data repository.

    Parameters
    ----------
    root : `str`
        Root of the Gen2 data repository.
    universe : `lsst.daf.butler.DimensionUniverse`
        Object containing all dimension definitions.
    baseDataId : `dict`
        Key-value pairs that may need to appear in the Gen3 data ID, but can
        never be inferred from a Gen2 filename.  This should always include
        the instrument name (even Gen3 data IDs that don't involve the
        instrument dimension have instrument-dependent Gen2 filenames) and
        should also include the skymap name in order to process any data IDs
        that involve tracts or patches.
    mapper : `lsst.obs.base.CameraMapper`, optional
        Object that defines Gen2 filename templates.  Will be identified,
        imported, and constructed from ``root`` if not provided.
    skyMap : `lsst.skymap.BaseSkyMap`, optional
        SkyMap that defines tracts and patches.  Must be provided in order to
        provess datasets with a ``patch`` key in their data IDs.
    """

    COADD_NAMES = ("deep", "goodSeeing", "dcr")
    REPO_ROOT_FILES = ("registry.sqlite3", "_mapper", "repositoryCfg.yaml",
                       "calibRegistry.sqlite3", "_parent")

    def __init__(self, root, *, universe, baseDataId, mapper=None, skyMap=None):
        self.root = root
        if mapper is None:
            # Shush spurious log messages from Gen2 Mapper classes.
            with loggedAt("CameraMapper", Log.ERROR):
                with loggedAt("HscMapper", Log.ERROR):
                    cls = findMapperClass(root)
                    mapper = cls(root=root)
        self.mapper = mapper
        self.universe = universe
        self.baseDataId = baseDataId
        self.extractors = OrderedDict()  # for move_to_end
        if "skymap" in baseDataId:
            if skyMap is None:
                for name in self.COADD_NAMES:
                    mapping = self.mapper.mappings.get(f"{name}Coadd_skyMap", None)
                    if mapping is None:
                        continue
                    filename = os.path.join(self.root, mapping.template)
                    if os.path.exists(filename):
                        if skyMap is not None:
                            raise ValueError("Multiple SkyMaps found in repository; please use multiple "
                                             "RepoConverters with an explicit skyMap argument for each.")
                        with open(filename, "rb") as f:
                            skyMap = pickle.load(f, encoding="latin1")
        self.skyMap = skyMap

    def addDatasetType(self, datasetTypeName, storageClass):
        """Add a dataset type to those recognized by the converter.

        Parameters
        ----------
        datasetTypeName : `str`
            String name of the dataset type.
        storageClass : `str` or `lsst.daf.butler.StorageClass`
            Gen3 storage class of the dataset type.

        Returns
        -------
        extractor : `DataIdExtractor`
            The object that will be used to extract data IDs for instances of
            this dataset type (also held internally, so the return value can
            usually be ignored).
        """
        r = DataIdExtractor(datasetTypeName, storageClass, mapper=self.mapper,
                            universe=self.universe, baseDataId=self.baseDataId, skyMap=self.skyMap)
        self.extractors[datasetTypeName] = r
        return r

    def extractDatasetRef(self, fileNameInRoot):
        """Extract a Gen3 `~lsst.daf.butler.DatasetRef` from a filename in a
        Gen2 data repository.

        Parameters
        ----------
        fileNameInRoot : `str`
            Name of the file, relative to the root of its Gen2 repository.

        Return
        ------
        ref : `lsst.daf.butler.DatasetRef` or `None`
            Reference to the Gen3 dataset that would be created by converting
            this file, or `None` if the file is not recognized as an instance
            of a dataset type known to this converter.
        """
        for datasetTypeName, extractor in self.extractors.items():
            dataId = extractor.apply(fileNameInRoot)
            if dataId is not None:
                # Move the extractor that matched to the front of the
                # dictionary, as we're likely to see instances of the
                # same DatasetType together.
                self.extractors.move_to_end(datasetTypeName, last=False)
                return DatasetRef(extractor.datasetType, dataId=dataId)
        return None

    def walkRepo(self, directory=None, skipDirs=()):
        """Recursively a (subset of) a Gen2 data repository, yielding files
        that may be convertible.

        Parameters
        ----------
        directory : `str`, optional
            A subdirectory of the repository root to process, instead of
            processing the entire repository.
        skipDirs : sequence of `str`
            Subdirectories that should be skipped.

        Yields
        ------
        fileNameInRoot : `str`
            Name of a file in the repository, relative to the root of the
            repository.
        """
        if directory is None:
            directory = self.root
        for dirPath, subdirNamesInDir, fileNamesInDir in os.walk(directory, followlinks=True):
            # Remove subdirectories that appear to be repositories themselves
            # from the walking
            def isRepoRoot(dirName):
                return any(os.path.exists(os.path.join(dirPath, dirName, f))
                           for f in self.REPO_ROOT_FILES)
            subdirNamesInDir[:] = [d for d in subdirNamesInDir if not isRepoRoot(d) and d not in skipDirs]
            # Loop over files in this directory, and ask per-DatasetType
            # extractors if they recognize them and can extract a data ID;
            # if so, ingest.
            dirPathInRoot = dirPath[len(self.root) + len(os.path.sep):]
            for fileNameInDir in fileNamesInDir:
                fileNameInRoot = os.path.join(dirPathInRoot, fileNameInDir)
                if fileNameInRoot in self.REPO_ROOT_FILES:
                    continue
                yield fileNameInRoot

    def convertRepo(self, butler, *, directory=None, transfer=None, formatter=None, skipDirs=()):
        """Ingest all recognized files into a Gen3 repository.

        Parameters
        ----------
        butler : `lsst.daf.butler.Butler`
            Gen3 butler that files should be ingested into.
        directory : `str`, optional
            A subdirectory of the repository root to process, instead of
            processing the entire repository.
        transfer : str, optional
            If not `None`, must be one of 'move', 'copy', 'hardlink', or
            'symlink' indicating how to transfer the file.
        formatter : `lsst.daf.butler.Formatter`, optional
            Formatter that should be used to retreive the Dataset.  If not
            provided, the formatter will be constructed according to
            Datastore configuration.  This should only be used when converting
            only a single dataset type multiple dataset types of the same
            storage class.
        skipDirs : sequence of `str`
            Subdirectories that should be skipped.
        """
        log = Log.getLogger("RepoConverter")
        for extractor in self.extractors.values():
            butler.registry.registerDatasetType(extractor.datasetType)
        skipped = {}
        for file in self.walkRepo(directory=directory, skipDirs=skipDirs):
            ref = self.extractDatasetRef(file)
            if ref is not None:
                try:
                    butler.ingest(os.path.join(self.root, file), ref, transfer=transfer, formatter=formatter)
                except Exception as err:
                    skipped.setdefault(type(err), []).append(str(err))
        if skipped:
            for cls, messages in skipped.items():
                log.warn("Skipped %s files due to exceptions of type %s.", len(messages), cls.__name__)
                if log.isDebugEnabled():
                    for message in messages:
                        log.debug(message)
