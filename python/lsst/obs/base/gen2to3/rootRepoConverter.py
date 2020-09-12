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

__all__ = ["RootRepoConverter"]

import os
import re
import itertools
from typing import TYPE_CHECKING, Iterator, Optional, Tuple, List, Set

from lsst.skymap import BaseSkyMap
from lsst.daf.butler import DatasetType, DatasetRef, DimensionGraph, FileDataset
from .standardRepoConverter import StandardRepoConverter

SKYMAP_DATASET_TYPES = {
    coaddName: f"{coaddName}Coadd_skyMap" for coaddName in ("deep", "goodSeeing", "dcr")
}

if TYPE_CHECKING:
    from lsst.daf.butler import SkyPixDimension


def getDataPaths(dataRefs):
    """Strip HDU identifiers from paths and return a unique set of paths.

    Parameters
    ----------
    dataRefs : `lsst.daf.persistence.ButlerDataRef`
        The gen2 datarefs to strip "[HDU]" values from.

    Returns
    -------
    paths : `set` [`str`]
        The unique file paths without appended "[HDU]".
    """
    paths = set()
    for dataRef in dataRefs:
        path = dataRef.getUri()
        # handle with FITS files with multiple HDUs (e.g. decam raw)
        paths.add(path.split('[')[0])
    return paths


class RootRepoConverter(StandardRepoConverter):
    """A specialization of `RepoConverter` for root data repositories.

    `RootRepoConverter` adds support for raw images (mostly delegated to the
    parent task's `RawIngestTask` subtask) and reference catalogs.

    Parameters
    ----------
    kwds
        Keyword arguments are forwarded to (and required by) `RepoConverter`.
    """

    def __init__(self, **kwds):
        super().__init__(run=None, **kwds)
        self._refCats: List[Tuple[str, SkyPixDimension]] = []
        if self.task.config.rootSkyMapName is not None:
            self._rootSkyMap = self.task.config.skyMaps[self.task.config.rootSkyMapName].skyMap.apply()
        else:
            self._rootSkyMap = None  # All access to _rootSkyMap is guarded
        self._chain = {}
        self._rawRefs = []

    def isDatasetTypeSpecial(self, datasetTypeName: str) -> bool:
        # Docstring inherited from RepoConverter.
        return (
            super().isDatasetTypeSpecial(datasetTypeName)
            or datasetTypeName in ("raw", "ref_cat", "ref_cat_config")
            # in Gen2, some of these are in the root repo, not a calib repo
            or datasetTypeName in self.instrument.getCuratedCalibrationNames()
        )

    def getSpecialDirectories(self) -> List[str]:
        # Docstring inherited from RepoConverter.
        return super().getSpecialDirectories() + ["CALIB", "ref_cats", "rerun"]

    def findMatchingSkyMap(self, datasetTypeName: str) -> Tuple[Optional[BaseSkyMap], Optional[str]]:
        # Docstring inherited from StandardRepoConverter.findMatchingSkyMap.
        skyMap, name = super().findMatchingSkyMap(datasetTypeName)
        if skyMap is None and self.task.config.rootSkyMapName is not None:
            self.task.log.debug(
                ("Assuming configured root skymap with name '%s' for dataset %s."),
                self.task.config.rootSkyMapName, datasetTypeName
            )
            skyMap = self._rootSkyMap
            name = self.task.config.rootSkyMapName
        return skyMap, name

    def runRawIngest(self):
        if self.task.raws is None:
            return
        self.task.log.info(f"Finding raws in root {self.root}.")
        if self.subset is not None:
            dataRefs = itertools.chain.from_iterable(
                self.butler2.subset(self.task.config.rawDatasetType,
                                    visit=visit) for visit in self.subset.visits
            )
        else:
            dataRefs = self.butler2.subset(self.task.config.rawDatasetType)
        dataPaths = getDataPaths(dataRefs)
        self.task.log.info("Ingesting raws from root %s into run %s.", self.root, self.task.raws.butler.run)
        self._rawRefs.extend(self.task.raws.run(dataPaths))
        self._chain = {self.task.raws.butler.run: {self.task.raws.datasetType.name}}

    def runDefineVisits(self):
        if self.task.defineVisits is None:
            return
        dimensions = DimensionGraph(self.task.universe, names=["exposure"])
        exposureDataIds = set(ref.dataId.subset(dimensions) for ref in self._rawRefs)
        self.task.log.info("Defining visits from exposures.")
        self.task.defineVisits.run(exposureDataIds)

    def prep(self):
        # Docstring inherited from RepoConverter.
        # Gather information about reference catalogs.
        if self.task.isDatasetTypeIncluded("ref_cat") and len(self.task.config.refCats) != 0:
            from lsst.meas.algorithms import DatasetConfig as RefCatDatasetConfig
            for refCat in os.listdir(os.path.join(self.root, "ref_cats")):
                path = os.path.join(self.root, "ref_cats", refCat)
                configFile = os.path.join(path, "config.py")
                if not os.path.exists(configFile):
                    continue
                if refCat not in self.task.config.refCats:
                    continue
                self.task.log.info(f"Preparing ref_cat {refCat} from root {self.root}.")
                onDiskConfig = RefCatDatasetConfig()
                onDiskConfig.load(configFile)
                if onDiskConfig.indexer.name != "HTM":
                    raise ValueError(f"Reference catalog '{refCat}' uses unsupported "
                                     f"pixelization '{onDiskConfig.indexer.name}'.")
                level = onDiskConfig.indexer["HTM"].depth
                try:
                    dimension = self.task.universe[f"htm{level}"]
                except KeyError as err:
                    raise ValueError(f"Reference catalog {refCat} uses HTM level {level}, but no htm{level} "
                                     f"skypix dimension is configured for this registry.") from err
                self.task.useSkyPix(dimension)
                self._refCats.append((refCat, dimension))
        if self.task.isDatasetTypeIncluded("brightObjectMask") and self.task.config.rootSkyMapName:
            self.task.useSkyMap(self._rootSkyMap, self.task.config.rootSkyMapName)
        super().prep()

    def iterDatasets(self) -> Iterator[FileDataset]:
        # Docstring inherited from RepoConverter.
        # Iterate over reference catalog files.
        for refCat, dimension in self._refCats:
            datasetType = DatasetType(refCat, dimensions=[dimension], universe=self.task.universe,
                                      storageClass="SimpleCatalog")
            if self.subset is None:
                regex = re.compile(r"(\d+)\.fits")
                for fileName in os.listdir(os.path.join(self.root, "ref_cats", refCat)):
                    m = regex.match(fileName)
                    if m is not None:
                        htmId = int(m.group(1))
                        dataId = self.task.registry.expandDataId({dimension: htmId})
                        yield FileDataset(path=os.path.join(self.root, "ref_cats", refCat, fileName),
                                          refs=DatasetRef(datasetType, dataId))
            else:
                for begin, end in self.subset.skypix[dimension]:
                    for htmId in range(begin, end):
                        dataId = self.task.registry.expandDataId({dimension: htmId})
                        yield FileDataset(path=os.path.join(self.root, "ref_cats", refCat, f"{htmId}.fits"),
                                          refs=DatasetRef(datasetType, dataId))
        yield from super().iterDatasets()

    def getRun(self, datasetTypeName: str, calibDate: Optional[str] = None) -> str:
        # Docstring inherited from RepoConverter.
        run = self.task.config.runs[datasetTypeName]
        self._chain.setdefault(run, set()).add(datasetTypeName)
        return run

    def getCollectionChain(self) -> List[Tuple[str, Set[str]]]:
        """Return tuples of run name and associated dataset type names that
        can be used to construct a chained collection that refers to the
        converted root repository (`list` [ `tuple` ]).
        """
        return list(self._chain.items())
