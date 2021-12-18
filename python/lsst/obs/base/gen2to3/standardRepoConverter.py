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

__all__ = ["StandardRepoConverter"]

from dataclasses import dataclass
from typing import TYPE_CHECKING, Dict, Iterator, List, Mapping, Optional, Tuple

from lsst.daf.butler import DataCoordinate, DatasetRef, DatasetType, FileDataset
from lsst.daf.persistence import Butler as Butler2
from lsst.log import Log
from lsst.log.utils import temporaryLogLevel
from lsst.skymap import BaseSkyMap

from .repoConverter import RepoConverter
from .repoWalker import RepoWalker

SKYMAP_DATASET_TYPES = {coaddName: f"{coaddName}Coadd_skyMap" for coaddName in ("deep", "goodSeeing", "dcr")}

if TYPE_CHECKING:
    from lsst.daf.butler import FormatterParameter, StorageClass

    from ..mapping import Mapping as CameraMapperMapping  # disambiguate from collections.abc.Mapping
    from .cameraMapper import CameraMapper
    from .repoWalker.scanner import PathElementHandler


@dataclass
class FoundSkyMap:
    """Struct containing information about a SkyMap in a Gen2 repository."""

    name: str
    """Name of the skymap used in Gen3 data IDs.
    """

    instance: BaseSkyMap
    """An instance of the actual skymap class.
    """

    coaddName: str
    """The coadd name used as a prefix for the dataset type this skymap was
    found in.
    """

    ref: DatasetRef
    """A `DatasetRef` that can be used to ingest the skymap dataset into a
    Gen3 repository.
    """

    filename: str
    """Name of the file containing the skymap dataset, relative to the
    repository root.
    """


class StandardRepoConverter(RepoConverter):
    """A specialization of `RepoConverter` for non-calibration repositories.

    Parameters
    ----------
    kwds
        Keyword arguments are forwarded to (and required by) `RepoConverter`.
    """

    def __init__(self, **kwds):
        super().__init__(**kwds)
        # Shush noisy log messages from Gen2 Mapper classes.
        # These are currently lsst.log loggers.
        with temporaryLogLevel("CameraMapper", Log.ERROR):
            with temporaryLogLevel("HscMapper", Log.ERROR):
                self.butler2 = Butler2(self.root)
                self.mapper = self.butler2.getMapperClass(self.root)(root=self.root)
        self._foundSkyMapsByCoaddName = {}
        self._chain = []

    def isDatasetTypeSpecial(self, datasetTypeName: str) -> bool:
        # Docstring inherited from RepoConverter.
        return datasetTypeName in SKYMAP_DATASET_TYPES.values()

    def prep(self):
        # Docstring inherited from RepoConverter.
        self.task.log.info("Looking for skymaps in root %s", self.root)
        for coaddName, datasetTypeName in SKYMAP_DATASET_TYPES.items():
            if not self.task.isDatasetTypeIncluded(datasetTypeName):
                continue
            try:
                exists = self.butler2.datasetExists(datasetTypeName)
            except AttributeError:
                # This mapper doesn't even define this dataset type.
                continue
            if not exists:
                continue
            instance = self.butler2.get(datasetTypeName)
            name = self.task.useSkyMap(instance, datasetTypeName)
            datasetType = DatasetType(
                datasetTypeName, dimensions=["skymap"], storageClass="SkyMap", universe=self.task.universe
            )
            dataId = DataCoordinate.standardize(skymap=name, universe=self.task.universe)
            struct = FoundSkyMap(
                name=name,
                instance=instance,
                coaddName=coaddName,
                ref=DatasetRef(datasetType, dataId),
                filename=self.butler2.getUri(datasetTypeName),
            )
            self._foundSkyMapsByCoaddName[coaddName] = struct
            self.task.log.info("Found skymap %s in %s in %s.", name, datasetTypeName, self.root)
        super().prep()

    def iterMappings(self) -> Iterator[Tuple[str, CameraMapperMapping]]:
        # Docstring inherited from RepoConverter.
        for datasetTypeName, mapping in self.mapper.mappings.items():
            if datasetTypeName not in self.mapper.calibrations:
                yield datasetTypeName, mapping

    def findMatchingSkyMap(self, datasetTypeName: str) -> Tuple[Optional[BaseSkyMap], Optional[str]]:
        """Return the appropriate SkyMap for the given dataset type.

        Parameters
        ----------
        datasetTypeName : `str`
            Name of the dataset type for which a skymap is sought.

        Returns
        -------
        skyMap : `BaseSkyMap` or `None`
            The `BaseSkyMap` instance, or `None` if there was no match.
        skyMapName : `str` or `None`
            The Gen3 name for the SkyMap, or `None` if there was no match.
        """
        # Use deepCoadd_skyMap by default; there are some dataset types
        # that use it but don't have "deep" anywhere in their name.
        struct = self._foundSkyMapsByCoaddName.get("deep")
        for coaddName in SKYMAP_DATASET_TYPES.keys():
            if coaddName in datasetTypeName:
                try:
                    struct = self._foundSkyMapsByCoaddName[coaddName]
                    break
                except KeyError:
                    # Don't use the default, since we did find a specific
                    # coaddName.
                    struct = None
                    self.task.log.debug(
                        "Dataset %s looks like it might need a skymap, but no %sCoadd_skyMap "
                        "found in repo %s.",
                        datasetTypeName,
                        coaddName,
                        self.root,
                    )
        if struct is not None:
            return struct.instance, struct.name
        else:
            return None, None

    def makeRepoWalkerTarget(
        self,
        datasetTypeName: str,
        template: str,
        keys: Dict[str, type],
        storageClass: StorageClass,
        formatter: FormatterParameter = None,
        targetHandler: Optional[PathElementHandler] = None,
    ) -> RepoWalker.Target:
        # Docstring inherited from RepoConverter.
        skyMap, skyMapName = self.findMatchingSkyMap(datasetTypeName)
        return RepoWalker.Target(
            datasetTypeName=datasetTypeName,
            storageClass=storageClass,
            template=template,
            keys=keys,
            universe=self.task.registry.dimensions,
            instrument=self.task.instrument.getName(),
            skyMap=skyMap,
            skyMapName=skyMapName,
            formatter=formatter,
            targetHandler=targetHandler,
            translatorFactory=self.task.translatorFactory,
        )

    def iterDatasets(self) -> Iterator[FileDataset]:
        # Docstring inherited from RepoConverter.
        yield from super().iterDatasets()

    def getRun(self, datasetTypeName: str, calibDate: Optional[str] = None) -> str:
        # Docstring inherited from RepoConverter.
        run = self.task.config.runsForced.get(datasetTypeName)
        if run is None:
            if self._run is not None:
                run = self._run
            else:
                run = self.task.config.runs.get(datasetTypeName)
        if run is None:
            raise ValueError(
                f"No default run for repo at {self.root}, and no override for dataset {datasetTypeName}."
            )
        if run not in self._chain:
            self._chain.append(run)
        return run

    def getCollectionChain(self) -> List[str]:
        """Return run names that can be used to construct a chained collection
        that refers to the converted repository (`list` [ `str` ]).
        """
        return self._chain

    def _finish(
        self, datasets: Mapping[DatasetType, Mapping[Optional[str], List[FileDataset]]], count: int
    ) -> None:
        # Docstring inherited from RepoConverter.
        super()._finish(datasets, count)
        if self._foundSkyMapsByCoaddName:
            self._chain.append(BaseSkyMap.SKYMAP_RUN_COLLECTION_NAME)

    # Class attributes that will be shadowed by public instance attributes;
    # defined here only for documentation purposes.

    butler2: Butler2
    """Gen2 butler associated with this repository.
    """

    mapper: CameraMapper
    """Gen2 mapper associated with this repository.
    """
