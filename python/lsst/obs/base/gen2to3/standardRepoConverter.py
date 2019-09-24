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
from __future__ import annotations

__all__ = ["StandardRepoConverter"]

from dataclasses import dataclass
from typing import TYPE_CHECKING, Iterator, Tuple

from lsst.log import Log
from lsst.log.utils import temporaryLogLevel
from lsst.daf.persistence import Butler as Butler2
from lsst.daf.butler import DatasetType, DatasetRef, DataCoordinate
from .repoConverter import RepoConverter
from .filePathParser import FilePathParser
from .dataIdExtractor import DataIdExtractor

SKYMAP_DATASET_TYPES = {
    coaddName: f"{coaddName}Coadd_skyMap" for coaddName in ("deep", "goodSeeing", "dcr")
}

if TYPE_CHECKING:
    from lsst.skymap import BaseSkyMap
    from lsst.daf.butler import StorageClass
    from .cameraMapper import CameraMapper
    from ..mapping import Mapping as CameraMapperMapping  # disambiguate from collections.abc.Mapping


@dataclass
class FoundSkyMap:
    """Struct containing information about a SkyMap found in a Gen2 repository.
    """

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
        with temporaryLogLevel("CameraMapper", Log.ERROR):
            with temporaryLogLevel("HscMapper", Log.ERROR):
                self.butler2 = Butler2(self.root)
                self.mapper = self.butler2.getMapperClass(self.root)(root=self.root)
        self._foundSkyMapsByCoaddName = {}

    def isDatasetTypeSpecial(self, datasetTypeName: str) -> bool:
        # Docstring inherited from RepoConverter.
        return datasetTypeName in SKYMAP_DATASET_TYPES.values()

    def isDirectorySpecial(self, subdirectory: str) -> bool:
        # Docstring inherited from RepoConverter.
        return False

    def prep(self):
        # Docstring inherited from RepoConverter.
        self.task.log.info(f"Looking for skymaps in root {self.root}.")
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
            name = self.task.useSkyMap(instance)
            datasetType = DatasetType(datasetTypeName, dimensions=["skymap"],
                                      storageClass="SkyMap", universe=self.task.universe)
            dataId = DataCoordinate.standardize(skymap=name, universe=self.task.universe)
            struct = FoundSkyMap(name=name, instance=instance, coaddName=coaddName,
                                 ref=DatasetRef(datasetType, dataId),
                                 filename=self.butler2.getUri(datasetTypeName))
            self._foundSkyMapsByCoaddName[coaddName] = struct
            self.task.log.info("Found skymap %s in %s in %s.", name, datasetTypeName, self.root)
        super().prep()

    def iterMappings(self) -> Iterator[Tuple[str, CameraMapperMapping]]:
        # Docstring inherited from RepoConverter.
        for datasetTypeName, mapping in self.mapper.mappings.items():
            if datasetTypeName not in self.mapper.calibrations:
                yield datasetTypeName, mapping

    def makeDataIdExtractor(self, datasetTypeName: str, parser: FilePathParser,
                            storageClass: StorageClass) -> DataIdExtractor:
        # Docstring inherited from RepoConverter.
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
                        ("Dataset %s looks like it might need a skymap, but no %sCoadd_skyMap "
                         "found in repo %s."),
                        datasetTypeName, coaddName, self.root
                    )
        return DataIdExtractor(
            datasetTypeName,
            storageClass,
            filePathParser=parser,
            universe=self.task.universe,
            instrument=self.task.instrument.getName(),
            skyMap=struct.instance if struct is not None else None,
            skyMapName=struct.name if struct is not None else None,
        )

    def iterDatasets(self) -> Iterator[Tuple[str, DatasetRef]]:
        # Docstring inherited from RepoConverter.
        for struct in self._foundSkyMapsByCoaddName.values():
            if self.task.isDatasetTypeIncluded(struct.ref.datasetType.name):
                yield struct.filename, struct.ref
        yield from super().iterDatasets()

    # Class attributes that will be shadowed by public instance attributes;
    # defined here only for documentation purposes.

    butler2: Butler2
    """Gen2 butler associated with this repository.
    """

    mapper: CameraMapper
    """Gen2 mapper associated with this repository.
    """
