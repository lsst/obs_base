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

__all__ = ["CalibRepoConverter"]

import os
import sqlite3
from datetime import datetime, timedelta
from typing import TYPE_CHECKING, Dict, Iterator, Tuple

from lsst.daf.butler import Butler as Butler3

from .repoConverter import RepoConverter
from .repoWalker import RepoWalker
from .translators import makeCalibrationLabel

if TYPE_CHECKING:
    from lsst.daf.butler import StorageClass
    from ..cameraMapper import CameraMapper
    from ..mapping import Mapping as CameraMapperMapping  # disambiguate from collections.abc.Mapping

CURATED_CALIBRATION_DATASET_TYPES = (
    "camera",
    "transmission_sensor",
    "transmission_filter",
    "transmission_optics",
    "transmission_atmosphere",
    "bfKernel"
)


class CalibRepoConverter(RepoConverter):
    """A specialization of `RepoConverter` for calibration repositories.

    Parameters
    ----------
    mapper : `CameraMapper`
        Gen2 mapper for the data repository.  The root associated with the
        mapper is ignored and need not match the root of the repository.
    kwds
        Additional keyword arguments are forwarded to (and required by)
        `RepoConverter`.
    """

    def __init__(self, *, mapper: CameraMapper, **kwds):
        super().__init__(**kwds)
        self.mapper = mapper
        self._datasetTypes = set()

    def isDatasetTypeSpecial(self, datasetTypeName: str) -> bool:
        # Docstring inherited from RepoConverter.
        return datasetTypeName in CURATED_CALIBRATION_DATASET_TYPES

    def iterMappings(self) -> Iterator[Tuple[str, CameraMapperMapping]]:
        # Docstring inherited from RepoConverter.
        yield from self.mapper.calibrations.items()

    def makeRepoWalkerTarget(self, datasetTypeName: str, template: str, keys: Dict[str, type],
                             storageClass: StorageClass) -> RepoWalker.Target:
        # Docstring inherited from RepoConverter.
        target = RepoWalker.Target(
            datasetTypeName=datasetTypeName,
            storageClass=storageClass,
            template=template,
            keys=keys,
            instrument=self.task.instrument.getName(),
            universe=self.task.registry.dimensions,
        )
        self._datasetTypes.add(target.datasetType)
        return target

    def insertDimensionData(self):
        # Docstring inherited from RepoConverter.
        # This has only been tested on HSC, and it's not clear how general it
        # is.  The catch is that it needs to generate calibration_label strings
        # consistent with those produced by the Translator system.
        db = sqlite3.connect(os.path.join(self.root, "calibRegistry.sqlite3"))
        db.row_factory = sqlite3.Row
        records = []
        for datasetType in self._datasetTypes:
            if "calibration_label" not in datasetType.dimensions:
                continue
            fields = ["validStart", "validEnd", "calibDate"]
            if "detector" in datasetType.dimensions.names:
                fields.append(self.task.config.ccdKey)
            else:
                fields.append(f"NULL AS {self.task.config.ccdKey}")
            if ("physical_filter" in datasetType.dimensions.names
                    or "abstract_filter" in datasetType.dimensions.names):
                fields.append("filter")
            else:
                fields.append("NULL AS filter")
            query = f"SELECT DISTINCT {', '.join(fields)} FROM {datasetType.name};"
            try:
                results = db.execute(query)
            except sqlite3.OperationalError as e:
                if (self.mapper.mappings[datasetType.name].tables is None
                        or len(self.mapper.mappings[datasetType.name].tables) == 0):
                    self.task.log.warn("Could not extract calibration ranges for %s in %s: %r",
                                       datasetType.name, self.root, e)
                    continue
                # Try using one of the alternate table names in the mapper (e.g. cpBias->bias for DECam).
                name = self.mapper.mappings[datasetType.name].tables[0]
                query = f"SELECT DISTINCT {', '.join(fields)} FROM {name};"
                try:
                    results = db.execute(query)
                except sqlite3.OperationalError as e:
                    self.task.log.warn("Could not extract calibration ranges for %s in %s: %r",
                                       datasetType.name, self.root, e)
                    continue
            for row in results:
                label = makeCalibrationLabel(datasetType.name, row["calibDate"],
                                             ccd=row[self.task.config.ccdKey], filter=row["filter"])
                records.append({
                    "instrument": self.task.instrument.getName(),
                    "name": label,
                    "datetime_begin": datetime.strptime(row["validStart"], "%Y-%m-%d"),
                    "datetime_end": datetime.strptime(row["validEnd"], "%Y-%m-%d") + timedelta(days=1),
                })
        if records:
            self.task.registry.insertDimensionData("calibration_label", *records)

    def ingest(self):
        # Docstring inherited from RepoConverter.
        if self.task.config.doWriteCuratedCalibrations:
            try:
                collections = self.getCollections(None)
            except LookupError as err:
                raise ValueError("Cannot ingest curated calibration into a calibration repo with no "
                                 "collections of its own; skipping.") from err
            # TODO: associate the curated calibrations with any other
            # collections and remove this assert; blocker is DM-23230.
            assert len(collections) == 1, \
                "Multiple collections for curated calibrations is not yet supported."
            butler3 = Butler3(butler=self.task.butler3, run=collections[0])
            self.task.instrument.writeCuratedCalibrations(butler3)
        super().ingest()

    # Class attributes that will be shadowed by public instance attributes;
    # defined here only for documentation purposes.

    mapper: CameraMapper
    """Gen2 mapper associated with this repository.
    """
