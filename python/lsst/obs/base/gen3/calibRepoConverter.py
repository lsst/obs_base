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

__all__ = ("CalibRepoConverter",)

import sqlite3
import os
from datetime import datetime

from lsst.daf.butler.gen2convert import makeCalibrationLabel
from .repoConverter import RepoConverter


class CalibRepoConverter(RepoConverter):
    """A helper class that ingests (some of) the contents of a Gen2 calibration
    data repository into a Gen3 data repository.

    This class must be used instead of the base `RepoConverter` to process
    dataset types that are associated with validity ranges.

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
        instrument dimension have instrument-dependent Gen2 filenames).
    mapper : `lsst.obs.base.CameraMapper`, optional
        Object that defines Gen2 filename templates.  Will be identified,
        imported, and constructed from ``root`` if not provided.
    """

    def __init__(self, root, *, universe, baseDataId, mapper=None):
        super().__init__(root, universe=universe, baseDataId=baseDataId, mapper=mapper)
        self.db = sqlite3.connect(os.path.join(root, "calibRegistry.sqlite3"))
        self.db.row_factory = sqlite3.Row
        self.dimensionEntries = []
        self.baseDataId = baseDataId

    def addDatasetType(self, datasetTypeName, storageClass):
        # Docstring inherited from RepoConverter.addDatasetType
        extractor = super().addDatasetType(datasetTypeName, storageClass)
        if "calibration_label" not in extractor.datasetType.dimensions:
            return extractor
        fields = ["validStart", "validEnd", "calibDate"]
        if "detector" in extractor.datasetType.dimensions.names:
            fields.append("ccd")
        else:
            fields.append("NULL AS ccd")
        if "physical_filter" in extractor.datasetType.dimensions.names:
            fields.append("filter")
        else:
            fields.append("NULL AS filter")
        query = f"SELECT DISTINCT {', '.join(fields)} FROM {datasetTypeName};"
        for row in self.db.execute(query):
            label = makeCalibrationLabel(datasetTypeName, row["calibDate"],
                                         ccd=row["ccd"], filter=row["filter"])
            self.dimensionEntries.append({
                "instrument": self.baseDataId["instrument"],
                "calibration_label": label,
                "valid_first": datetime.strptime(row["validStart"], "%Y-%m-%d"),
                "valid_last": datetime.strptime(row["validEnd"], "%Y-%m-%d"),
            })
        return extractor

    def convertRepo(self, butler, *, directory=None, transfer=None, formatter=None, skipDirs=()):
        # Docstring inherited from RepoConverter.convertRepo.
        butler.registry.addDimensionEntryList(butler.registry.dimensions["calibration_label"],
                                              self.dimensionEntries)
        super().convertRepo(butler, directory=directory, transfer=transfer, formatter=formatter,
                            skipDirs=tuple(skipDirs) + ("focalplane",))
