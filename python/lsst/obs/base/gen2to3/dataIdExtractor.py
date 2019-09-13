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

__all__ = ["DataIdExtractor"]

from typing import Union, Optional

from lsst.skymap import BaseSkyMap
from lsst.daf.butler import DataCoordinate, DatasetType, StorageClass, DimensionUniverse
from ..mapping import Mapping
from .filePathParser import FilePathParser
from .translators import Translator


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
    instrument : `str`
        Name of the Gen3 instrument for output data IDs that include that
        dimension.
    filePathParser : `lsst.daf.butler.gen2convert.FilePathParser`, optional
        Object responsible for reading a Gen2 data ID from a filename.  Will
        be created from ``mapper`` if not provided.
    translator : `lsst.daf.butler.gen2convert.Translator`, optional
        Object responsible for converting a Gen2 data ID into a Gen3 data ID.
        Will be created if not provided.
    mapping : `lsst.obs.base.Mapper`, optional
        Object that defines a Gen2 dataset type.  Must be provided if
        ``filePathParser`` is not.
    skyMap : `lsst.skymap.BaseSkyMap`, optional
        SkyMap that defines tracts and patches.  Must be provided for datasets
        with a ``patch`` key in their data IDs.
    skyMapName: `str`, optional
        Name of the Gen3 skymap for output data IDs that include that
        dimension.

    Raises
    ------
    RuntimeError
        Raised if the given mapping has no template.
    """

    def __init__(self, datasetTypeName: str, storageClass: Union[str, StorageClass], *,
                 universe: DimensionUniverse,
                 instrument: str,
                 filePathParser: Optional[FilePathParser] = None,
                 translator: Optional[Translator] = None,
                 mapping: Optional[Mapping] = None,
                 skyMap: Optional[BaseSkyMap] = None,
                 skyMapName: Optional[str] = None):
        if filePathParser is None:
            filePathParser = FilePathParser.fromMapping(mapping)
        self.filePathParser = filePathParser
        if translator is None:
            translator = Translator.makeMatching(datasetTypeName, filePathParser.keys,
                                                 instrument=instrument, skyMap=skyMap, skyMapName=skyMapName)
        self.translator = translator
        self.datasetType = DatasetType(datasetTypeName, dimensions=self.translator.dimensionNames,
                                       storageClass=storageClass, universe=universe)

    def apply(self, fileNameInRoot: str) -> Optional[DataCoordinate]:
        """Extract a Gen3 data ID from the given filename,

        Parameters
        ----------
        fileNameInRoot : `str`
            Filename relative to a Gen2 data repository root.

        Returns
        -------
        dataId : `lsst.daf.butler.DataCoordinate` or `None`
            The Gen3 data ID, or `None` if the file was not recognized as an
            instance of the extractor's dataset type.
        """
        gen2id = self.filePathParser(fileNameInRoot)
        if gen2id is None:
            return None
        return DataCoordinate.standardize(self.translator(gen2id), graph=self.datasetType.dimensions)
