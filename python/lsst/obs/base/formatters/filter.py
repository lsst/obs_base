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

__all__ = ("FilterFormatter",)

import yaml
from lsst.afw.image import Filter

from typing import (
    Any,
    Optional,
    Type,
)

from lsst.daf.butler.formatters.file import FileFormatter


class FilterFormatter(FileFormatter):
    """Read and write `~lsst.afw.image.Filter` filter information."""

    extension = ".yaml"

    unsupportedParameters = None
    """This formatter does not support any parameters."""

    def _readFile(self, path: str, pytype: Type[Any] = None) -> Any:
        """Read a file from the path in YAML format.

        Parameters
        ----------
        path : `str`
            Path to use to open the file.
        pytype : `class`, optional
            The type expected to be returned.

        Returns
        -------
        data : `object`
            Either data as Python object read from YAML file, or None
            if the file could not be opened.
        """
        try:
            with open(path, "rb") as fd:
                data = self._fromBytes(fd.read(), pytype)
        except FileNotFoundError:
            data = None

        return data

    def _fromBytes(self, serializedDataset: bytes, pytype: Optional[Type[Any]] = None) -> Any:
        """Read the bytes object as a python object.

        Parameters
        ----------
        serializedDataset : `bytes`
            Bytes object to unserialize.
        pytype : `type`, optional
            Expected python type to be returned.

        Returns
        -------
        inMemoryDataset : `lsst.afw.image.Filter`
            The requested data as an object.
        """
        data = yaml.load(serializedDataset, Loader=yaml.SafeLoader)

        if pytype is None:
            pytype = Filter

        # This will be a simple dict so we need to convert it to
        # the Filter type -- just needs the name
        filter = pytype(data["canonicalName"], force=True)

        return filter

    def _writeFile(self, inMemoryDataset: Any) -> None:
        """Write the in memory dataset to file on disk.

        Parameters
        ----------
        inMemoryDataset : `lsst.afw.image.Filter`
            Filter to serialize.

        Raises
        ------
        Exception
            Raised if the file could not be written or the dataset could not be
            serialized.
        """
        with open(self.fileDescriptor.location.path, "wb") as fd:
            fd.write(self._toBytes(inMemoryDataset))

    def _toBytes(self, inMemoryDataset: Any) -> bytes:
        """Write the in memory dataset to a bytestring.

        Parameters
        ----------
        inMemoryDataset : `lsst.afw.image.Filter`
            Object to serialize.

        Returns
        -------
        serializedDataset : `bytes`
            YAML string encoded to bytes.

        Raises
        ------
        Exception
            Raised if the object could not be serialized.
        """

        # Convert the Filter to a dict for dumping
        # Given the singleton situation, only the name is really
        # needed but it does not hurt to put some detail in the file
        # to aid debugging.
        filter = {}
        filter["canonicalName"] = inMemoryDataset.getCanonicalName()
        filter["name"] = inMemoryDataset.getName()
        filter["aliases"] = inMemoryDataset.getAliases()

        return yaml.dump(filter).encode()
