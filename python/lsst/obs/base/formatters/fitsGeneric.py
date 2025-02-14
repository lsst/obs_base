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

__all__ = ("FitsGenericFormatter",)

import importlib.metadata
from typing import Any

from packaging import version

from lsst.daf.butler import DatasetProvenance, FormatterV2
from lsst.resources import ResourcePath

from ..utils import add_provenance_to_fits_header

_astropy_version = version.Version(importlib.metadata.version("astropy"))
_ALLOW_LONG_FITS_STRINGS = True
if _astropy_version < version.Version("7.1.0"):
    # Astropy cannot handle long values for HIERARCH strings.
    _ALLOW_LONG_FITS_STRINGS = False


class FitsGenericFormatter(FormatterV2):
    """Interface for reading and writing objects that support the standard
    afw I/O readFits/writeFits methods.
    """

    supported_extensions = frozenset({".fits", ".fits.gz", ".fits.fz", ".fz", ".fit"})
    default_extension = ".fits"
    can_read_from_local_file = True

    def read_from_local_file(self, path: str, component: str | None = None, expected_size: int = -1) -> Any:
        pytype = self.file_descriptor.storageClass.pytype
        if self.file_descriptor.parameters:
            try:
                return pytype.readFitsWithOptions(  # type: ignore
                    path, options=self.file_descriptor.parameters
                )
            except AttributeError:
                pass

        return pytype.readFits(path)  # type: ignore

    def write_local_file(self, in_memory_dataset: Any, uri: ResourcePath) -> None:
        in_memory_dataset.writeFits(uri.ospath)

    def add_provenance(
        self, in_memory_dataset: Any, /, *, provenance: DatasetProvenance | None = None
    ) -> Any:
        if hasattr(in_memory_dataset, "getMetadata"):
            metadata = in_memory_dataset.getMetadata()
        elif hasattr(in_memory_dataset, "metadata"):
            metadata = in_memory_dataset.metadata
        else:
            # Unable to find compliant metadata attribute.
            return in_memory_dataset
        # Sometimes astropy is used to write FITS headers and astropy
        # cannot write long HIERARCH values. Therefore always elide some
        # of the provenance just in case astropy is used.
        add_provenance_to_fits_header(
            metadata, self.dataset_ref, provenance, allow_long_headers=_ALLOW_LONG_FITS_STRINGS
        )
        return in_memory_dataset
