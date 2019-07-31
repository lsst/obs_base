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

__all__ = ("FitsRawFormatterBase",)

from abc import ABCMeta, abstractmethod

from astro_metadata_translator import ObservationInfo

import lsst.afw.image
from lsst.daf.butler.formatters.fitsExposureFormatter import FitsExposureFormatter
from lsst.obs.base import MakeRawVisitInfoViaObsInfo


class FitsRawFormatterBase(FitsExposureFormatter, metaclass=ABCMeta):
    """Abstract base class for reading and writing raw data to and from
    FITS files.
    """

    @property
    @abstractmethod
    def translatorClass(self):
        """MetadataTranslator to translate metadata header to ObservationInfo.
        """
        return None

    _observationInfo = None
    _metadata = None

    def readImage(self):
        """Read just the image component of the Exposure.

        Returns
        -------
        image : `~lsst.afw.image.Image`
            In-memory image component.
        """
        return lsst.afw.image.ImageU(self.fileDescriptor.location.path)

    def readMask(self):
        """Read just the mask component of the Exposure.

        May return None (as the default implementation does) to indicate that
        there is no mask information to be extracted (at least not trivially)
        from the raw data.  This will prohibit direct reading of just the mask,
        and set the mask of the full Exposure to zeros.

        Returns
        -------
        mask : `~lsst.afw.image.Mask`
            In-memory mask component.
        """
        return None

    def readVariance(self):
        """Read just the variance component of the Exposure.

        May return None (as the default implementation does) to indicate that
        there is no variance information to be extracted (at least not
        trivially) from the raw data.  This will prohibit direct reading of
        just the variance, and set the variance of the full Exposure to zeros.

        Returns
        -------
        image : `~lsst.afw.image.Image`
            In-memory variance component.
        """
        return None

    def stripMetadata(self):
        """Remove metadata entries that are parsed into components.
        """
        self.makeVisitInfo(self.metadata)
        self.makeWcs(self.metadata)

    def makeVisitInfo(self):
        """Construct a VisitInfo from ObservationInfo.

        Returns
        -------
        visitInfo : `~lsst.afw.image.VisitInfo`
            Structured metadata about the observation.
        """
        return MakeRawVisitInfoViaObsInfo.observationInfo2visitInfo(self.observationInfo)

    def makeWcs(self):
        """Construct a SkyWcs from metadata.

        Returns
        -------
        wcs : `~lsst.afw.geom.SkyWcs`
            Reversible mapping from pixel coordinates to sky coordinates.
        """
        from lsst.afw.geom import makeSkyWcs
        return makeSkyWcs(self.metadata, strip=True)

    def makeFilter(self):
        """Construct a Filter from metadata.

        Returns
        -------
        filter : `~lsst.afw.image.Filter`
            Object that identifies the filter for this image.
        """
        raise NotImplementedError("Must be implemented by subclasses.")

    def readImageComponent(self, component):
        """Read the image, mask, or variance component of an Exposure.

        Parameters
        ----------
        component : `str`, optional
            Component to read from the file.  Always one of "image",
            "variance", or "mask".

        Returns
        -------
        image : `~lsst.afw.image.Image` or `~lsst.afw.image.Mask`
            In-memory image, variance, or mask component.
        """
        if component == "image":
            return self.readImage()
        elif component == "mask":
            return self.readMask()
        elif component == "variance":
            return self.readVariance()

    def readInfoComponent(self, component):
        """Read a component held by ExposureInfo.

        The implementation provided by FitsRawFormatter provides only "wcs"
        and "visitInfo".  When adding support for other components, subclasses
        should delegate to `super()` for those and update `readFull` with
        similar logic.

        Parameters
        ----------
        component : `str`, optional
            Component to read from the file.

        Returns
        -------
        obj : component-dependent
            In-memory component object.
        """
        if component == "filter":
            return self.makeFilter()
        elif component == "visitInfo":
            return self.makeVisitInfo()
        elif component == "wcs":
            return self.makeWcs()
        return None

    def readFull(self, parameters=None):
        """Read the full Exposure object.

        Parameters
        ----------
        parameters : `dict`, optional
            If specified a dictionary of slicing parameters that overrides
            those in ``self.fileDescriptor`.

        Returns
        -------
        exposure : `~lsst.afw.image.Exposure`
            Complete in-memory exposure.
        """
        from lsst.afw.image import makeExposure, makeMaskedImage
        full = makeExposure(makeMaskedImage(self.readImage()))
        mask = self.readMask()
        if mask is not None:
            full.setMask(mask)
        variance = self.readVariance()
        if variance is not None:
            full.setVariance(variance)
        info = full.getInfo()
        info.setWcs(self.makeWcs())
        info.setFilter(self.makeFilter())
        info.setVisitInfo(self.makeVisitInfo())
        # We shouldn't have to call stripMetadata() here because that should
        # have been done by makeVisitInfo and makeWcs (or by subclasses that
        # strip metadata for other components when constructing them).
        full.setMetadata(self.metadata)
        return full

    def write(self, inMemoryDataset):
        """Write a Python object to a file.

        Parameters
        ----------
        inMemoryDataset : `object`
            The Python object to store.

        Returns
        -------
        path : `str`
            The `URI` where the primary file is stored.
        """
        raise NotImplementedError("Raw data cannot be `put`.")

    @property
    def observationInfo(self):
        """The ObservationInfo extracted from this raw file's metadata.
        """
        if self._observationInfo is None:
            self._observationInfo = ObservationInfo(self._metadata, translator_class=self.translatorClass)
        return self._observationInfo

    @property
    def metadata(self):
        """The metadata read from this raw file. It will be stripped as
        components are extracted from it.
        """
        if self._metadata is None:
            self._metadata = self.readMetadata()
        return self._metadata
