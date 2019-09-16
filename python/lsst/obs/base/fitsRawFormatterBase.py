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
import lsst.afw.geom
from lsst.daf.butler.formatters.fitsExposureFormatter import FitsExposureFormatter
import lsst.log

from lsst.obs.base import MakeRawVisitInfoViaObsInfo
from lsst.obs.base.utils import createInitialSkyWcs, InitialSkyWcsError


class FitsRawFormatterBase(FitsExposureFormatter, metaclass=ABCMeta):
    """Abstract base class for reading and writing raw data to and from
    FITS files.
    """

    def __init__(self, *args, **kwargs):
        self.filterDefinitions.defineFilters()
        super().__init__(*args, **kwargs)

    @property
    @abstractmethod
    def translatorClass(self):
        """`~astro_metadata_translator.MetadataTranslator` to translate
        metadata header to `~astro_metadata_translator.ObservationInfo`.
        """
        return None

    _observationInfo = None

    @property
    @abstractmethod
    def filterDefinitions(self):
        """`~lsst.obs.base.FilterDefinitions`, defining the filters for this
        instrument.
        """
        return None

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
        # NOTE: makeVisitInfo() may not strip any metadata itself, but calling
        # it ensures that ObservationInfo is created from the metadata, which
        # will strip the VisitInfo keys and more.
        self.makeVisitInfo()
        self._createSkyWcsFromMetadata()

    def makeVisitInfo(self):
        """Construct a VisitInfo from metadata.

        Returns
        -------
        visitInfo : `~lsst.afw.image.VisitInfo`
            Structured metadata about the observation.
        """
        return MakeRawVisitInfoViaObsInfo.observationInfo2visitInfo(self.observationInfo)

    @abstractmethod
    def getDetector(self, id):
        """Return the detector that acquired this raw exposure.

        Parameters
        ----------
        id : `int`
            The identifying number of the detector to get.

        Returns
        -------
        detector : `~lsst.afw.cameraGeom.Detector`
            The detector associated with that ``id``.
        """
        raise NotImplementedError("Must be implemented by subclasses.")

    def makeWcs(self, visitInfo, detector):
        """Create a SkyWcs from information about the exposure.

        If VisitInfo is not None, use it and the detector to create a SkyWcs,
        otherwise return the metadata-based SkyWcs (always created, so that
        the relevant metadata keywords are stripped).

        Parameters
        ----------
        visitInfo : `~lsst.afw.image.VisitInfo`
            The information about the telescope boresight and camera
            orientation angle for this exposure.
        detector : `~lsst.afw.cameraGeom.Detector`
            The detector used to acquire this exposure.

        Returns
        -------
        skyWcs : `~lsst.afw.geom.SkyWcs`
            Reversible mapping from pixel coordinates to sky coordinates.

        Raises
        ------
        InitialSkyWcsError
            Raised if there is an error generating the SkyWcs, chained from the
            lower-level exception if available.
        """
        if self.observationInfo.tracking_radec is None:
            # This is not an on-sky observation
            return None

        skyWcs = self._createSkyWcsFromMetadata()

        log = lsst.log.Log.getLogger("fitsRawFormatter")
        if visitInfo is None:
            msg = "No VisitInfo; cannot access boresight information. Defaulting to metadata-based SkyWcs."
            log.warn(msg)
            if skyWcs is None:
                raise InitialSkyWcsError("Failed to create both metadata and boresight-based SkyWcs."
                                         "See warnings in log messages for details.")
            return skyWcs
        skyWcs = createInitialSkyWcs(visitInfo, detector)

        return skyWcs

    def _createSkyWcsFromMetadata(self):
        """Create a SkyWcs from the FITS header metadata in an Exposure.

        Returns
        -------
        skyWcs: `lsst.afw.geom.SkyWcs`, or None
            The WCS that was created from ``self.metadata``, or None if that
            creation fails due to invalid metadata.
        """
        if self.observationInfo.tracking_radec is None:
            # This is not an on-sky observation
            return None

        try:
            return lsst.afw.geom.makeSkyWcs(self.metadata, strip=True)
        except TypeError as e:
            log = lsst.log.Log.getLogger("fitsRawFormatter")
            log.warn("Cannot create a valid WCS from metadata: %s", e.args[0])
            return None

    def makeFilter(self):
        """Construct a Filter from metadata.

        Returns
        -------
        filter : `~lsst.afw.image.Filter`
            Object that identifies the filter for this image.

        Raises
        ------
        NotFoundError
            Raised if the physical filter was not registered via
            `~lsst.afw.image.utils.defineFilter`.
        """
        return lsst.afw.image.Filter(self.observationInfo.physical_filter)

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
            detector = self.getDetector(self.observationInfo.detector_num)
            visitInfo = self.makeVisitInfo()
            return self.makeWcs(visitInfo, detector)
        return None

    def readFull(self, parameters=None):
        """Read the full Exposure object.

        Parameters
        ----------
        parameters : `dict`, optional
            If specified, a dictionary of slicing parameters that overrides
            those in the `fileDescriptor` attribute.

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
        full.setDetector(self.getDetector(self.observationInfo.detector_num))
        info = full.getInfo()
        info.setFilter(self.makeFilter())
        info.setVisitInfo(self.makeVisitInfo())
        info.setWcs(self.makeWcs(info.getVisitInfo(), info.getDetector()))
        # We don't need to call stripMetadata() here because it has already
        # been stripped during creation of the ObservationInfo, WCS, etc.
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
        """The `~astro_metadata_translator.ObservationInfo` extracted from
        this file's metadata (`~astro_metadata_translator.ObservationInfo`,
        read-only).
        """
        if self._observationInfo is None:
            self._observationInfo = ObservationInfo(self.metadata, translator_class=self.translatorClass)
        return self._observationInfo
