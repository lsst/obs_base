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

__all__ = ("FitsRawFormatterBase",)

import logging
from abc import abstractmethod
from typing import TYPE_CHECKING, Any, ClassVar, Self

from astro_metadata_translator import ObservationInfo, fix_header

import lsst.afw.fits
import lsst.afw.geom
import lsst.afw.image
from lsst.daf.butler import FileDescriptor, FormatterNotImplementedError, Location, StorageClass
from lsst.resources import ResourcePath
from lsst.utils.classes import cached_getter

from .formatters.fitsExposure import FitsImageFormatterBase, standardizeAmplifierParameters
from .makeRawVisitInfoViaObsInfo import MakeRawVisitInfoViaObsInfo
from .utils import InitialSkyWcsError, createInitialSkyWcsFromBoresight

if TYPE_CHECKING:
    import astro_metadata_translator

    import lsst.daf.base

    from .filters import FilterDefinitionCollection

log = logging.getLogger(__name__)


class FitsRawFormatterBase(FitsImageFormatterBase):
    """Abstract base class for reading and writing raw data to and from
    FITS files.
    """

    ReaderClass = lsst.afw.image.ImageFitsReader

    wcsFlipX: ClassVar[bool] = False
    """Control whether the WCS is flipped in the X-direction (`bool`).
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self._metadata: lsst.daf.base.PropertyList | None = None
        self._observationInfo: ObservationInfo | None = None

    @classmethod
    def fromMetadata(
        cls,
        metadata: lsst.daf.base.PropertyList,
        obsInfo: ObservationInfo | None = None,
        storageClass: StorageClass | None = None,
        location: Location | None = None,
    ) -> Self:
        """Construct a possibly-limited formatter from known metadata.

        Parameters
        ----------
        metadata : `lsst.daf.base.PropertyList`
            Raw header metadata, with any fixes (see
            `astro_metadata_translator.fix_header`) applied but nothing
            stripped.
        obsInfo : `astro_metadata_translator.ObservationInfo`, optional
            Structured information already extracted from ``metadata``.
            If not provided, will be read from ``metadata`` on first use.
        storageClass : `lsst.daf.butler.StorageClass`, optional
            StorageClass for this file.  If not provided, the formatter will
            only support `makeWcs`, `makeVisitInfo`, `makeFilterLabel`, and
            other operations that operate purely on metadata and not the actual
            file.
        location : `lsst.daf.butler.Location`, optional.
            Location of the file.  If not provided, the formatter will only
            support `makeWcs`, `makeVisitInfo`, `makeFilterLabel`, and other
            operations that operate purely on metadata and not the actual file.

        Returns
        -------
        formatter : `FitsRawFormatterBase`
            An instance of ``cls``.
        """
        self = cls(
            FileDescriptor(
                location if location is not None else Location("", "<unspecified>"),
                storageClass if storageClass is not None else StorageClass(),
            )
        )
        self._metadata = metadata
        self._observationInfo = obsInfo
        return self

    @property
    @abstractmethod
    def translatorClass(self) -> type[astro_metadata_translator.MetadataTranslator] | None:
        """`~astro_metadata_translator.MetadataTranslator` to translate
        metadata header to `~astro_metadata_translator.ObservationInfo`.
        """
        return None

    @property
    @abstractmethod
    def filterDefinitions(self) -> FilterDefinitionCollection | None:
        """`~lsst.obs.base.FilterDefinitionCollection`, defining the filters
        for this instrument.
        """
        return None

    @property
    @cached_getter
    def checked_parameters(self) -> dict[str, Any]:
        # Docstring inherited.
        parameters = super().checked_parameters
        if "bbox" in parameters:
            raise TypeError(
                "Raw formatters do not support reading arbitrary subimages, as some "
                "implementations may be assembled on-the-fly."
            )
        return parameters

    def readImage(self) -> lsst.afw.image.Image:
        """Read just the image component of the Exposure.

        Returns
        -------
        image : `~lsst.afw.image.Image`
            In-memory image component.
        """
        return lsst.afw.image.ImageU(self._reader_path)

    def isOnSky(self) -> bool:
        """Boolean to determine if the exposure is thought to be on the sky.

        Returns
        -------
        onSky : `bool`
            Returns `True` if the observation looks like it was taken on the
            sky.  Returns `False` if this observation looks like a calibration
            observation.

        Notes
        -----
        If there is tracking RA/Dec information associated with the
        observation it is assumed that the observation is on sky.
        Currently the observation type is not checked.
        """
        if self.observationInfo.tracking_radec is None:
            return False
        return True

    @property
    def metadata(self) -> lsst.daf.base.PropertyList:
        """The metadata read from this file. It will be stripped as
        components are extracted from it
        (`lsst.daf.base.PropertyList`).
        """
        if self._metadata is None:
            self._metadata = self.readMetadata()
        return self._metadata

    def readMetadata(self) -> lsst.daf.base.PropertyList:
        """Read all header metadata directly into a PropertyList.

        Returns
        -------
        metadata : `~lsst.daf.base.PropertyList`
            Header metadata.
        """
        md = self.reader.readMetadata()
        translatorClass = self.translatorClass
        assert translatorClass is not None
        fix_header(md, translator_class=translatorClass)
        return md

    def stripMetadata(self) -> None:
        """Remove metadata entries that are parsed into components."""
        try:
            lsst.afw.geom.stripWcsMetadata(self.metadata)
        except TypeError as e:
            log.debug("Error caught and ignored while stripping metadata: %s", e.args[0])

    def makeVisitInfo(self) -> lsst.afw.image.VisitInfo:
        """Construct a VisitInfo from metadata.

        Returns
        -------
        visitInfo : `~lsst.afw.image.VisitInfo`
            Structured metadata about the observation.
        """
        visit_info = MakeRawVisitInfoViaObsInfo.observationInfo2visitInfo(self.observationInfo)
        if self.data_id and "exposure" in self.data_id:
            # Special case exposure existing for this dataset type.
            # Want to ensure that the id stored in the visitInfo matches that
            # from the dataId from butler, regardless of what may have come
            # from the ObservationInfo. In some edge cases they might differ.
            exposure_id = self.data_id["exposure"]
            if exposure_id != visit_info.id:
                visit_info = visit_info.copyWith(id=exposure_id)

        return visit_info

    @abstractmethod
    def getDetector(self, id: int) -> lsst.afw.cameraGeom.Detector:
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

    def makeWcs(
        self, visitInfo: lsst.afw.image.VisitInfo, detector: lsst.afw.cameraGeom.Detector
    ) -> lsst.afw.geom.SkyWcs:
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
        if not self.isOnSky():
            # This is not an on-sky observation
            return None

        if visitInfo is None:
            msg = "No VisitInfo; cannot access boresight information. Defaulting to metadata-based SkyWcs."
            log.warning(msg)
            skyWcs = self._createSkyWcsFromMetadata()
            if skyWcs is None:
                raise InitialSkyWcsError(
                    "Failed to create both metadata and boresight-based SkyWcs."
                    "See warnings in log messages for details."
                )
            return skyWcs

        return self.makeRawSkyWcsFromBoresight(
            visitInfo.getBoresightRaDec(), visitInfo.getBoresightRotAngle(), detector
        )

    @classmethod
    def makeRawSkyWcsFromBoresight(
        cls,
        boresight: lsst.geom.SpherePoint,
        orientation: lsst.geom.Angle,
        detector: lsst.afw.cameraGeom.Detector,
    ) -> lsst.afw.geom.SkyWcs:
        """Class method to make a raw sky WCS from boresight and detector.

        Parameters
        ----------
        boresight : `lsst.geom.SpherePoint`
            The ICRS boresight RA/Dec
        orientation : `lsst.geom.Angle`
            The rotation angle of the focal plane on the sky.
        detector : `lsst.afw.cameraGeom.Detector`
            Where to get the camera geometry from.

        Returns
        -------
        skyWcs : `~lsst.afw.geom.SkyWcs`
            Reversible mapping from pixel coordinates to sky coordinates.
        """
        return createInitialSkyWcsFromBoresight(boresight, orientation, detector, flipX=cls.wcsFlipX)

    def _createSkyWcsFromMetadata(self) -> lsst.afw.geom.SkyWcs | None:
        """Create a SkyWcs from the FITS header metadata in an Exposure.

        Returns
        -------
        skyWcs: `lsst.afw.geom.SkyWcs`, or None
            The WCS that was created from ``self.metadata``, or None if that
            creation fails due to invalid metadata.
        """
        if not self.isOnSky():
            # This is not an on-sky observation
            return None

        try:
            return lsst.afw.geom.makeSkyWcs(self.metadata, strip=True)
        except TypeError as e:
            log.warning("Cannot create a valid WCS from metadata: %s", e.args[0])
            return None

    def makeFilterLabel(self) -> lsst.afw.image.FilterLabel:
        """Construct a FilterLabel from metadata.

        Returns
        -------
        filter : `~lsst.afw.image.FilterLabel`
            Object that identifies the filter for this image.
        """
        physical = self.observationInfo.physical_filter
        filter_definitions = self.filterDefinitions
        assert filter_definitions is not None
        band = filter_definitions.physical_to_band[physical]
        return lsst.afw.image.FilterLabel(physical=physical, band=band)

    def readComponent(self, component: str) -> Any:
        # Docstring inherited.
        _ = self.checked_parameters  # just for checking; no supported parameters.
        if component == "image":
            return self.readImage()
        elif component == "filter":
            return self.makeFilterLabel()
        elif component == "visitInfo":
            return self.makeVisitInfo()
        elif component == "detector":
            return self.getDetector(self.observationInfo.detector_num)
        elif component == "wcs":
            detector = self.getDetector(self.observationInfo.detector_num)
            visitInfo = self.makeVisitInfo()
            return self.makeWcs(visitInfo, detector)
        elif component == "metadata":
            self.stripMetadata()
            return self.metadata
        return None

    def readFull(self) -> Any:
        # Docstring inherited.
        amplifier, detector, _ = standardizeAmplifierParameters(
            self.checked_parameters,
            self.getDetector(self.observationInfo.detector_num),
        )
        if amplifier is not None:
            amplifier_isolator = lsst.afw.cameraGeom.AmplifierIsolator(
                amplifier,
                self.reader.readBBox(),
                detector,
            )
            subimage = amplifier_isolator.transform_subimage(
                self.reader.read(bbox=amplifier_isolator.subimage_bbox)
            )
            exposure = lsst.afw.image.makeExposure(lsst.afw.image.makeMaskedImage(subimage))
            exposure.setDetector(amplifier_isolator.make_detector())
        else:
            exposure = lsst.afw.image.makeExposure(lsst.afw.image.makeMaskedImage(self.readImage()))
            exposure.setDetector(detector)
        self.attachComponentsFromMetadata(exposure)
        return exposure

    def write_local_file(self, in_memory_dataset: Any, uri: ResourcePath) -> None:
        raise FormatterNotImplementedError("Raw data cannot be `put`.")

    @property
    def observationInfo(self) -> ObservationInfo:
        """The `~astro_metadata_translator.ObservationInfo` extracted from
        this file's metadata (`~astro_metadata_translator.ObservationInfo`,
        read-only).
        """
        if self._observationInfo is None:
            # Use the primary path rather than any local variant that the
            # formatter might be using.
            location = self.file_descriptor.location
            path = location.path if location is not None else None
            self._observationInfo = ObservationInfo(
                self.metadata, translator_class=self.translatorClass, filename=path
            )
        return self._observationInfo

    def attachComponentsFromMetadata(self, exposure: lsst.afw.image.Exposure) -> None:
        """Attach all `lsst.afw.image.Exposure` components derived from
        metadata (including the stripped metadata itself).

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to attach components to (modified in place).  Must already
            have a detector attached.
        """
        info = exposure.getInfo()
        info.id = self.observationInfo.detector_exposure_id
        info.setFilter(self.makeFilterLabel())
        info.setVisitInfo(self.makeVisitInfo())
        info.setWcs(self.makeWcs(info.getVisitInfo(), info.getDetector()))

        self.stripMetadata()
        exposure.setMetadata(self.metadata)
