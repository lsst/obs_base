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

__all__ = (
    "FitsExposureFormatter",
    "FitsImageFormatter",
    "FitsImageFormatterBase",
    "FitsMaskFormatter",
    "FitsMaskedImageFormatter",
    "standardizeAmplifierParameters",
)

import hashlib
import json
import logging
import threading
import uuid
import warnings
from abc import abstractmethod
from collections.abc import Mapping, Set
from io import BytesIO
from typing import TYPE_CHECKING, Any, ClassVar, NamedTuple, Protocol

import astropy.io.fits
import numpy as np

import lsst.geom
from lsst.afw.cameraGeom import AmplifierGeometryComparison, AmplifierIsolator
from lsst.afw.fits import CompressionOptions, MemFileManager
from lsst.afw.geom.wcsUtils import getImageXY0FromMetadata
from lsst.afw.image import (
    ExposureFitsReader,
    ExposureInfo,
    FilterLabel,
    ImageFitsReader,
    MaskedImageFitsReader,
    MaskFitsReader,
)

# Needed for ApCorrMap to resolve properly
from lsst.afw.math import BoundedField  # noqa: F401
from lsst.daf.base import PropertyList
from lsst.daf.butler import DatasetProvenance, FormatterV2
from lsst.resources import ResourcePath
from lsst.utils.classes import cached_getter
from lsst.utils.introspection import find_outside_stacklevel

from ..utils import add_provenance_to_fits_header

if TYPE_CHECKING:
    import lsst.afw.cameraGeom


_LOG = logging.getLogger(__name__)

_ALWAYS_USE_ASTROPY_FOR_COMPONENT_READ = False
"""If True, the astropy code will always be used to read component and cutouts
even if the file is local, the cutout is too large, or the dataset type is
wrong. This should mostly be used for testing.
"""


class _ReaderClassLike(Protocol):
    def __init__(self, path: str) -> None: ...
    def readBBox(self) -> lsst.geom.Box2I: ...
    def read(self, bbox: lsst.geom.Box2I = lsst.geom.Box2I(), dtype: Any = None) -> Any: ...
    def readImage(self, bbox: lsst.geom.Box2I = lsst.geom.Box2I(), dtype: Any = None) -> Any: ...
    def readMask(self, bbox: lsst.geom.Box2I = lsst.geom.Box2I(), dtype: Any = None) -> Any: ...
    def readVariance(self, bbox: lsst.geom.Box2I = lsst.geom.Box2I(), dtype: Any = None) -> Any: ...
    def readDetector(self) -> lsst.afw.cameraGeom.Detector: ...
    def readComponent(self, component: str) -> Any: ...
    def readMetadata(self) -> PropertyList: ...
    def readSerializationVersion(self) -> int: ...


class FitsImageFormatterBase(FormatterV2):
    """Base class formatter for image-like storage classes stored via FITS.

    Notes
    -----
    This class makes no assumptions about how many HDUs are used to represent
    the image on disk, and includes no support for writing.  It's really just a
    collection of miscellaneous boilerplate common to all FITS image
    formatters.

    Concrete subclasses must implement `readComponent`, `readFull`, and `write`
    (even if just to disable them by raising an exception).
    """

    can_read_from_local_file = True
    default_extension = ".fits"
    supported_extensions: ClassVar[Set[str]] = frozenset({".fits", ".fits.gz", ".fits.fz", ".fz", ".fit"})

    unsupported_parameters: ClassVar[Set[str]] = frozenset()
    """Support all parameters."""

    _reader = None
    _reader_path: str | None = None

    ReaderClass: type[_ReaderClassLike]  # must be set by concrete subclasses

    @property
    def reader(self) -> _ReaderClassLike:
        """The reader object that backs this formatter's read operations.

        This is computed on first use and then cached.  It should never be
        accessed when writing. Currently assumes a local file.
        """
        if self._reader is None:
            if self._reader_path is None:
                raise RuntimeError("Internal error in formatter; failing to set path.")
            self._reader = self.ReaderClass(self._reader_path)
        return self._reader

    @property
    @cached_getter
    def checked_parameters(self) -> dict[str, Any]:
        """The parameters passed by the butler user, after checking them
        against the storage class and transforming `None` into an empty `dict`
        (`dict`).

        This is computed on first use and then cached.  It should never be
        accessed when writing.  Subclasses that need additional checking should
        delegate to `super` and then check the result before returning it.
        """
        parameters = self.file_descriptor.parameters
        if parameters is None:
            parameters = {}
        self.file_descriptor.storageClass.validateParameters(parameters)
        return parameters

    @property
    def storageClass_dtype(self) -> np.dtype | None:
        """The numpy data type associated with the storage class."""
        dtype: np.dtype | None = None
        try:
            # lsst.afw.image.Exposure is generic base class and does not have
            # the dtype attribute.
            dtype = np.dtype(self.file_descriptor.storageClass.pytype.dtype)  # type: ignore[attr-defined]
        except AttributeError:
            pass
        return dtype

    def read_from_local_file(self, path: str, component: str | None = None, expected_size: int = -1) -> Any:
        # Docstring inherited.
        # The methods doing the reading all currently assume local file
        # and assume that the file descriptor refers to a local file.
        # With FormatterV2 that file descriptor does not refer to a local
        # file.
        self._reader_path = path
        self._reader = None  # Ensure the reader class is reset.
        try:
            if component is not None:
                in_memory_dataset = self.readComponent(component)
            else:
                in_memory_dataset = self.readFull()
        finally:
            self._reader = None  # Release the file handle.
        return in_memory_dataset

    @abstractmethod
    def readComponent(self, component: str) -> Any:
        """Read a component dataset.

        Parameters
        ----------
        component : `str`, optional
            Component to read from the file.

        Returns
        -------
        obj : component-dependent
            In-memory component object.

        Raises
        ------
        KeyError
            Raised if the requested component cannot be handled.
        """
        raise NotImplementedError()

    @abstractmethod
    def readFull(self) -> Any:
        """Read the full dataset (while still accounting for parameters).

        Returns
        -------
        obj : component-dependent
            In-memory component object.

        """
        raise NotImplementedError()


class ReaderFitsImageFormatterBase(FitsImageFormatterBase):
    """Base class formatter for image-like storage classes stored via FITS
    backed by a "reader" object similar to `lsst.afw.image.ImageFitsReader`.

    Notes
    -----
    This class includes no support for writing.

    Concrete subclasses must provide at least the `ReaderClass` attribute
    and a `write` implementation (even just to disable writing by raising).

    The provided implementation of `readComponent` handles only the 'bbox',
    'dimensions', and 'xy0' components common to all image-like storage
    classes.  Subclasses with additional components should handle them first,
    then delegate to ``super()`` for these (or, if necessary, delegate first
    and catch `KeyError`).

    The provided implementation of `readFull` handles only parameters that
    can be forwarded directly to the reader class (usually ``bbox`` and
    ``origin``).  Concrete subclasses that need to handle additional parameters
    should generally reimplement without delegating (the implementation is
    trivial).
    """


class StandardFitsImageFormatterBase(ReaderFitsImageFormatterBase):
    """Base class interface for image-like storage stored via FITS,
    written using LSST code.

    Notes
    -----
    Concrete subclasses must provide at least the `ReaderClass` attribute.

    The provided implementation of `readComponent` handles only the 'bbox',
    'dimensions', and 'xy0' components common to all image-like storage
    classes.  Subclasses with additional components should handle them first,
    then delegate to ``super()`` for these (or, if necessary, delegate first
    and catch `KeyError`).

    The provided implementation of `readFull` handles only parameters that
    can be forwarded directly to the reader class (usually ``bbox`` and
    ``origin``).  Concrete subclasses that need to handle additional parameters
    should generally reimplement without delegating (the implementation is
    trivial).

    This Formatter supports write recipes, and assumes its in-memory type has
    ``writeFits`` and (for write recipes) ``writeFitsWithOptions`` methods.

    Each ``StandardFitsImageFormatterBase`` recipe for FITS compression should
    define ``image``, ``mask`` and ``variance`` entries, each of which may
    contain entries supported by
    `lsst.afw.fits.CompressionOptions.from_mapping` (``null`` disables
    compression).

    A very simple example YAML recipe (for the ``Exposure`` specialization):

    .. code-block:: yaml

        lsst.obs.base.fitsExposureFormatter.FitsExposureFormatter:
          default:
            image: &default
              algorithm: GZIP_2
            mask: *default
            variance: *default

    """

    supported_write_parameters = frozenset({"recipe"})

    def readComponent(self, component: str) -> Any:
        # Docstring inherited.
        if component in ("bbox", "dimensions", "xy0"):
            bbox = self.reader.readBBox()
            if component == "dimensions":
                return bbox.getDimensions()
            elif component == "xy0":
                return bbox.getMin()
            else:
                return bbox
        else:
            raise KeyError(f"Unknown component requested: {component}")

    def readFull(self) -> Any:
        # Docstring inherited.
        return self.reader.read(**self.checked_parameters, dtype=self.storageClass_dtype)

    def write_local_file(self, in_memory_dataset: Any, uri: ResourcePath) -> None:
        # check to see if we have a recipe requested
        recipeName = self.write_parameters.get("recipe")
        recipe = self.get_image_compression_settings(recipeName)
        if recipe:
            in_memory_dataset.writeFitsWithOptions(uri.ospath, options=recipe)
        else:
            in_memory_dataset.writeFits(uri.ospath)

    def get_image_compression_settings(self, recipeName: str | None) -> dict:
        """Retrieve the relevant compression settings for this recipe.

        Parameters
        ----------
        recipeName : `str` or `None`
            Label associated with the collection of compression parameters
            to select.

        Returns
        -------
        settings : `dict`
            The selected settings.
        """
        # if no recipe has been provided and there is no default
        # return immediately
        if not recipeName:
            if "default" not in self.write_recipes:
                return {}
            recipeName = "default"

        if recipeName not in self.write_recipes:
            raise RuntimeError(f"Unrecognized recipe option given for compression: {recipeName}")

        recipe = self.write_recipes[recipeName]
        if recipe is None:
            return {}
        seed: int | None = None
        for plane in ("image", "mask", "variance"):
            if plane in recipe and (quantization := recipe[plane].get("quantization")) is not None:
                if quantization.get("seed", 0) == 0:
                    if seed is None:
                        # Set the seed based on data ID.  We can't just use
                        # 'hash', since like 'set' that's not deterministic.
                        # And we can't rely on a DimensionPacker because those
                        # are only defined for certain combinations of
                        # dimensions.  Doing an MD5 of the JSON feels like
                        # overkill but I don't really see anything much
                        # simpler.
                        hash_bytes = hashlib.md5(
                            json.dumps(list(self.data_id.required_values)).encode(),
                            usedforsecurity=False,
                        ).digest()
                        # And it *really* feels like overkill when we squash
                        # that into the [1, 10000] range allowed by FITS.
                        seed = 1 + int.from_bytes(hash_bytes) % 9999
                    _LOG.debug(
                        "Setting compression quantization seed for %s %s %s to %s.",
                        self.data_id,
                        self.dataset_ref.datasetType.name,
                        plane,
                        seed,
                    )
                    quantization["seed"] = seed
                else:
                    _LOG.warning(
                        "Compression quantization seed for %s %s %s was set explicitly to %s.",
                        self.dataset_ref.datasetType.name,
                        self.data_id,
                        plane,
                        quantization["seed"],
                    )
            else:
                _LOG.debug(
                    "No quantization found for %s %s %s.",
                    self.dataset_ref.datasetType.name,
                    self.data_id,
                    plane,
                )
        return recipe

    @classmethod
    def validate_write_recipes(cls, recipes: Mapping[str, Any] | None) -> Mapping[str, Any] | None:
        """Validate supplied recipes for this formatter.

        The recipes are supplemented with default values where appropriate.

        Parameters
        ----------
        recipes : `dict` or `None`
            Recipes to validate. Can be empty dict or `None`.

        Returns
        -------
        validated : `dict`
            Validated recipes. Returns what was given if there are no
            recipes listed.

        Raises
        ------
        RuntimeError
            Raised if validation fails.
        """
        if not recipes:
            # We can not insist on recipes being specified.
            return recipes

        validated: dict[str, Any] = {}
        for name, recipe in recipes.items():
            if recipe is not None:
                validated[name] = {}
                for plane in ["image", "mask", "variance"]:
                    try:
                        options = CompressionOptions.from_mapping(recipe[plane])
                    except Exception as err:
                        err.add_note(f"Validating write recipe {name!r} ({plane!r} section).")
                        raise
                    validated[name][plane] = options.to_dict()
            else:
                validated[name] = None
        return validated


class FitsImageFormatter(StandardFitsImageFormatterBase):
    """Concrete formatter for reading/writing `~lsst.afw.image.Image`
    from/to FITS.
    """

    ReaderClass = ImageFitsReader


class FitsMaskFormatter(StandardFitsImageFormatterBase):
    """Concrete formatter for reading/writing `~lsst.afw.image.Mask`
    from/to FITS.
    """

    ReaderClass = MaskFitsReader


class FitsMaskedImageFormatter(StandardFitsImageFormatterBase):
    """Concrete formatter for reading/writing `~lsst.afw.image.MaskedImage`
    from/to FITS.
    """

    ReaderClass = MaskedImageFitsReader

    def readComponent(self, component: str) -> Any:
        # Docstring inherited.
        if component == "image":
            return self.reader.readImage(**self.checked_parameters, dtype=self.storageClass_dtype)
        elif component == "mask":
            return self.reader.readMask(**self.checked_parameters)
        elif component == "variance":
            return self.reader.readVariance(**self.checked_parameters, dtype=self.storageClass_dtype)
        else:
            # Delegate to base for bbox, dimensions, xy0.
            return super().readComponent(component)


def standardizeAmplifierParameters(
    parameters: dict[str, Any], on_disk_detector: lsst.afw.cameraGeom.Detector | None
) -> tuple[lsst.afw.cameraGeom.Amplifier, lsst.afw.cameraGeom.Detector, bool]:
    """Preprocess the Exposure storage class's "amp" and "detector" parameters.

    This checks the given objects for consistency with the on-disk geometry and
    converts amplifier IDs/names to Amplifier instances.

    Parameters
    ----------
    parameters : `dict`
        Dictionary of parameters passed to formatter.  See the Exposure storage
        class definition in daf_butler for allowed keys and values.
    on_disk_detector : `lsst.afw.cameraGeom.Detector` or `None`
        Detector that represents the on-disk image being loaded, or `None` if
        this is unknown (and hence the user must provide one in
        ``parameters`` if "amp" is in ``parameters``).

    Returns
    -------
    amplifier : `lsst.afw.cameraGeom.Amplifier` or `None`
        An amplifier object that defines a subimage to load, or `None` if there
        was no "amp" parameter.
    detector : `lsst.afw.cameraGeom.Detector` or `None`
        A detector object whose amplifiers are in the same s/orientation
        state as the on-disk image.  If there is no "amp" parameter,
        ``on_disk_detector`` is simply passed through.
    regions_differ : `bool`
        `True` if the on-disk detector and the detector given in the parameters
        had different bounding boxes for one or more regions.  This can happen
        if the true overscan region sizes can only be determined when the image
        is actually read, but otherwise it should be considered user error.
    """
    if (amplifier := parameters.get("amp")) is None:
        return None, on_disk_detector, False
    if "bbox" in parameters or "origin" in parameters:
        raise ValueError("Cannot pass 'amp' with 'bbox' or 'origin'.")
    if isinstance(amplifier, int | str):
        amp_key = amplifier
        target_amplifier = None
    else:
        amp_key = amplifier.getName()
        target_amplifier = amplifier
    if (detector := parameters.get("detector")) is not None:
        if on_disk_detector is not None:
            # User passed a detector and we also found one on disk.  Check them
            # for consistency.  Note that we are checking the amps we'd get
            # from the two detectors against each other, not the amplifier we
            # got directly from the user, as the latter is allowed to differ in
            # assembly/orientation state.
            comparison = on_disk_detector[amp_key].compareGeometry(detector[amp_key])
            if comparison & comparison.ASSEMBLY_DIFFERS:
                raise ValueError(
                    "The given 'detector' has a different assembly state and/or orientation from "
                    f"the on-disk one for amp {amp_key}."
                )
    else:
        if on_disk_detector is None:
            raise ValueError(
                f"No on-disk detector and no detector given; cannot load amplifier from key {amp_key}. "
                "Please provide either a 'detector' parameter or an Amplifier instance in the "
                "'amp' parameter."
            )
        comparison = AmplifierGeometryComparison.EQUAL
        detector = on_disk_detector
    if target_amplifier is None:
        target_amplifier = detector[amp_key]
    return target_amplifier, detector, comparison & comparison.REGIONS_DIFFER


class _ComponentCache(NamedTuple):
    id_: uuid.UUID | None = None
    reader: ExposureFitsReader | None = None
    bbox: lsst.geom.Box2I | None = None
    mem: MemFileManager | None = None


class FitsExposureFormatter(FitsMaskedImageFormatter):
    """Concrete formatter for reading/writing `~lsst.afw.image.Exposure`
    from/to FITS.

    Notes
    -----
    This class inherits from `FitsMaskedImageFormatter` even though
    `lsst.afw.image.Exposure` doesn't inherit from
    `lsst.afw.image.MaskedImage`; this is just an easy way to be able to
    delegate to `FitsMaskedImageFormatter.super()` for component-handling, and
    should be replaced with e.g. both calling a free function if that slight
    type covariance violation ever becomes a practical problem.
    """

    can_read_from_uri = True
    ReaderClass = ExposureFitsReader
    # TODO: Remove MemFileManager from cache when DM-49640 is fixed.
    _lock = threading.Lock()
    _cached_fits: _ComponentCache = _ComponentCache()

    def read_from_uri(self, uri: ResourcePath, component: str | None = None, expected_size: int = -1) -> Any:
        # For now only support small non-pixel components. In future
        # could work with cutouts.
        self._reader = None  # Guarantee things are reset.

        # Full read, always use local file read.
        if not component and not self.checked_parameters:
            return NotImplemented

        if not _ALWAYS_USE_ASTROPY_FOR_COMPONENT_READ and uri.isLocal:
            # For a local URI allow afw to read it directly.
            return NotImplemented
        pixel_components = ("mask", "image", "variance")

        if component in pixel_components:
            # For pixel access currently this can not be cached in memory
            # and the performance gains are unclear. Assume local file
            # read with file caching for now.
            return NotImplemented

        # With current file layouts the non-pixel extensions account for 1/3
        # of the file size and it is more efficient to download the entire
        # file.
        if not (
            _ALWAYS_USE_ASTROPY_FOR_COMPONENT_READ
            or self._dataset_ref.dataId.mapping.keys().isdisjoint({"tract", "patch"})
        ):
            return NotImplemented

        # Cutouts can be optimized. For now only use this optimization
        # if bbox is the only parameter and the number of pixels in the
        # bounding box is reasonable.
        bbox = None
        origin = lsst.afw.image.PARENT
        if not component:
            # Try to support PARENT and LOCAL origin but if there are any
            # other parameters do not attempt a cutout.
            if self.checked_parameters.keys() - {"bbox", "origin"}:
                return NotImplemented
            bbox = self.checked_parameters["bbox"]
            origin = self.checked_parameters.get("origin", lsst.afw.image.PARENT)
            # For larger cutouts use the full file.
            max_cutout_size = 500 * 500
            if not _ALWAYS_USE_ASTROPY_FOR_COMPONENT_READ and bbox.width * bbox.height > max_cutout_size:
                return NotImplemented

        # We only cache component reads since those are small.
        if component:
            with self._lock:
                cache = type(self)._cached_fits
            if self.dataset_ref.id == cache.id_:
                if component in {"xy0", "dimensions", "bbox"} and cache.bbox is not None:
                    match component:
                        case "xy0":
                            return cache.bbox.getMin()
                        case "dimensions":
                            return cache.bbox.getDimensions()
                        case "bbox":
                            return cache.bbox
                else:
                    self._reader = cache.reader
                    return self.readComponent(component)

        try:
            fs, fspath = uri.to_fsspec()
        except Exception:
            # fsspec cannot be initialized, fall back to downloading the file.
            return NotImplemented

        bbox_component = None
        try:
            hdul = []
            with fs.open(fspath) as f, astropy.io.fits.open(f) as fits_obj:
                # Read all non-pixel components and cache.
                for hdu in fits_obj:
                    hdr = hdu.header
                    extname = hdr.get("EXTNAME")
                    # Older files have IMAGE in EXTNAME in PRIMARY so check
                    # for EXTEND=T.
                    extend = hdr.get("EXTEND")
                    if not extend and extname and extname.lower() in pixel_components:
                        # Calculate the dimensional components for later
                        # caching. Do not derive from cached FITS reader
                        # because they depend on the dimensionality of the
                        # pixel data and we do not want to cache the pixel
                        # data.
                        if bbox_component is None:
                            shape = hdu.shape
                            dimensions = lsst.geom.Extent2I(shape[1], shape[0])

                            # XY0 is defined in the A WCS.
                            pl = PropertyList()
                            pl.update(hdr)
                            xy0 = getImageXY0FromMetadata(pl, "A", strip=False)

                            # This is the PARENT bbox.
                            bbox_component = lsst.geom.Box2I(xy0, dimensions)

                        # Handle cutout request.
                        if bbox:
                            if origin == lsst.afw.image.PARENT:
                                full_bbox = bbox_component
                            else:
                                full_bbox = lsst.geom.Box2I(
                                    lsst.geom.Point2I(0, 0), bbox_component.getDimensions
                                )
                            minX = bbox.getBeginX() - full_bbox.getBeginX()
                            maxX = bbox.getEndX() - full_bbox.getBeginX()
                            minY = bbox.getBeginY() - full_bbox.getBeginY()
                            maxY = bbox.getEndY() - full_bbox.getBeginY()
                            data = hdu.section[minY:maxY, minX:maxX]

                            # Must correct the header WCS to take into
                            # account the offset.
                            if (k := "CRPIX1") in hdr:
                                hdr[k] -= minX
                            if (k := "CRPIX2") in hdr:
                                hdr[k] -= minY
                            if (k := "LTV1") in hdr:
                                hdr[k] = -bbox.getBeginX()
                            if (k := "LTV2") in hdr:
                                hdr[k] = -bbox.getBeginY()
                            if (k := "CRVAL1A") in hdr:
                                hdr[k] = bbox.getBeginX()
                            if (k := "CRVAL2A") in hdr:
                                hdr[k] = bbox.getBeginY()
                        else:
                            data = np.zeros([1, 1], dtype=np.int32)

                        # Construct a new HDU and copy the header.
                        stripped_hdu = astropy.io.fits.ImageHDU(data=data, header=hdr)
                        hdul.append(stripped_hdu)
                    else:
                        hdul.append(hdu)
                stripped_fits = astropy.io.fits.HDUList(hdus=hdul)
                # Write the FITS file to in-memory FITS.
                buffer = BytesIO()
                stripped_fits.writeto(buffer)
        except Exception as e:
            # For some reason we can't open the remote file so fall back.
            _LOG.debug(
                "Attempted remote read of components but encountered an error. "
                "Falling back to file download. Error: %s",
                str(e),
            )
            return NotImplemented

        # Pass the new FITS buffer to the reader class without going through
        # a temporary file. We can assume this is relatively small for
        # components.
        fits_data = buffer.getvalue()
        mem = MemFileManager(len(fits_data))
        mem.setData(fits_data, len(fits_data))
        self._reader = self.ReaderClass(mem)

        if component:
            with self._lock:
                type(self)._cached_fits = _ComponentCache(
                    id_=self.dataset_ref.id,
                    reader=self._reader,
                    mem=mem,
                    bbox=bbox_component,
                )
            match component:
                case "xy0":
                    if bbox_component is None:  # For mypy.
                        return None
                    return bbox_component.getMin()
                case "dimensions":
                    if bbox_component is None:
                        return None
                    return bbox_component.getDimensions()
                case "bbox":
                    return bbox_component
                case _:
                    return self.readComponent(component)
        else:
            # Must be a cutout. We have applied the bbox parameter so no
            # parameters should be passed here.
            cutout = self.reader.read(dtype=self.storageClass_dtype)
            cutout.getInfo().setFilter(self._fixFilterLabels(cutout.getInfo().getFilter()))
            return cutout

    def add_provenance(
        self, in_memory_dataset: Any, /, *, provenance: DatasetProvenance | None = None
    ) -> Any:
        # Add provenance via FITS headers.
        add_provenance_to_fits_header(in_memory_dataset.metadata, self.dataset_ref, provenance)
        return in_memory_dataset

    def readComponent(self, component: str) -> Any:
        # Docstring inherited.
        # Generic components can be read via a string name; DM-27754 will make
        # this mapping larger at the expense of the following one.
        genericComponents = {
            "summaryStats": ExposureInfo.KEY_SUMMARY_STATS,
        }
        if (genericComponentName := genericComponents.get(component)) is not None:
            return self.reader.readComponent(genericComponentName)
        # Other components have hard-coded method names, but don't take
        # parameters.
        standardComponents = {
            "id": "readExposureId",
            "metadata": "readMetadata",
            "wcs": "readWcs",
            "coaddInputs": "readCoaddInputs",
            "psf": "readPsf",
            "photoCalib": "readPhotoCalib",
            "filter": "readFilter",
            "validPolygon": "readValidPolygon",
            "apCorrMap": "readApCorrMap",
            "visitInfo": "readVisitInfo",
            "transmissionCurve": "readTransmissionCurve",
            "detector": "readDetector",
            "exposureInfo": "readExposureInfo",
        }
        if (methodName := standardComponents.get(component)) is not None:
            result = getattr(self.reader, methodName)()
            if component == "filter":
                return self._fixFilterLabels(result)
            return result
        # Delegate to MaskedImage and ImageBase implementations for the rest.
        return super().readComponent(component)

    def readFull(self) -> Any:
        # Docstring inherited.
        amplifier, detector, _ = standardizeAmplifierParameters(
            self.checked_parameters,
            self.reader.readDetector(),
        )
        if amplifier is not None:
            amplifier_isolator = AmplifierIsolator(
                amplifier,
                self.reader.readBBox(),
                detector,
            )
            result = amplifier_isolator.transform_subimage(
                self.reader.read(bbox=amplifier_isolator.subimage_bbox, dtype=self.storageClass_dtype)
            )
            result.setDetector(amplifier_isolator.make_detector())
        else:
            result = self.reader.read(**self.checked_parameters, dtype=self.storageClass_dtype)
        result.getInfo().setFilter(self._fixFilterLabels(result.getInfo().getFilter()))
        return result

    def _fixFilterLabels(
        self, file_filter_label: lsst.afw.image.FilterLabel, should_be_standardized: bool | None = None
    ) -> lsst.afw.image.FilterLabel:
        """Compare the filter label read from the file with the one in the
        data ID.

        Parameters
        ----------
        file_filter_label : `lsst.afw.image.FilterLabel` or `None`
            Filter label read from the file, if there was one.
        should_be_standardized : `bool`, optional
            If `True`, expect ``file_filter_label`` to be consistent with the
            data ID and warn only if it is not.  If `False`, expect it to be
            inconsistent and warn only if the data ID is incomplete and hence
            the `FilterLabel` cannot be fixed.  If `None` (default) guess
            whether the file should be standardized by looking at the
            serialization version number in file, which requires this method to
            have been run after `readFull` or `readComponent`.

        Returns
        -------
        filter_label : `lsst.afw.image.FilterLabel` or `None`
            The preferred filter label; may be the given one or one built from
            the data ID.  `None` is returned if there should never be any
            filters associated with this dataset type.

        Notes
        -----
        Most test coverage for this method is in ci_hsc_gen3, where we have
        much easier access to test data that exhibits the problems it attempts
        to solve.
        """
        # Remember filter data ID keys that weren't in this particular data ID,
        # so we can warn about them later.
        missing = []
        band = None
        physical_filter = None
        if "band" in self.data_id.dimensions.names:
            band = self.data_id.get("band")
            # band isn't in the data ID; is that just because this data ID
            # hasn't been filled in with everything the Registry knows, or
            # because this dataset is never associated with a band?
            if band is None and not self.data_id.hasFull() and "band" in self.data_id.dimensions.implied:
                missing.append("band")
        if "physical_filter" in self.data_id.dimensions.names:
            physical_filter = self.data_id.get("physical_filter")
            # Same check as above for band, but for physical_filter.
            if (
                physical_filter is None
                and not self.data_id.hasFull()
                and "physical_filter" in self.data_id.dimensions.implied
            ):
                missing.append("physical_filter")
        if should_be_standardized is None:
            version = self.reader.readSerializationVersion()
            should_be_standardized = version >= 2
        if missing:
            # Data ID identifies a filter but the actual filter label values
            # haven't been fetched from the database; we have no choice but
            # to use the one in the file.
            # Warn if that's more likely than not to be bad, because the file
            # predates filter standardization.
            if not should_be_standardized:
                warnings.warn(
                    f"Data ID {self.data_id} is missing (implied) value(s) for {missing}; "
                    "the correctness of this Exposure's FilterLabel cannot be guaranteed. "
                    "Call Registry.expandDataId before Butler.get to avoid this.",
                    # Report the warning from outside of middleware or the
                    # relevant runQuantum method.
                    stacklevel=find_outside_stacklevel(
                        "lsst.obs.base", "lsst.pipe.base", "lsst.daf.butler", allow_methods={"runQuantum"}
                    ),
                )
            return file_filter_label
        if band is None and physical_filter is None:
            data_id_filter_label = None
        else:
            data_id_filter_label = FilterLabel(band=band, physical=physical_filter)
        if data_id_filter_label != file_filter_label and should_be_standardized:
            # File was written after FilterLabel and standardization, but its
            # FilterLabel doesn't agree with the data ID: this indicates a bug
            # in whatever code produced the Exposure (though it may be one that
            # has been fixed since the file was written).
            warnings.warn(
                f"Reading {self.file_descriptor.location} with data ID {self.data_id}: "
                f"filter label mismatch (file is {file_filter_label}, data ID is "
                f"{data_id_filter_label}).  This is probably a bug in the code that produced it.",
                stacklevel=find_outside_stacklevel(
                    "lsst.obs.base", "lsst.pipe.base", "lsst.daf.butler", allow_methods={"runQuantum"}
                ),
            )
        return data_id_filter_label
