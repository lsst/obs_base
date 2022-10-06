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
    "FitsMaskFormatter",
    "FitsMaskedImageFormatter",
    "standardizeAmplifierParameters",
)

import warnings
from abc import abstractmethod
from typing import AbstractSet, ClassVar

from lsst.afw.cameraGeom import AmplifierGeometryComparison, AmplifierIsolator
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
from lsst.daf.base import PropertySet
from lsst.daf.butler import Formatter
from lsst.utils.classes import cached_getter


class FitsImageFormatterBase(Formatter):
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

    extension = ".fits"
    supportedExtensions: ClassVar[AbstractSet[str]] = frozenset(
        {".fits", ".fits.gz", ".fits.fz", ".fz", ".fit"}
    )

    unsupportedParameters: ClassVar[AbstractSet[str]] = frozenset()
    """Support all parameters."""

    @property
    @cached_getter
    def checked_parameters(self):
        """The parameters passed by the butler user, after checking them
        against the storage class and transforming `None` into an empty `dict`
        (`dict`).

        This is computed on first use and then cached.  It should never be
        accessed when writing.  Subclasses that need additional checking should
        delegate to `super` and then check the result before returning it.
        """
        parameters = self.fileDescriptor.parameters
        if parameters is None:
            parameters = {}
        self.fileDescriptor.storageClass.validateParameters(parameters)
        return parameters

    def read(self, component=None):
        # Docstring inherited.
        if component is not None:
            return self.readComponent(component)
        return self.readFull()

    @abstractmethod
    def readComponent(self, component):
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
    def readFull(self):
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
    contain ``compression`` and ``scaling`` entries. Defaults will be
    provided for any missing elements under ``compression`` and
    ``scaling``.

    The allowed entries under ``compression`` are:

    * ``algorithm`` (`str`): compression algorithm to use
    * ``rows`` (`int`): number of rows per tile (0 = entire dimension)
    * ``columns`` (`int`): number of columns per tile (0 = entire dimension)
    * ``quantizeLevel`` (`float`): cfitsio quantization level

    The allowed entries under ``scaling`` are:

    * ``algorithm`` (`str`): scaling algorithm to use
    * ``bitpix`` (`int`): bits per pixel (0,8,16,32,64,-32,-64)
    * ``fuzz`` (`bool`): fuzz the values when quantising floating-point values?
    * ``seed`` (`int`): seed for random number generator when fuzzing
    * ``maskPlanes`` (`list` of `str`): mask planes to ignore when doing
      statistics
    * ``quantizeLevel`` (`float`): divisor of the standard deviation for
      ``STDEV_*`` scaling
    * ``quantizePad`` (`float`): number of stdev to allow on the low side (for
      ``STDEV_POSITIVE``/``NEGATIVE``)
    * ``bscale`` (`float`): manually specified ``BSCALE``
      (for ``MANUAL`` scaling)
    * ``bzero`` (`float`): manually specified ``BSCALE``
      (for ``MANUAL`` scaling)

    A very simple example YAML recipe (for the ``Exposure`` specialization):

    .. code-block:: yaml

        lsst.obs.base.fitsExposureFormatter.FitsExposureFormatter:
          default:
            image: &default
              compression:
                algorithm: GZIP_SHUFFLE
            mask: *default
            variance: *default

    """

    supportedWriteParameters = frozenset({"recipe"})
    ReaderClass: type  # must be set by concrete subclasses

    @property
    @cached_getter
    def reader(self):
        """The reader object that backs this formatter's read operations.

        This is computed on first use and then cached.  It should never be
        accessed when writing.
        """
        return self.ReaderClass(self.fileDescriptor.location.path)

    def readComponent(self, component):
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

    def readFull(self):
        # Docstring inherited.
        return self.reader.read(**self.checked_parameters)

    def write(self, inMemoryDataset):
        """Write a Python object to a file.

        Parameters
        ----------
        inMemoryDataset : `object`
            The Python object to store.
        """
        # Update the location with the formatter-preferred file extension
        self.fileDescriptor.location.updateExtension(self.extension)
        outputPath = self.fileDescriptor.location.path

        # check to see if we have a recipe requested
        recipeName = self.writeParameters.get("recipe")
        recipe = self.getImageCompressionSettings(recipeName)
        if recipe:
            # Can not construct a PropertySet from a hierarchical
            # dict but can update one.
            ps = PropertySet()
            ps.update(recipe)
            inMemoryDataset.writeFitsWithOptions(outputPath, options=ps)
        else:
            inMemoryDataset.writeFits(outputPath)

    def getImageCompressionSettings(self, recipeName):
        """Retrieve the relevant compression settings for this recipe.

        Parameters
        ----------
        recipeName : `str`
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
            if "default" not in self.writeRecipes:
                return {}
            recipeName = "default"

        if recipeName not in self.writeRecipes:
            raise RuntimeError(f"Unrecognized recipe option given for compression: {recipeName}")

        recipe = self.writeRecipes[recipeName]

        # Set the seed based on dataId
        seed = hash(tuple(self.dataId.items())) % 2**31
        for plane in ("image", "mask", "variance"):
            if plane in recipe and "scaling" in recipe[plane]:
                scaling = recipe[plane]["scaling"]
                if "seed" in scaling and scaling["seed"] == 0:
                    scaling["seed"] = seed

        return recipe

    @classmethod
    def validateWriteRecipes(cls, recipes):
        """Validate supplied recipes for this formatter.

        The recipes are supplemented with default values where appropriate.

        TODO: replace this custom validation code with Cerberus (DM-11846)

        Parameters
        ----------
        recipes : `dict`
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
        # Schemas define what should be there, and the default values (and by
        # the default value, the expected type).
        compressionSchema = {
            "algorithm": "NONE",
            "rows": 1,
            "columns": 0,
            "quantizeLevel": 0.0,
        }
        scalingSchema = {
            "algorithm": "NONE",
            "bitpix": 0,
            "maskPlanes": ["NO_DATA"],
            "seed": 0,
            "quantizeLevel": 4.0,
            "quantizePad": 5.0,
            "fuzz": True,
            "bscale": 1.0,
            "bzero": 0.0,
        }

        if not recipes:
            # We can not insist on recipes being specified
            return recipes

        def checkUnrecognized(entry, allowed, description):
            """Check to see if the entry contains unrecognised keywords"""
            unrecognized = set(entry) - set(allowed)
            if unrecognized:
                raise RuntimeError(
                    f"Unrecognized entries when parsing image compression recipe {description}: "
                    f"{unrecognized}"
                )

        validated = {}
        for name in recipes:
            checkUnrecognized(recipes[name], ["image", "mask", "variance"], name)
            validated[name] = {}
            for plane in ("image", "mask", "variance"):
                checkUnrecognized(recipes[name][plane], ["compression", "scaling"], f"{name}->{plane}")

                np = {}
                validated[name][plane] = np
                for settings, schema in (("compression", compressionSchema), ("scaling", scalingSchema)):
                    np[settings] = {}
                    if settings not in recipes[name][plane]:
                        for key in schema:
                            np[settings][key] = schema[key]
                        continue
                    entry = recipes[name][plane][settings]
                    checkUnrecognized(entry, schema.keys(), f"{name}->{plane}->{settings}")
                    for key in schema:
                        value = type(schema[key])(entry[key]) if key in entry else schema[key]
                        np[settings][key] = value
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

    def readComponent(self, component):
        # Docstring inherited.
        if component == "image":
            return self.reader.readImage(**self.checked_parameters)
        elif component == "mask":
            return self.reader.readMask(**self.checked_parameters)
        elif component == "variance":
            return self.reader.readVariance(**self.checked_parameters)
        else:
            # Delegate to base for bbox, dimensions, xy0.
            return super().readComponent(component)


def standardizeAmplifierParameters(parameters, on_disk_detector):
    """Preprocess the Exposure storage class's "amp" and "detector" parameters

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
    if isinstance(amplifier, (int, str)):
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

    ReaderClass = ExposureFitsReader

    def readComponent(self, component):
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

    def readFull(self):
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
                self.reader.read(bbox=amplifier_isolator.subimage_bbox)
            )
            result.setDetector(amplifier_isolator.make_detector())
        else:
            result = self.reader.read(**self.checked_parameters)
        result.getInfo().setFilter(self._fixFilterLabels(result.getInfo().getFilter()))
        return result

    def _fixFilterLabels(self, file_filter_label, should_be_standardized=None):
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
        if "band" in self.dataId.graph.dimensions.names:
            band = self.dataId.get("band")
            # band isn't in the data ID; is that just because this data ID
            # hasn't been filled in with everything the Registry knows, or
            # because this dataset is never associated with a band?
            if band is None and not self.dataId.hasFull() and "band" in self.dataId.graph.implied.names:
                missing.append("band")
        if "physical_filter" in self.dataId.graph.dimensions.names:
            physical_filter = self.dataId.get("physical_filter")
            # Same check as above for band, but for physical_filter.
            if (
                physical_filter is None
                and not self.dataId.hasFull()
                and "physical_filter" in self.dataId.graph.implied.names
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
                    f"Data ID {self.dataId} is missing (implied) value(s) for {missing}; "
                    "the correctness of this Exposure's FilterLabel cannot be guaranteed. "
                    "Call Registry.expandDataId before Butler.get to avoid this."
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
                f"Reading {self.fileDescriptor.location} with data ID {self.dataId}: "
                f"filter label mismatch (file is {file_filter_label}, data ID is "
                f"{data_id_filter_label}).  This is probably a bug in the code that produced it."
            )
        return data_id_filter_label
