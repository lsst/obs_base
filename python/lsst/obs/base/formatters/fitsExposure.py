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

__all__ = ("FitsExposureFormatter", "FitsImageFormatter", "FitsMaskFormatter",
           "FitsMaskedImageFormatter")

from astro_metadata_translator import fix_header
from lsst.daf.base import PropertySet
from lsst.daf.butler import Formatter
# Do not use ExposureFitsReader.readMetadata because that strips
# out lots of headers and there is no way to recover them
from lsst.afw.fits import readMetadata
from lsst.afw.image import ExposureFitsReader, ImageFitsReader, MaskFitsReader, MaskedImageFitsReader
from lsst.afw.image import ExposureInfo
# Needed for ApCorrMap to resolve properly
from lsst.afw.math import BoundedField  # noqa: F401

from ..exposureAssembler import fixFilterLabels


class FitsImageFormatterBase(Formatter):
    """Base class interface for reading and writing afw images to and from
    FITS files.

    This Formatter supports write recipes.

    Each ``FitsImageFormatterBase`` recipe for FITS compression should
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
    supportedExtensions = frozenset({".fits", ".fits.gz", ".fits.fz", ".fz", ".fit"})
    extension = ".fits"
    _metadata = None
    supportedWriteParameters = frozenset({"recipe"})
    _readerClass: type  # must be set by concrete subclasses

    unsupportedParameters = {}
    """Support all parameters."""

    @property
    def metadata(self):
        """The metadata read from this file. It will be stripped as
        components are extracted from it
        (`lsst.daf.base.PropertyList`).
        """
        if self._metadata is None:
            self._metadata = self.readMetadata()
        return self._metadata

    def readMetadata(self):
        """Read all header metadata directly into a PropertyList.

        Returns
        -------
        metadata : `~lsst.daf.base.PropertyList`
            Header metadata.
        """
        md = readMetadata(self.fileDescriptor.location.path)
        fix_header(md)
        return md

    def stripMetadata(self):
        """Remove metadata entries that are parsed into components.

        This is only called when just the metadata is requested; stripping
        entries there forces code that wants other components to ask for those
        components directly rather than trying to extract them from the
        metadata manually, which is fragile.  This behavior is an intentional
        change from Gen2.

        Parameters
        ----------
        metadata : `~lsst.daf.base.PropertyList`
            Header metadata, to be modified in-place.
        """
        # TODO: make sure this covers everything, by delegating to something
        # that doesn't yet exist in afw.image.ExposureInfo.
        from lsst.afw.image import bboxFromMetadata
        from lsst.afw.geom import makeSkyWcs

        # Protect against the metadata being missing
        try:
            bboxFromMetadata(self.metadata)  # always strips
        except LookupError:
            pass
        try:
            makeSkyWcs(self.metadata, strip=True)
        except Exception:
            pass

    def readComponent(self, component, parameters=None):
        """Read a component held by the Exposure.

        Parameters
        ----------
        component : `str`, optional
            Component to read from the file.
        parameters : `dict`, optional
            If specified, a dictionary of slicing parameters that
            overrides those in ``fileDescriptor``.

        Returns
        -------
        obj : component-dependent
            In-memory component object.

        Raises
        ------
        KeyError
            Raised if the requested component cannot be handled.
        """

        # Metadata is handled explicitly elsewhere
        componentMap = {'wcs': ('readWcs', False, None),
                        'coaddInputs': ('readCoaddInputs', False, None),
                        'psf': ('readPsf', False, None),
                        'image': ('readImage', True, None),
                        'mask': ('readMask', True, None),
                        'variance': ('readVariance', True, None),
                        'photoCalib': ('readPhotoCalib', False, None),
                        'bbox': ('readBBox', True, None),
                        'dimensions': ('readBBox', True, None),
                        'xy0': ('readXY0', True, None),
                        # TODO: deprecate in DM-27170, remove in DM-27177
                        'filter': ('readFilter', False, None),
                        # TODO: deprecate in DM-27177, remove in DM-27811
                        'filterLabel': ('readFilterLabel', False, None),
                        'validPolygon': ('readValidPolygon', False, None),
                        'apCorrMap': ('readApCorrMap', False, None),
                        'visitInfo': ('readVisitInfo', False, None),
                        'transmissionCurve': ('readTransmissionCurve', False, None),
                        'detector': ('readDetector', False, None),
                        'exposureInfo': ('readExposureInfo', False, None),
                        'summaryStats': ('readComponent', False, ExposureInfo.KEY_SUMMARY_STATS),
                        }
        method, hasParams, componentName = componentMap.get(component, (None, False, None))

        if method:
            # This reader can read standalone Image/Mask files as well
            # when dealing with components.
            self._reader = self._readerClass(self.fileDescriptor.location.path)
            caller = getattr(self._reader, method, None)

            if caller:
                if parameters is None:
                    parameters = self.fileDescriptor.parameters
                if parameters is None:
                    parameters = {}
                self.fileDescriptor.storageClass.validateParameters(parameters)

                if componentName is None:
                    if hasParams and parameters:
                        thisComponent = caller(**parameters)
                    else:
                        thisComponent = caller()
                else:
                    thisComponent = caller(componentName)

                if component == "dimensions" and thisComponent is not None:
                    thisComponent = thisComponent.getDimensions()

                return thisComponent
        else:
            raise KeyError(f"Unknown component requested: {component}")

    def readFull(self, parameters=None):
        """Read the full Exposure object.

        Parameters
        ----------
        parameters : `dict`, optional
            If specified a dictionary of slicing parameters that overrides
            those in ``fileDescriptor``.

        Returns
        -------
        exposure : `~lsst.afw.image.Exposure`
            Complete in-memory exposure.
        """
        fileDescriptor = self.fileDescriptor
        if parameters is None:
            parameters = fileDescriptor.parameters
        if parameters is None:
            parameters = {}
        fileDescriptor.storageClass.validateParameters(parameters)
        self._reader = self._readerClass(fileDescriptor.location.path)
        return self._reader.read(**parameters)

    def read(self, component=None):
        """Read data from a file.

        Parameters
        ----------
        component : `str`, optional
            Component to read from the file. Only used if the `StorageClass`
            for reading differed from the `StorageClass` used to write the
            file.

        Returns
        -------
        inMemoryDataset : `object`
            The requested data as a Python object. The type of object
            is controlled by the specific formatter.

        Raises
        ------
        ValueError
            Component requested but this file does not seem to be a concrete
            composite.
        KeyError
            Raised when parameters passed with fileDescriptor are not
            supported.
        """
        fileDescriptor = self.fileDescriptor
        if fileDescriptor.readStorageClass != fileDescriptor.storageClass:
            if component == "metadata":
                self.stripMetadata()
                return self.metadata
            elif component is not None:
                return self.readComponent(component)
            else:
                raise ValueError("Storage class inconsistency ({} vs {}) but no"
                                 " component requested".format(fileDescriptor.readStorageClass.name,
                                                               fileDescriptor.storageClass.name))
        return self.readFull()

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
                    f"{unrecognized}")

        validated = {}
        for name in recipes:
            checkUnrecognized(recipes[name], ["image", "mask", "variance"], name)
            validated[name] = {}
            for plane in ("image", "mask", "variance"):
                checkUnrecognized(recipes[name][plane], ["compression", "scaling"],
                                  f"{name}->{plane}")

                np = {}
                validated[name][plane] = np
                for settings, schema in (("compression", compressionSchema),
                                         ("scaling", scalingSchema)):
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


class FitsExposureFormatter(FitsImageFormatterBase):
    """Specialization for `~lsst.afw.image.Exposure` reading.
    """

    _readerClass = ExposureFitsReader

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
        if should_be_standardized is None:
            version = self._reader.readSerializationVersion()
            should_be_standardized = (version >= 2)

        return fixFilterLabels(file_filter_label, self.dataId, should_be_standardized=should_be_standardized,
                               msg=f"Reading file {self.fileDescriptor.location}")

    def readComponent(self, component, parameters=None):
        # Docstring inherited.
        obj = super().readComponent(component, parameters)
        if component == "filterLabel":
            return self._fixFilterLabels(obj)
        else:
            return obj

    def readFull(self, parameters=None):
        # Docstring inherited.
        full = super().readFull(parameters)
        full.getInfo().setFilterLabel(self._fixFilterLabels(full.getInfo().getFilterLabel()))
        return full


class FitsImageFormatter(FitsImageFormatterBase):
    """Specialisation for `~lsst.afw.image.Image` reading.
    """

    _readerClass = ImageFitsReader


class FitsMaskFormatter(FitsImageFormatterBase):
    """Specialisation for `~lsst.afw.image.Mask` reading.
    """

    _readerClass = MaskFitsReader


class FitsMaskedImageFormatter(FitsImageFormatterBase):
    """Specialisation for `~lsst.afw.image.MaskedImage` reading.
    """

    _readerClass = MaskedImageFitsReader
