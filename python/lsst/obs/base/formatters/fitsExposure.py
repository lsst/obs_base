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
# Needed for ApCorrMap to resolve properly
from lsst.afw.math import BoundedField  # noqa: F401


class FitsExposureFormatter(Formatter):
    """Interface for reading and writing Exposures to and from FITS files.

    This Formatter supports write recipes.

    Each ``FitsExposureFormatter`` recipe for FITS compression should
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

    A very simple example YAML recipe:

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
    _readerClass = ExposureFitsReader

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
        componentMap = {'wcs': ('readWcs', False),
                        'coaddInputs': ('readCoaddInputs', False),
                        'psf': ('readPsf', False),
                        'image': ('readImage', True),
                        'mask': ('readMask', True),
                        'variance': ('readVariance', True),
                        'photoCalib': ('readPhotoCalib', False),
                        'bbox': ('readBBox', True),
                        'dimensions': ('readBBox', True),
                        'xy0': ('readXY0', True),
                        # TODO: deprecate in DM-27170, remove in DM-27177
                        'filter': ('readFilter', False),
                        # TODO: deprecate in DM-27177, remove in DM-27811
                        'filterLabel': ('readFilterLabel', False),
                        'validPolygon': ('readValidPolygon', False),
                        'apCorrMap': ('readApCorrMap', False),
                        'visitInfo': ('readVisitInfo', False),
                        'transmissionCurve': ('readTransmissionCurve', False),
                        'detector': ('readDetector', False),
                        'exposureInfo': ('readExposureInfo', False),
                        }
        method, hasParams = componentMap.get(component, (None, False))

        if method:
            # This reader can read standalone Image/Mask files as well
            # when dealing with components.
            reader = self._readerClass(self.fileDescriptor.location.path)
            caller = getattr(reader, method, None)

            if caller:
                if parameters is None:
                    parameters = self.fileDescriptor.parameters
                if parameters is None:
                    parameters = {}
                self.fileDescriptor.storageClass.validateParameters(parameters)

                if hasParams and parameters:
                    thisComponent = caller(**parameters)
                else:
                    thisComponent = caller()
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
        reader = self._readerClass(fileDescriptor.location.path)
        return reader.read(**parameters)

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


class FitsImageFormatter(FitsExposureFormatter):
    """Specialisation for `~lsst.afw.image.Image` reading.
    """

    _readerClass = ImageFitsReader


class FitsMaskFormatter(FitsExposureFormatter):
    """Specialisation for `~lsst.afw.image.Mask` reading.
    """

    _readerClass = MaskFitsReader


class FitsMaskedImageFormatter(FitsExposureFormatter):
    """Specialisation for `~lsst.afw.image.MaskedImage` reading.
    """

    _readerClass = MaskedImageFitsReader
