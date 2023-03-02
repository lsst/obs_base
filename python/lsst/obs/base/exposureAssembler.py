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

"""Support for assembling and disassembling afw Exposures."""

import logging
from typing import Any, Dict, Iterable, Mapping, Optional, Set, Tuple, Type

# Need to enable PSFs to be instantiated
import lsst.afw.detection  # noqa: F401
from lsst.afw.image import Exposure, makeExposure, makeMaskedImage
from lsst.daf.butler import DatasetComponent, StorageClassDelegate

log = logging.getLogger(__name__)


class ExposureAssembler(StorageClassDelegate):
    EXPOSURE_COMPONENTS = set(("image", "variance", "mask", "wcs", "psf"))
    EXPOSURE_INFO_COMPONENTS = set(
        (
            "apCorrMap",
            "coaddInputs",
            "photoCalib",
            "metadata",
            "filter",
            "transmissionCurve",
            "visitInfo",
            "detector",
            "validPolygon",
            "summaryStats",
            "id",
        )
    )
    EXPOSURE_READ_COMPONENTS = {
        "bbox",
        "dimensions",
        "xy0",
    }

    COMPONENT_MAP = {"bbox": "BBox", "xy0": "XY0"}
    """Map component name to actual getter name."""

    def _groupRequestedComponents(self) -> Tuple[Set[str], Set[str]]:
        """Group requested components into top level and ExposureInfo.

        Returns
        -------
        expComps : `set` [`str`]
            Components associated with the top level Exposure.
        expInfoComps : `set` [`str`]
            Components associated with the ExposureInfo

        Raises
        ------
        ValueError
            There are components defined in the storage class that are not
            expected by this assembler.
        """
        requested = set(self.storageClass.components.keys())

        # Check that we are requesting something that we support
        unknown = requested - (self.EXPOSURE_COMPONENTS | self.EXPOSURE_INFO_COMPONENTS)
        if unknown:
            raise ValueError("Asking for unrecognized component: {}".format(unknown))

        expItems = requested & self.EXPOSURE_COMPONENTS
        expInfoItems = requested & self.EXPOSURE_INFO_COMPONENTS
        return expItems, expInfoItems

    def getComponent(self, composite: lsst.afw.image.Exposure, componentName: str) -> Any:
        """Get a component from an Exposure

        Parameters
        ----------
        composite : `~lsst.afw.image.Exposure`
            `Exposure` to access component.
        componentName : `str`
            Name of component to retrieve.

        Returns
        -------
        component : `object`
            The component. Can be None.

        Raises
        ------
        AttributeError
            The component can not be found.
        """
        if componentName in self.EXPOSURE_COMPONENTS or componentName in self.EXPOSURE_READ_COMPONENTS:
            # Use getter translation if relevant or the name itself
            return super().getComponent(composite, self.COMPONENT_MAP.get(componentName, componentName))
        elif componentName in self.EXPOSURE_INFO_COMPONENTS:
            if hasattr(composite, "getInfo"):
                # it is possible for this method to be called with
                # an ExposureInfo composite so trap for that and only get
                # the ExposureInfo if the method is supported
                composite = composite.getInfo()
            return super().getComponent(composite, self.COMPONENT_MAP.get(componentName, componentName))
        else:
            raise AttributeError(
                "Do not know how to retrieve component {} from {}".format(componentName, type(composite))
            )

    def getValidComponents(self, composite: Exposure) -> Dict[str, Any]:
        """Extract all non-None components from a composite.

        Parameters
        ----------
        composite : `object`
            Composite from which to extract components.

        Returns
        -------
        comps : `dict`
            Non-None components extracted from the composite, indexed by the
            component name as derived from the `self.storageClass`.
        """
        # For Exposure we call the generic version twice: once for top level
        # components, and again for ExposureInfo.
        expItems, expInfoItems = self._groupRequestedComponents()

        components = super().getValidComponents(composite)
        infoComps = super().getValidComponents(composite.getInfo())
        components.update(infoComps)
        return components

    def disassemble(
        self, composite: Any, subset: Optional[Iterable] = None, override: Optional[Any] = None
    ) -> Dict[str, DatasetComponent]:
        """Disassemble an afw Exposure.

        This implementation attempts to extract components from the parent
        by looking for attributes of the same name or getter methods derived
        from the component name.

        Parameters
        ----------
        composite : `~lsst.afw.image.Exposure`
            `Exposure` composite object consisting of components to be
            extracted.
        subset : iterable, optional
            Not supported by this assembler.
        override : `object`, optional
            Not supported by this assembler.

        Returns
        -------
        components : `dict`
            `dict` with keys matching the components defined in
            `self.storageClass` and values being `DatasetComponent` instances
            describing the component.

        Raises
        ------
        ValueError
            A requested component can not be found in the parent using generic
            lookups.
        TypeError
            The parent object does not match the supplied `self.storageClass`.

        Notes
        -----
        If a PSF is present but is not persistable, the PSF will not be
        included in the returned components.
        """
        if subset is not None:
            raise NotImplementedError(
                "ExposureAssembler does not support the 'subset' argument to disassemble."
            )
        if override is not None:
            raise NotImplementedError(
                "ExposureAssembler does not support the 'override' argument to disassemble."
            )
        if not self.storageClass.validateInstance(composite):
            raise TypeError(
                "Unexpected type mismatch between parent and StorageClass"
                " ({} != {})".format(type(composite), self.storageClass.pytype)
            )

        # Only look for components that are defined by the StorageClass
        components: Dict[str, DatasetComponent] = {}
        expItems, expInfoItems = self._groupRequestedComponents()

        fromExposure = super().disassemble(composite, subset=expItems)
        assert fromExposure is not None, "Base class implementation guarantees this, but ABC does not."
        components.update(fromExposure)

        fromExposureInfo = super().disassemble(composite, subset=expInfoItems, override=composite.getInfo())
        assert fromExposureInfo is not None, "Base class implementation guarantees this, but ABC does not."
        components.update(fromExposureInfo)

        if "psf" in components and not components["psf"].component.isPersistable():
            log.warning(
                "PSF of type %s is not persistable and has been ignored.",
                type(components["psf"].component).__name__,
            )
            del components["psf"]

        return components

    def assemble(self, components: Dict[str, Any], pytype: Optional[Type] = None) -> Exposure:
        """Construct an Exposure from components.

        Parameters
        ----------
        components : `dict`
            All the components from which to construct the Exposure.
            Some can be missing.
        pytype : `type`, optional
            Not supported by this assembler.

        Returns
        -------
        exposure : `~lsst.afw.image.Exposure`
            Assembled exposure.

        Raises
        ------
        ValueError
            Some supplied components are not recognized.
        """
        if pytype is not None:
            raise NotImplementedError("ExposureAssembler does not support the 'pytype' argument to assemble.")
        components = components.copy()
        maskedImageComponents = {}
        hasMaskedImage = False
        for component in ("image", "variance", "mask"):
            value = None
            if component in components:
                hasMaskedImage = True
                value = components.pop(component)
            maskedImageComponents[component] = value

        wcs = None
        if "wcs" in components:
            wcs = components.pop("wcs")

        pytype = self.storageClass.pytype
        if hasMaskedImage:
            maskedImage = makeMaskedImage(**maskedImageComponents)
            exposure = makeExposure(maskedImage, wcs=wcs)

            if not isinstance(exposure, pytype):
                raise RuntimeError(
                    "Unexpected type created in assembly;"
                    " was {} expected {}".format(type(exposure), pytype)
                )

        else:
            exposure = pytype()
            if wcs is not None:
                exposure.setWcs(wcs)

        # Set other components
        exposure.setPsf(components.pop("psf", None))
        exposure.setPhotoCalib(components.pop("photoCalib", None))

        info = exposure.getInfo()
        if "visitInfo" in components:
            info.setVisitInfo(components.pop("visitInfo"))
        # Until DM-32138, "visitInfo" and "id" can both set the exposure ID.
        # While they should always be consistent unless a component is
        # corrupted, handle "id" second to ensure it takes precedence.
        if "id" in components:
            info.id = components.pop("id")
        info.setApCorrMap(components.pop("apCorrMap", None))
        info.setCoaddInputs(components.pop("coaddInputs", None))
        info.setMetadata(components.pop("metadata", None))
        info.setValidPolygon(components.pop("validPolygon", None))
        info.setDetector(components.pop("detector", None))
        info.setTransmissionCurve(components.pop("transmissionCurve", None))
        info.setSummaryStats(components.pop("summaryStats", None))

        info.setFilter(components.pop("filter", None))

        # If we have some components left over that is a problem
        if components:
            raise ValueError(
                "The following components were not understood: {}".format(list(components.keys()))
            )

        return exposure

    def handleParameters(self, inMemoryDataset: Any, parameters: Optional[Mapping[str, Any]] = None) -> Any:
        """Modify the in-memory dataset using the supplied parameters,
        returning a possibly new object.

        Parameters
        ----------
        inMemoryDataset : `object`
            Object to modify based on the parameters.
        parameters : `dict`, optional
            Parameters to apply. Values are specific to the parameter.
            Supported parameters are defined in the associated
            `StorageClass`.  If no relevant parameters are specified the
            inMemoryDataset will be return unchanged.

        Returns
        -------
        inMemoryDataset : `object`
            Updated form of supplied in-memory dataset, after parameters
            have been used.
        """
        if parameters is None:
            return inMemoryDataset
        # Understood by *this* subset command
        understood = ("bbox", "origin")
        use = self.storageClass.filterParameters(parameters, subset=understood)
        if use:
            inMemoryDataset = inMemoryDataset.subset(**use)

        return inMemoryDataset

    @classmethod
    def selectResponsibleComponent(cls, readComponent: str, fromComponents: Set[Optional[str]]) -> str:
        # Docstring inherited.
        imageComponents = ["mask", "image", "variance"]
        forwarderMap = {
            "bbox": imageComponents,
            "dimensions": imageComponents,
            "xy0": imageComponents,
        }
        forwarder = forwarderMap.get(readComponent)
        if forwarder is not None:
            for c in forwarder:
                if c in fromComponents:
                    return c
        raise ValueError(f"Can not calculate read component {readComponent} from {fromComponents}")
