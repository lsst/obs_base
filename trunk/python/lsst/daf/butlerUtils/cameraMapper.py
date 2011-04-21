#!/bin/env python
# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import os
import re

import lsst.daf.persistence as dafPersist
from lsst.daf.butlerUtils import ExposureMapping, CalibrationMapping, DatasetMapping, Registry
import lsst.daf.base as dafBase
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.cameraGeom as afwCameraGeom
import lsst.afw.cameraGeom.utils as cameraGeomUtils
import lsst.afw.image.utils as imageUtils
import lsst.pex.logging as pexLog
import lsst.pex.policy as pexPolicy

"""This module defines the CameraMapper base class."""

class CameraMapper(dafPersist.Mapper):
    
    """CameraMapper is a base class for mappers that handle images from a
    camera and products derived from them.  This provides an abstraction layer
    between the data on disk and the code.

    Public methods: keys, queryMetadata, getDatasetTypes, map,
    canStandardize, standardize

    Mappers for specific data sources (e.g., CFHT Megacam, LSST
    simulations, etc.) should inherit this class.

    A camera is assumed to consist of one or more rafts, each composed of
    multiple CCDs.  Each CCD is in turn composed of one or more amplifiers
    (amps).  A camera is also assumed to have a camera geometry description
    (CameraGeom object) as a policy file, a filter description (Filter class
    static configuration) as another policy file, and an optional defects
    description directory.

    Information from the camera geometry and defects are inserted into all
    Exposure objects returned.

    The mapper uses one or two registries to retrieve metadata about the
    images.  The first is a registry of all raw exposures.  This must contain
    the time of the observation.  One or more tables (or the equivalent)
    within the registry are used to look up data identifier components that
    are not specified by the user (e.g. filter) and to return results for
    metadata queries.  The second is an optional registry of all calibration
    data.  This should contain validity start and end entries for each
    calibration dataset in the same timescale as the observation time.

    The following method must be provided by the subclass:

    _extractDetectorName(self, dataId): returns the detector name for a CCD
    (e.g., "CFHT 21", "R:1,2 S:3,4") as used in the AFW CameraGeom class given
    a dataset identifier referring to that CCD or a subcomponent of it.

    Other methods that the subclass may wish to override include:

    _transformId(self, dataId): transformation of a data identifier
    from colloquial usage (e.g., "ccdname") to proper/actual usage
    (e.g., "ccd").  The default implementation does nothing.  Note that this
    method should not modify its input parameter.

    _mapActualToPath(self, template, actualId): convert a template path to an
    actual path, using the actual dataset identifier.

    _extractAmpId(self, dataId): extract the amplifier identifer from a
    dataset identifier.  The amplifier identifier has three parts: the
    detector name for the CCD containing the amplifier and x and y coordinates
    for the amplifier within the CCD.  The default implementation just uses
    the amplifier number, which is assumed to be the same as the x coordinate.

    The mapper's behaviors are largely specified by the policy file,
    which consists of:

    calibRoot (string, optional): Path to calibration files; must be provided
    to the mapper if none is specified and calibration files exist

    registryPath (string, optional): Path to the image registry -- absolute,
    relative to the root path provided to the mapper, or relative to the
    current directory if no root path is provided; must be provided to the
    mapper if none is specified

    camera (string): Path to camera geometry policy file from subclassing
    module's repository

    defects (string): Path to defects directory from subclassing module's
    repository

    filters (string): Path to filters policy file from subclassing module's
    repository

    exposures (policy): Exposure mappings (e.g., "raw", "postISR")

    needCalibRegistry (bool): True if a calibration registry is used

    If needCalibRegistry is true, the following is used:

    calibRegistryPath (string, optional): Path to calibration registry --
    absolute or relative to the calibration root path; must be provided to the
    mapper if none is specified

    calibrations (policy): Calibration mappings (e.g., "bias", "flat")

    datasets (policy): Other, non-Exposure mappings

    The 'exposures' and 'calibrations' policies consist of mappings
    (see Mappings class).

    Functions to map (provide a path to the data given a dataset
    identifier dictionary) and standardize (convert data into some standard
    format or type) may be provided in the subclass as "map_{dataset type}"
    and "std_{dataset type}", respectively.

    If non-Exposure datasets cannot be retrieved using standard
    daf_persistence methods alone, a "bypass_{dataset type}" function may be
    provided in the subclass to return the dataset instead of using the
    "datasets" subpolicy.

    Implementations of map_camera and std_camera that should typically be
    sufficient are provided in this base class.
    """

    def __init__(self, policy, repositoryDir,
                 root=None, registry=None, calibRoot=None, calibRegistry=None,
                 provided=None):
        """Initialize the CameraMapper.
        @param policy        (pexPolicy.Policy) Policy with per-camera defaults
                             already merged
        @param repositoryDir (string) Policy repository for the subclassing
                             module (obtained with getRepositoryPath() on the
                             per-camera default dictionary)
        @param root          (string) Root directory for data
        @param registry      (string) Path to registry with data's metadata
        @param calibRoot     (string) Root directory for calibrations
        @param calibRegistry (string) Path to registry with calibrations'
                             metadata
        @param provided      (list of strings) Keys provided by the mapper
        """

        dafPersist.Mapper.__init__(self)

        self.log = pexLog.Log(pexLog.getDefaultLog(), "CameraMapper")

        # Dictionary
        dictFile = pexPolicy.DefaultPolicyFile("daf_butlerUtils",
                "MapperDictionary.paf", "policy")
        dictPolicy = pexPolicy.Policy.createPolicy(dictFile,
                dictFile.getRepositoryPath())
        policy.mergeDefaults(dictPolicy)

        # Root directories
        if root is None:
            root = "."
        if calibRoot is None:
            if policy.exists('calibRoot'):
                calibRoot = policy.getString('calibRoot')
            else:
                calibRoot = root
        # Do any location substitutions
        root = dafPersist.LogicalLocation(root).locString()
        calibRoot = dafPersist.LogicalLocation(calibRoot).locString()
        if not os.path.exists(root):
            self.log.log(pexLog.Log.WARN,
                    "Root directory not found: %s" % (root,))
        if not os.path.exists(calibRoot):
            self.log.log(pexLog.Log.WARN,
                    "Calibration root directory not found: %s" % (calibRoot,))

        # Registries
        self.registry = self._setupRegistry(
                "registry", registry, policy, "registryPath", root)
        if policy.exists('needCalibRegistry') and \
                policy.getBool('needCalibRegistry'):
            calibRegistry = self._setupRegistry(
                    "calibRegistry", calibRegistry,
                    policy, "calibRegistryPath", calibRoot)
        else:
            calibRegistry = None

        # Sub-dictionaries (for exposure/calibration/dataset types)
        expMappingFile = pexPolicy.DefaultPolicyFile("daf_butlerUtils",
                "ExposureMappingDictionary.paf", "policy")
        expMappingPolicy = pexPolicy.Policy.createPolicy(expMappingFile,
                expMappingFile.getRepositoryPath())
        calMappingFile = pexPolicy.DefaultPolicyFile("daf_butlerUtils",
                "CalibrationMappingDictionary.paf", "policy")
        calMappingPolicy = pexPolicy.Policy.createPolicy(calMappingFile,
                calMappingFile.getRepositoryPath())
        dsMappingFile = pexPolicy.DefaultPolicyFile("daf_butlerUtils",
                "DatasetMappingDictionary.paf", "policy")
        dsMappingPolicy = pexPolicy.Policy.createPolicy(dsMappingFile,
                dsMappingFile.getRepositoryPath())

        # Set of valid keys
        self.keySet = set()

        # Mappings
        mappingList = (
                ("exposures", expMappingPolicy, ExposureMapping),
                ("calibrations", calMappingPolicy, CalibrationMapping),
                ("datasets", dsMappingPolicy, DatasetMapping)
                )
        for name, defPolicy, cls in mappingList:
            if policy.exists(name):
                datasets = policy.getPolicy(name)
                mappings = dict()
                setattr(self, name, mappings)
                for datasetType in datasets.names(True):
                    subPolicy = datasets.getPolicy(datasetType)
                    subPolicy.mergeDefaults(defPolicy)
                    if name == "calibrations":
                        mapping = cls(datasetType, subPolicy,
                                self.registry, calibRegistry, calibRoot, provided=provided)
                    else:
                        mapping = cls(datasetType, subPolicy,
                                self.registry, root, provided=provided)
                    self.keySet.update(mapping.keys())
                    mappings[datasetType] = mapping
                    if not hasattr(self, "map_" + datasetType):
                        def mapClosure(dataId, mapper=self, mapping=mapping):
                            return mapping.map(mapper, dataId)
                        setattr(self, "map_" + datasetType, mapClosure)
                        def queryClosure(key, format, dataId, mapping=mapping):
                            return mapping.lookup(format, dataId)
                        setattr(self, "query_" + datasetType, queryClosure)
                    if hasattr(mapping, "standardize") and \
                            not hasattr(self, "std_" + datasetType):
                        def stdClosure(item, dataId,
                                mapper=self, mapping=mapping):
                            return mapping.standardize(mapper, item, dataId)
                        setattr(self, "std_" + datasetType, stdClosure)

                    mapFunc = "map_" + datasetType + "_filename"
                    bypassFunc = "bypass_" + datasetType + "_filename"
                    if not hasattr(self, mapFunc):
                        setattr(self, mapFunc, getattr(self, "map_" + datasetType))
                    if not hasattr(self, bypassFunc):
                        setattr(self, bypassFunc,
                                lambda datasetType, pythonType, location, dataId: location.getLocations())

                    # Set up metadata versions
                    if name == "exposures":
                        expFunc = "map_" + datasetType # Function name to map exposure
                        mdFunc = expFunc + "_md"       # Function name to map metadata
                        bypassFunc = "bypass_" + datasetType + "_md" # Function name to bypass daf_persistence
                        if not hasattr(self, mdFunc):
                            setattr(self, mdFunc, getattr(self, expFunc))
                        if not hasattr(self, bypassFunc):
                            setattr(self, bypassFunc,
                                    lambda datasetType, pythonType, location, dataId:
                                    afwImage.readMetadata(location.getLocations()[0]))
                        if not hasattr(self, "query_" + datasetType + "_md"):
                            setattr(self, "query_" + datasetType + "_md",
                                    getattr(self, "query_" + datasetType))

                        subFunc = expFunc + "_sub" # Function name to map subimage
                        if not hasattr(self, subFunc):
                            def mapSubClosure(dataId, mapper=self, mapping=mapping):
                                subId = dataId.copy()
                                del subId['bbox']
                                loc = mapping.map(mapper, subId)
                                bbox = dataId['bbox']
                                llcX = bbox.getMinX()
                                llcY = bbox.getMinY()
                                width = bbox.getWidth()
                                height = bbox.getHeight()
                                loc.additionalData.set('llcX', llcX)
                                loc.additionalData.set('llcY', llcY)
                                loc.additionalData.set('width', width)
                                loc.additionalData.set('height', height)
                                return loc
                            setattr(self, subFunc, mapSubClosure)
                        if not hasattr(self, "query_" + datasetType + "_sub"):
                            def querySubClosure(key, format, dataId, mapping=mapping):
                                subId = dataId.copy()
                                del subId['bbox']
                                return mapping.lookup(format, subId)
                            setattr(self, "query_" + datasetType + "_sub", querySubClosure)

        # Camera geometry
        self.cameraPolicyLocation = None
        self.camera = None
        if policy.exists('camera'):
            cameraPolicyLocation = policy.getString('camera')

            # must be explicit for ButlerLocation later
            self.cameraPolicyLocation = os.path.join(repositoryDir, cameraPolicyLocation)
            cameraPolicy = pexPolicy.Policy.createPolicy(self.cameraPolicyLocation)
            cameraPolicy = cameraGeomUtils.getGeomPolicy(cameraPolicy)
            self.camera = cameraGeomUtils.makeCamera(cameraPolicy)

        # Defect registry and root
        self.defectRegistry = None
        if policy.exists('defects'):
            self.defectPath = os.path.join(
                    repositoryDir, policy.getString('defects'))
            defectRegistryLocation = os.path.join(
                    self.defectPath, "defectRegistry.sqlite3")
            self.defectRegistry = \
                    Registry.create(defectRegistryLocation)

        # Filters
        if policy.exists('filters'):
            filterPolicyLocation = os.path.join(repositoryDir, policy.getString('filters'))
            filterPolicy = pexPolicy.Policy.createPolicy(filterPolicyLocation)
            imageUtils.defineFiltersFromPolicy(filterPolicy, reset=True)

        # Filter translation table
        self.filters = None

        # Skytile policy
        self.skypolicy = policy.getPolicy("skytiles")

    def keys(self):
        """Return supported keys.
        @return (iterable) List of keys usable in a dataset identifier"""
        return self.keySet

    def map_camera(self, dataId):
        """Map a camera dataset."""
        if self.cameraPolicyLocation is None:
            raise RuntimeError, "No camera dataset available."
        actualId = self._transformId(dataId)
        return ButlerLocation("lsst.afw.cameraGeom.Camera", "Camera",
                "PafStorage", self.cameraPolicyLocation, actualId)

    def std_camera(self, item, dataId):
        """Standardize a camera dataset by converting it to a camera
        object."""
        return cameraGeomUtils.makeCamera(cameraGeomUtils.getGeomPolicy(item))

    def std_raw(self, item, dataId):
        """Standardize a raw dataset by converting it to an Exposure instead of an Image"""
        item = exposureFromImage(item)
        return self._standardizeExposure(self.exposures['raw'], item, dataId)

    def map_skypolicy(self, dataId):
        """Map a sky policy."""
        return ButlerLocation("lsst.pex.policy.Policy", "Policy",
                "Internal", None, None)

    def std_skypolicy(self, item, dataId):
        """Standardize a sky policy by returning the one we use."""
        return self.skypolicy

###############################################################################
#
# Utility functions
#
###############################################################################

    def _setupRegistry(self, name, path, policy, policyKey, root):
        """Set up a registry (usually SQLite3), trying a number of possible
        paths.
        @param name       (string) Name of registry
        @param path       (string) Path for registry
        @param policyKey  (string) Key in policy for registry path
        @param root       (string) Root directory to look in
        @return (lsst.daf.butlerUtils.Registry) Registry object"""

        if path is None and policy.exists(policyKey):
            path = dafPersist.LogicalLocation(
                    policy.getString(policyKey)).locString()
            if not os.path.exists(path):
                if not os.path.isabs(path) and root is not None:
                    rootedPath = os.path.join(root, path)
                    if not os.path.exists(rootedPath):
                        self.log.log(pexLog.Log.WARN,
                                "Unable to locate registry at policy path (also looked in root): %s" % path)
                        path = None
                    else:
                        path = rootedPath
                else:
                    self.log.log(pexLog.Log.WARN,
                            "Unable to locate registry at policy path: %s" % path)
                    path = None
        if path is None and root is not None:
            path = os.path.join(root, "%s.sqlite3" % name)
            if not os.path.exists(path):
                self.log.log(pexLog.Log.WARN,
                        "Unable to locate %s registry in root: %s" % (name, path))
                path = None
        if path is None:
            path = "%s.sqlite3" % name
            if not os.path.exists(path):
                self.log.log(pexLog.Log.WARN,
                        "Unable to locate %s registry in current dir: %s" % (name, path))
                path = None
        if path is not None:
            self.log.log(pexLog.Log.INFO,
                    "Loading %s registry from %s" % (name, path))
            registry = Registry.create(path)
            if registry is None:
                raise RuntimeError, "Unable to load %s registry from %s" % (name, path)
            return registry
        else:
            # TODO Try a FsRegistry(root)
            self.log.log(pexLog.Log.WARN,
                         "No registry loaded; proceeding without one")
            return None

    def _transformId(self, dataId):
        """Transform an id from camera-specific usage to standard form (e.g,.
        ccdname --> ccd).  The default implementation merely copies its input.
        @param dataId[in] (dict) Dataset identifier
        @return (dict) Transformed dataset identifier"""

        return dataId.copy()

    def _mapActualToPath(self, template, actualId):
        """Convert a template path to an actual path, using the actual data
        identifier.  This implementation is usually sufficient but can be
        overridden by the subclass.
        @param template (string) Template path
        @param actualId (dict) Dataset identifier
        @return (string) Pathname"""

        return template % self._transformId(actualId)

    def _extractDetectorName(self, dataId):
        """Extract the detector (CCD) name (used with the afw camera geometry
        class) from the dataset identifier.  This function must be provided by
        the subclass.
        @param dataId (dict) Dataset identifier
        @return (string) Detector name"""

        raise RuntimeError, "No _extractDetectorName() function specified"

    def _extractAmpId(self, dataId):
        """Extract the amplifier identifer from a dataset identifier.  The
        amplifier identifier has three parts: the detector name for the CCD
        containing the amplifier and x and y coordinates for the amplifier
        within the CCD.  The default implementation just uses the amplifier
        number, which is assumed to be the same as the x coordinate.
        @param dataId (dict) Dataset identifer
        @return (tuple) Amplifier identifier"""

        return (self._extractDetectorName(dataId),
                int(dataId['amp']), 0)

    def _setAmpDetector(self, item, dataId):
        """Set the detector object in an Exposure for an amplifier.
        Defects are also added to the Exposure based on the detector object.
        @param[in,out] item (lsst.afw.image.Exposure)
        @param dataId (dict) Dataset identifier"""

        ampId = self._extractAmpId(dataId)
        detector = cameraGeomUtils.findAmp(
                self.camera, afwCameraGeom.Id(ampId[0]), ampId[1], ampId[2])
        self._addDefects(dataId, amp=detector)
        item.setDetector(detector)

    def _setCcdDetector(self, item, dataId):
        """Set the detector object in an Exposure for a CCD.
        Defects are also added to the Exposure based on the detector object.
        @param[in,out] item (lsst.afw.image.Exposure)
        @param dataId (dict) Dataset identifier"""

        ccdId = self._extractDetectorName(dataId)
        detector = cameraGeomUtils.findCcd(self.camera, afwCameraGeom.Id(ccdId))
        self._addDefects(dataId, ccd=detector)
        item.setDetector(detector)

    def _setFilter(self, mapping, item, dataId):
        """Set the filter object in an Exposure.  If the Exposure had a FILTER
        keyword, this was already processed during load.  But if it didn't,
        use the filter from the registry.
        @param mapping (lsst.daf.butlerUtils.Mapping)
        @param[in,out] item (lsst.afw.image.Exposure)
        @param dataId (dict) Dataset identifier"""

        if not (isinstance(item, afwImage.ExposureU) or isinstance(item, afwImage.ExposureI) or
                isinstance(item, afwImage.ExposureF) or isinstance(item, afwImage.ExposureD)):
            return

        actualId = mapping.need(['filter'], dataId)
        filterName = actualId['filter']
        if self.filters is not None and self.filters.has_key(filterName):
            filterName = self.filters[filterName]
        item.setFilter(afwImage.Filter(filterName))

    def _setTimes(self, mapping, item, dataId):
        """Set the exposure time and exposure midpoint in the calib object in
        an Exposure.  Use the EXPTIME and MJD-OBS keywords (and strip out
        EXPTIME).
        @param mapping (lsst.daf.butlerUtils.Mapping)
        @param[in,out] item (lsst.afw.image.Exposure)
        @param dataId (dict) Dataset identifier"""

        md = item.getMetadata()
        calib = item.getCalib()
        if md.exists("EXPTIME"):
            expTime = md.get("EXPTIME")
            calib.setExptime(expTime)
            md.remove("EXPTIME")
        else:
            expTime = calib.getExptime()
        if md.exists("MJD-OBS"):
            obsStart = dafBase.DateTime(md.get("MJD-OBS"),
                    dafBase.DateTime.MJD, dafBase.DateTime.UTC)
            obsMidpoint = obsStart.nsecs() + long(expTime * 1000000000L / 2)
            calib.setMidTime(dafBase.DateTime(obsMidpoint))


    # Default standardization function for exposures
    def _standardizeExposure(self, mapping, item, dataId, filter=True):
        """Default standardization function for images.
        @param mapping (lsst.daf.butlerUtils.Mapping)
        @param[in,out] item (lsst.afw.image.Exposure)
        @param dataId (dict) Dataset identifier
        @param filter (bool) Set filter?
        @return (lsst.afw.image.Exposure) the standardized Exposure"""

        if (re.search(r'Exposure', mapping.python) and re.search(r'Image',mapping.persistable)):
            item = exposureFromImage(item)

        if mapping.level.lower() == "amp":
            self._setAmpDetector(item, dataId)
        elif mapping.level.lower() == "ccd":
            self._setCcdDetector(item, dataId)

        if filter:
            self._setFilter(mapping, item, dataId)
        if not isinstance(mapping, CalibrationMapping):
            self._setTimes(mapping, item, dataId)

        return item

    def _defectLookup(self, dataId, ccdSerial):
        """Find the defects for a given CCD.
        @param dataId (dict) Dataset identifier
        @param ccdSerial (string) CCD serial number
        @return (string) path to the defects file or None if not available"""

        if self.defectRegistry is None:
            return None
        if self.registry is None:
            raise RuntimeError, "No registry for defect lookup"

        rows = self.registry.executeQuery(("taiObs",), ("raw_visit",),
                [("visit", "?")], None, (dataId['visit'],))
        if len(rows) == 0:
            return None
        assert len(rows) == 1
        taiObs = rows[0][0]

        # Lookup the defects for this CCD serial number that are valid at the
        # exposure midpoint.
        rows = self.defectRegistry.executeQuery(("path",), ("defect",),
                [("ccdSerial", "?")],
                ("DATETIME(?)", "DATETIME(validStart)", "DATETIME(validEnd)"),
                (ccdSerial, taiObs))
        if not rows or len(rows) == 0:
            return None
        assert len(rows) == 1
        return os.path.join(self.defectPath, rows[0][0])

    def _addDefects(self, dataId, amp=None, ccd=None):
        """Add the defects for an amplifier or a CCD to the detector object
        for that amplifier/CCD.
        @param dataId (dict) Dataset identifier
        @param[in,out] amp (lsst.afw.cameraGeom.Detector)
        @param[in,out] ccd (lsst.afw.cameraGeom.Detector)
        Exactly one of amp or ccd should be specified."""

        if ccd is None:
            ccd = afwCameraGeom.cast_Ccd(amp.getParent())
        if len(ccd.getDefects()) > 0:
            # Assume we have loaded them properly already
            return
        defectFits = self._defectLookup(dataId, ccd.getId().getSerial())
        if defectFits is not None:
            defectDict = cameraGeomUtils.makeDefectsFromFits(defectFits)
            ccdDefects = None
            for k in defectDict.keys():
                if k == ccd.getId():
                    ccdDefects = defectDict[k]
                    break
            if ccdDefects is None:
                raise RuntimeError, "No defects for ccd %s in %s" % \
                        (str(ccd.getId()), defectFits)
            ccd.setDefects(ccdDefects)


def exposureFromImage(image):
    """Generate an exposure from a DecoratedImage or similar
    @param[in] image Image of interest
    @return (lsst.afw.image.Exposure) Exposure containing input image
    """
    if isinstance(image, afwImage.DecoratedImageU) or isinstance(image, afwImage.DecoratedImageI) or \
        isinstance(image, afwImage.DecoratedImageF) or isinstance(image, afwImage.DecoratedImageD):
        exposure = afwImage.makeExposure(afwImage.makeMaskedImage(image.getImage()))
    else:
        exposure = image
    md = image.getMetadata()
    exposure.setMetadata(md)
    wcs = afwImage.makeWcs(md)
    exposure.setWcs(wcs)
    wcsMetadata = wcs.getFitsMetadata()
    for kw in wcsMetadata.paramNames():
        md.remove(kw)

    return exposure
