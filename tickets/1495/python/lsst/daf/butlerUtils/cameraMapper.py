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
from lsst.daf.butlerUtils import Mapping, CalibrationMapping, Registry
import lsst.daf.base as dafBase
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

    Public methods: getKeys, queryMetadata, getDatasetTypes, map,
    canStandardize, standardize, getMapping

    Mappers for specific data sources (e.g., CFHT Megacam, LSST
    simulations, etc.) should inherit this class.  Such subclasses
    should set in __init__:

    keys: List of keys that can be used in data IDs.

    filterMap: Dict with mapping from data's filter name (e.g.,
    "r.12345") to the code's filter name (e.g., "r")

    filterIdMap: Mapping from the code's filter name (e.g., "r") to an
    integer identifier.

    The following method must be provided by the subclass:

    _extractDetectorName(self, dataId): returns the detector name for a CCD
    (e.g., "CFHT 21", "R:1,2 S:3,4") given a dataset identifier referring to
    that CCD or a subcomponent of it.

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

    camera (string): Path to camera geometry policy file from subclassing
    module's repository

    defects (string): Path to defects directory from subclassing module's
    repository

    filters (string): Path to filters policy file from subclassing module's
    repository

    exposures (policy): Exposure mappings (e.g., "raw", "postISR")

    calibrations (policy): Calibration mappings (e.g., "bias", "flat")

    The 'exposures' and 'calibrations' policies consist of mappings
    (see Mappings class).

    Functions to map (provide a path to the data given a dataset
    identifier dictionary) and standardize (convert data into some standard
    format or type) may be provided in the subclass as "map_{dataset type}"
    and "std_{dataset type}", respectively.

    Implementations of map_camera and std_camera that should typically be
    sufficient are provided in this base class.
    """

    def __init__(self, policy, repositoryDir,
            root=None, registry=None, calibRoot=None, calibRegistry=None):
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
                             metadata"""

        dafPersist.Mapper.__init__(self)

        self.log = pexLog.Log(pexLog.getDefaultLog(), "CameraMapper")

        self.policy = policy

        # Dictionary
        dictFile = pexPolicy.DefaultPolicyFile("daf_butlerUtils",
                "MapperDictionary.paf", "policy")
        dictPolicy = pexPolicy.Policy.createPolicy(dictFile,
                dictFile.getRepositoryPath())
        self.policy.mergeDefaults(dictPolicy)

        # Root directories
        self.root = root
        if self.root is None:
            if self.policy.exists('root'):
                self.root = self.policy.getString('root')
            else:
                self.root = "."
        self.calibRoot = calibRoot
        if self.calibRoot is None:
            if self.policy.exists('calibRoot'):
                self.calibRoot = self.policy.getString('calibRoot')
            else:
                self.calibRoot = self.root
        # Do any location substitutions
        self.root = dafPersist.LogicalLocation(self.root).locString()
        self.calibRoot = dafPersist.LogicalLocation(self.calibRoot).locString()
        if not os.path.exists(self.root):
            self.log.log(pexLog.Log.WARN,
                    "Root directory not found: %s" % (root,))
        if not os.path.exists(self.calibRoot):
            self.log.log(pexLog.Log.WARN,
                    "Calibration root directory not found: %s" % (calibRoot,))

        # Registries
        self.registry = self._setupRegistry("registry", registry, "registryPath", root)
        if self.policy.exists('needCalibRegistry') and \
                self.policy.getBool('needCalibRegistry'):
            self.calibRegistry = self._setupRegistry(
                    "calibRegistry", calibRegistry,
                    "calibRegistryPath", calibRoot)
        else:
            self.calibRegistry = None

        # Sub-dictionary (for exposure/calibration types)
        mappingFile = pexPolicy.DefaultPolicyFile("daf_butlerUtils",
                "MappingDictionary.paf", "policy")
        mappingPolicy = pexPolicy.Policy.createPolicy(mappingFile,
                mappingFile.getRepositoryPath())

        # Mappings
        self.mappings = dict()
        if self.policy.exists("exposures"):
            exposures = self.policy.getPolicy("exposures") # List of exposure types
            for datasetType in exposures.names(True):
                subPolicy = exposures.getPolicy(datasetType)
                subPolicy.mergeDefaults(mappingPolicy)
                self.mappings[datasetType] = Mapping(
                        mapper=self, policy=subPolicy, datasetType=datasetType,
                        registry=self.registry, root=self.root)
        if self.policy.exists("calibrations"):
            calibs = self.policy.getPolicy("calibrations") # List of calibration types
            for datasetType in calibs.names(True):
                subPolicy = calibs.getPolicy(datasetType)
                subPolicy.mergeDefaults(mappingPolicy)
                self.mappings[datasetType] = CalibrationMapping(
                        mapper=self, policy=subPolicy, datasetType=datasetType,
                        registry=self.calibRegistry, root=self.calibRoot)

        # Subclass should override these!
        self.keys = []
        self.filterMap = {}
        self.filterIdMap = {}

        # Camera geometry
        self.cameraPolicyLocation = None
        self.camera = None
        if self.policy.exists('camera'):
            cameraPolicyLocation = self.policy.getString('camera')

            # must be explicit for ButlerLocation later
            self.cameraPolicyLocation = os.path.join(
                    repositoryDir, cameraPolicyLocation)

            cameraPolicy = pexPolicy.Policy.createPolicy(cameraPolicyLocation,
                    repositoryDir)
            cameraPolicy = cameraGeomUtils.getGeomPolicy(cameraPolicy)
            self.camera = cameraGeomUtils.makeCamera(cameraPolicy)

        # Defect registry and root
        self.defectRegistry = None
        if self.policy.exists('defects'):
            self.defectPath = os.path.join(
                    repositoryDir, self.policy.getString('defects'))
            defectRegistryLocation = os.path.join(
                    self.defectPath, "defectRegistry.sqlite3")
            self.defectRegistry = \
                    butlerUtils.Registry.create(defectRegistryLocation)

        # Filters
        if self.policy.exists('filters'):
            filterPolicyLocation = self.policy.getString('filters')
            filterPolicy = pexPolicy.Policy.createPolicy(
                    filterPolicyLocation, repositoryDir)
            imageUtils.defineFiltersFromPolicy(filterPolicy, reset=True)


    def getKeys(self):
        """Return supported keys.
        @return (iterable) List of keys usable in a dataset identifier"""
        return self.keys

    def getMapping(self, datasetType):
        """Return the appropriate mapping for the dataset type.
        @param datasetType (string)
        @return (daf.butlerUtils.Mapping) Mapping object
        """
        return self.mappings[datasetType]

    def map_camera(self, datasetType, dataId):
        """Map a camera dataset."""
        if self.cameraPolicyLocation is None:
            raise RuntimeError, "No camera dataset available."
        actualId = self._transformId(dataId)
        return ButlerLocation("lsst.afw.cameraGeom.Camera", "Camera",
                "PafStorage", self.cameraPolicyLocation, actualId)

    def std_camera(self, datasetType, item, dataId):
        """Standardize a camera dataset by converting it to a camera
        object."""
        return cameraGeomUtils.makeCamera(cameraGeomUtils.getGeomPolicy(item))

###############################################################################
#
# Utility functions
#
###############################################################################

    def _setupRegistry(self, name, path, policyKey, root):
        """Set up a registry (usually SQLite3), trying a number of possible
        paths.
        @param name       (string) Name of registry
        @param path       (string) Path for registry
        @param policyKey  (string) Key in policy for registry path
        @param root       (string) Root directory to look in
        @return (lsst.daf.butlerUtils.Registry) Registry object"""

        if path is None and self.policy.exists(policyName):
            path = dafPersist.LogicalLocation(
                    self.policy.getString(policyName)).locString()
            if not os.path.exists(path):
                self.log.log(pexLog.Log.WARN,
                        "Unable to locate registry at path: %s" % path)
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

        return template % actualId

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
        """Set the filter object in an Exposure.  Use the filter specified by
        the FILTER keyword (and strip it out) or the filter obtained from the
        registry.
        @param mapping (lsst.daf.butlerUtils.Mapping)
        @param[in,out] item (lsst.afw.image.Exposure)
        @param dataId (dict) Dataset identifier"""

        md = item.getMetadata()
        filterName = None
        if md.exists("FILTER"):
            filterName = item.getMetadata().get("FILTER").strip()
            if self.filterMap.has_key(filterName):
                filterName = self.filterMap[filterName]
        if filterName is None:
            actualId = mapping.need(self, ['filter'], dataId)
            filterName = actualId['filter']
        filter = afwImage.Filter(filterName)
        item.setFilter(filter)

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

    # Default standardization function
    def _standardize(self, mapping, item, dataId):
        """Default standardization function for images.
        @param mapping (lsst.daf.butlerUtils.Mapping)
        @param[in,out] item (lsst.afw.image.Exposure)
        @param dataId (dict) Dataset identifier
        @return (lsst.afw.image.Exposure) the standardized Exposure"""

        stripFits(item.getMetadata())

        if mapping.level.lower() == "amp":
            self._setAmpDetector(item, dataId)
        elif mapping.level.lower() == "ccd":
            self._setCcdDetector(item, dataId)

        if not isinstance(mapping, CalibrationMapping):
            self._setTimes(mapping, item, dataId)
            self._setFilter(mapping, item, dataId)
        elif mapping.type in ['flat', 'fringe']:
            self._setFilter(mapping, item, dataId)
        return item

    def _defectLookup(self, dataId, ccdSerial):
        """Find the defects for a given CCD.
        @param dataId (dict) Dataset identifier
        @param ccdSerial (string) CCD serial number
        @return (string) path to the defects file or None if not available"""

        if self.defectRegistry is None:
            return None

        rows = self.registry.executeQuery(("taiObs",), ("raw_visit",),
                {"visit": "?"}, None, (dataId['visit'],))
        if len(rows) == 0:
            return None
        assert len(rows) == 1
        taiObs = rows[0][0]

        # Lookup the defects for this CCD serial number that are valid at the
        # exposure midpoint.
        rows = self.defectRegistry.executeQuery(("path",), ("defect",),
                {"ccdSerial": "?"},
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

def stripFits(propertySet):
    """Remove FITS-specific keywords from a flexible metadata PropertySet
    that have values maintained in an Exposure's fixed metadata.
    @param[in,out] propertySet (lsst.daf.base.PropertySet) Flexible metadata"""

    for kw in ("SIMPLE", "BITPIX", "EXTEND", "NAXIS", "NAXIS1", "NAXIS2",
            "BSCALE", "BZERO"):
        if propertySet.exists(kw):
            propertySet.remove(kw)
