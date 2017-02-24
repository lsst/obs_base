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

from builtins import str
import copy
import os
import pyfits  # required by _makeDefectsDict until defects are written as AFW tables
import re
import weakref
import lsst.daf.persistence as dafPersist
from . import ImageMapping, ExposureMapping, CalibrationMapping, DatasetMapping
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.afw.cameraGeom as afwCameraGeom
import lsst.log as lsstLog
import lsst.pex.policy as pexPolicy
from .exposureIdInfo import ExposureIdInfo
from .makeRawVisitInfo import MakeRawVisitInfo
from lsst.utils import getPackageDir

"""This module defines the CameraMapper base class."""


class CameraMapper(dafPersist.Mapper):

    """CameraMapper is a base class for mappers that handle images from a
    camera and products derived from them.  This provides an abstraction layer
    between the data on disk and the code.

    Public methods: keys, queryMetadata, getDatasetTypes, map,
    canStandardize, standardize

    Mappers for specific data sources (e.g., CFHT Megacam, LSST
    simulations, etc.) should inherit this class.

    The CameraMapper manages datasets within a "root" directory. Note that
    writing to a dataset present in the input root will hide the existing
    dataset but not overwrite it.  See #2160 for design discussion.

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

    Subclasses will typically set MakeRawVisitInfoClass:

    MakeRawVisitInfoClass: a class variable that points to a subclass of
    MakeRawVisitInfo, a functor that creates an
    lsst.afw.image.VisitInfo from the FITS metadata of a raw image.

    Subclasses must provide the following methods:

    _extractDetectorName(self, dataId): returns the detector name for a CCD
    (e.g., "CFHT 21", "R:1,2 S:3,4") as used in the AFW CameraGeom class given
    a dataset identifier referring to that CCD or a subcomponent of it.

    _computeCcdExposureId(self, dataId): see below

    _computeCoaddExposureId(self, dataId, singleFilter): see below

    Subclasses may also need to override the following methods:

    _transformId(self, dataId): transformation of a data identifier
    from colloquial usage (e.g., "ccdname") to proper/actual usage (e.g., "ccd"),
    including making suitable for path expansion (e.g. removing commas).
    The default implementation does nothing.  Note that this
    method should not modify its input parameter.

    getShortCcdName(self, ccdName): a static method that returns a shortened name
    suitable for use as a filename. The default version converts spaces to underscores.

    _getCcdKeyVal(self, dataId): return a CCD key and value
    by which to look up defects in the defects registry.
    The default value returns ("ccd", detector name)

    _mapActualToPath(self, template, actualId): convert a template path to an
    actual path, using the actual dataset identifier.

    The mapper's behaviors are largely specified by the policy file.
    See the MapperDictionary.paf for descriptions of the available items.

    The 'exposures', 'calibrations', and 'datasets' subpolicies configure
    mappings (see Mappings class).

    Common default mappings for all subclasses can be specified in the
    "policy/{images,exposures,calibrations,datasets}.yaml" files. This provides
    a simple way to add a product to all camera mappers.

    Functions to map (provide a path to the data given a dataset
    identifier dictionary) and standardize (convert data into some standard
    format or type) may be provided in the subclass as "map_{dataset type}"
    and "std_{dataset type}", respectively.

    If non-Exposure datasets cannot be retrieved using standard
    daf_persistence methods alone, a "bypass_{dataset type}" function may be
    provided in the subclass to return the dataset instead of using the
    "datasets" subpolicy.

    Implementations of map_camera and bypass_camera that should typically be
    sufficient are provided in this base class.

    @todo
    * Handle defects the same was as all other calibration products, using the calibration registry
    * Instead of auto-loading the camera at construction time, load it from the calibration registry
    * Rewrite defects as AFW tables so we don't need pyfits to unpersist them; then remove all mention
      of pyfits from this package.
    """
    packageName = None

    # a class or subclass of MakeRawVisitInfo, a functor that makes an
    # lsst.afw.image.VisitInfo from the FITS metadata of a raw image
    MakeRawVisitInfoClass = MakeRawVisitInfo

    def __init__(self, policy, repositoryDir,
                 root=None, registry=None, calibRoot=None, calibRegistry=None,
                 provided=None, parentRegistry=None, repositoryCfg=None):
        """Initialize the CameraMapper.

        Parameters
        ----------
        policy : daf_persistence.Policy,
            Can also be pexPolicy.Policy, only for backward compatibility.
            Policy with per-camera defaults already merged.
        repositoryDir : string
            Policy repository for the subclassing module (obtained with
            getRepositoryPath() on the per-camera default dictionary).
        root : string, optional
            Path to the root directory for data.
        registry : string, optional
            Path to registry with data's metadata.
        calibRoot : string, optional
            Root directory for calibrations.
        calibRegistry : string, optional
            Path to registry with calibrations' metadata.
        provided : list of string, optional
            Keys provided by the mapper.
        parentRegistry : Registry subclass, optional
            Registry from a parent repository that may be used to look up
            data's metadata.
        repositoryCfg : daf_persistence.RepositoryCfg or None, optional
            The configuration information for the repository this mapper is
            being used with.
        """

        dafPersist.Mapper.__init__(self)

        self.log = lsstLog.Log.getLogger("CameraMapper")

        self.root = root if root else repositoryCfg.root
        if isinstance(policy, pexPolicy.Policy):
            policy = dafPersist.Policy(policy)

        repoPolicy = repositoryCfg.policy if repositoryCfg else None
        if repoPolicy is not None:
            policy.update(repoPolicy)

        defaultPolicyFile = dafPersist.Policy.defaultPolicyFile("obs_base",
                                                                "MapperDictionary.paf",
                                                                "policy")
        dictPolicy = dafPersist.Policy(defaultPolicyFile)
        policy.merge(dictPolicy)

        # Levels
        self.levels = dict()
        if 'levels' in policy:
            levelsPolicy = policy['levels']
            for key in levelsPolicy.names(True):
                self.levels[key] = set(levelsPolicy.asArray(key))
        self.defaultLevel = policy['defaultLevel']
        self.defaultSubLevels = dict()
        if 'defaultSubLevels' in policy:
            self.defaultSubLevels = policy['defaultSubLevels']

        # Root directories
        if root is None:
            root = "."
        root = dafPersist.LogicalLocation(root).locString()

        self.rootStorage = dafPersist.Storage.makeFromURI(uri=root)

        # If the calibRoot is passed in, use that. If not and it's indicated in the policy, use that. And
        # otherwise, the calibs are in the regular root.
        if calibRoot is not None:
            calibStorage = dafPersist.Storage.makeFromURI(uri=calibRoot)
        elif 'calibRoot' in policy:
            calibRoot = policy['calibRoot']
            calibRoot = dafPersist.LogicalLocation(calibRoot).locString()
            calibStorage = dafPersist.Storage.makeFromURI(uri=calibRoot)
        else:
            calibStorage = self.rootStorage

        self.root = root

        # Registries
        self.registry = self._setupRegistry("registry", registry, policy, "registryPath", self.rootStorage,
                                            searchParents=False, posixIfNoSql=(not parentRegistry))
        if not self.registry:
            self.registry = parentRegistry
        needCalibRegistry = policy.get('needCalibRegistry', None)
        if needCalibRegistry:
            if calibStorage:
                self.calibRegistry = self._setupRegistry("calibRegistry", calibRegistry, policy,
                                                         "calibRegistryPath", calibStorage)
            else:
                raise RuntimeError(
                    "'needCalibRegistry' is true in Policy, but was unable to locate a repo at " +
                    "calibRoot ivar:%s or policy['calibRoot']:%s" %
                    (calibRoot, policy.get('calibRoot', None)))
        else:
            self.calibRegistry = None

        # Dict of valid keys and their value types
        self.keyDict = dict()

        self._initMappings(policy, self.rootStorage, calibStorage, provided=None)

        # Camera geometry
        self.cameraDataLocation = None  # path to camera geometry config file
        self.camera = self._makeCamera(policy=policy, repositoryDir=repositoryDir)

        # Defect registry and root. Defects are stored with the camera and the registry is loaded from the
        # camera package, which is on the local filesystem.
        self.defectRegistry = None
        if 'defects' in policy:
            self.defectPath = os.path.join(repositoryDir, policy['defects'])
            defectRegistryLocation = os.path.join(self.defectPath, "defectRegistry.sqlite3")
            self.defectRegistry = dafPersist.Registry.create(defectRegistryLocation)

        # Filter translation table
        self.filters = None

        # Skytile policy
        self.skypolicy = policy['skytiles']

        # verify that the class variable packageName is set before attempting
        # to instantiate an instance
        if self.packageName is None:
            raise ValueError('class variable packageName must not be None')

        self.makeRawVisitInfo = self.MakeRawVisitInfoClass(log=self.log)

    def _initMappings(self, policy, rootStorage=None, calibStorage=None, provided=None):
        """Initialize mappings

        For each of the dataset types that we want to be able to read, there are
        methods that can be created to support them:
        * map_<dataset> : determine the path for dataset
        * std_<dataset> : standardize the retrieved dataset
        * bypass_<dataset> : retrieve the dataset (bypassing the usual retrieval machinery)
        * query_<dataset> : query the registry

        Besides the dataset types explicitly listed in the policy, we create
        additional, derived datasets for additional conveniences, e.g., reading
        the header of an image, retrieving only the size of a catalog.

        @param policy        (Policy) Policy with per-camera defaults already merged
        @param rootStorage   (Storage subclass instance) Interface to persisted repository data
        @param calibRoot     (Storage subclass instance) Interface to persisted calib repository data
        @param provided      (list of strings) Keys provided by the mapper
        """
        # Sub-dictionaries (for exposure/calibration/dataset types)
        imgMappingPolicy = dafPersist.Policy(dafPersist.Policy.defaultPolicyFile(
            "obs_base", "ImageMappingDictionary.paf", "policy"))
        expMappingPolicy = dafPersist.Policy(dafPersist.Policy.defaultPolicyFile(
            "obs_base", "ExposureMappingDictionary.paf", "policy"))
        calMappingPolicy = dafPersist.Policy(dafPersist.Policy.defaultPolicyFile(
            "obs_base", "CalibrationMappingDictionary.paf", "policy"))
        dsMappingPolicy = dafPersist.Policy(dafPersist.Policy.defaultPolicyFile(
            "obs_base", "DatasetMappingDictionary.paf", "policy"))

        # Mappings
        mappingList = (
            ("images", imgMappingPolicy, ImageMapping),
            ("exposures", expMappingPolicy, ExposureMapping),
            ("calibrations", calMappingPolicy, CalibrationMapping),
            ("datasets", dsMappingPolicy, DatasetMapping)
        )
        self.mappings = dict()
        for name, defPolicy, cls in mappingList:
            if name in policy:
                datasets = policy[name]

                # Centrally-defined datasets
                defaultsPath = os.path.join(getPackageDir("obs_base"), "policy", name + ".yaml")
                if os.path.exists(defaultsPath):
                    datasets.merge(dafPersist.Policy(defaultsPath))

                mappings = dict()
                setattr(self, name, mappings)
                for datasetType in datasets.names(True):
                    subPolicy = datasets[datasetType]
                    subPolicy.merge(defPolicy)

                    if not hasattr(self, "map_" + datasetType) and 'composite' in subPolicy:
                        def compositeClosure(dataId, write=False, mapper=None, mapping=None, subPolicy=subPolicy):
                            components = subPolicy.get('composite')
                            assembler = subPolicy['assembler'] if 'assembler' in subPolicy else None
                            disassembler = subPolicy['disassembler'] if 'disassembler' in subPolicy else None
                            python = subPolicy['python']
                            butlerComposite = dafPersist.ButlerComposite(assembler=assembler,
                                                                         disassembler=disassembler,
                                                                         python=python,
                                                                         dataId=dataId,
                                                                         mapper=self)
                            for name, component in components.items():
                                butlerComposite.add(id=name,
                                                    datasetType=component.get('datasetType'),
                                                    setter=component.get('setter', None),
                                                    getter=component.get('getter', None),
                                                    subset=component.get('subset', False),
                                                    inputOnly=component.get('inputOnly', False))
                            return butlerComposite
                        setattr(self, "map_" + datasetType, compositeClosure)
                        # for now at least, don't set up any other handling for this dataset type.
                        continue

                    if name == "calibrations":
                        mapping = cls(datasetType, subPolicy, self.registry, self.calibRegistry, calibStorage,
                                      provided=provided)
                    else:
                        mapping = cls(datasetType, subPolicy, self.registry, rootStorage, provided=provided)
                    self.keyDict.update(mapping.keys())
                    mappings[datasetType] = mapping
                    self.mappings[datasetType] = mapping
                    if not hasattr(self, "map_" + datasetType):
                        def mapClosure(dataId, write=False, mapper=weakref.proxy(self), mapping=mapping):
                            return mapping.map(mapper, dataId, write)
                        setattr(self, "map_" + datasetType, mapClosure)
                    if not hasattr(self, "query_" + datasetType):
                        def queryClosure(format, dataId, mapping=mapping):
                            return mapping.lookup(format, dataId)
                        setattr(self, "query_" + datasetType, queryClosure)
                    if hasattr(mapping, "standardize") and not hasattr(self, "std_" + datasetType):
                        def stdClosure(item, dataId, mapper=weakref.proxy(self), mapping=mapping):
                            return mapping.standardize(mapper, item, dataId)
                        setattr(self, "std_" + datasetType, stdClosure)

                    def setMethods(suffix, mapImpl=None, bypassImpl=None, queryImpl=None):
                        """Set convenience methods on CameraMapper"""
                        mapName = "map_" + datasetType + "_" + suffix
                        bypassName = "bypass_" + datasetType + "_" + suffix
                        queryName = "query_" + datasetType + "_" + suffix
                        if not hasattr(self, mapName):
                            setattr(self, mapName, mapImpl or getattr(self, "map_" + datasetType))
                        if not hasattr(self, bypassName):
                            if bypassImpl is None and hasattr(self, "bypass_" + datasetType):
                                bypassImpl = getattr(self, "bypass_" + datasetType)
                            if bypassImpl is not None:
                                setattr(self, bypassName, bypassImpl)
                        if not hasattr(self, queryName):
                            setattr(self, queryName, queryImpl or getattr(self, "query_" + datasetType))

                    # Filename of dataset
                    setMethods("filename", bypassImpl=lambda datasetType, pythonType, location, dataId:
                        [os.path.join(location.getStorage().root, p) for p in location.getLocations()])

                    # Metadata from FITS file
                    if subPolicy["storage"] == "FitsStorage":  # a FITS image
                        setMethods("md", bypassImpl=lambda datasetType, pythonType, location, dataId:
                                   afwImage.readMetadata(location.getLocationsWithRoot()[0]))
                        if name == "exposures":
                            setMethods("wcs", bypassImpl=lambda datasetType, pythonType, location, dataId:
                                       afwImage.makeWcs(
                                           afwImage.readMetadata(location.getLocationsWithRoot()[0])))
                            setMethods("calib", bypassImpl=lambda datasetType, pythonType, location, dataId:
                                       afwImage.Calib(
                                           afwImage.readMetadata(location.getLocationsWithRoot()[0])))
                            setMethods("visitInfo",
                                       bypassImpl=lambda datasetType, pythonType, location, dataId:
                                       afwImage.VisitInfo(
                                           afwImage.readMetadata(location.getLocationsWithRoot()[0])))
                    if subPolicy["storage"] == "FitsCatalogStorage":  # a FITS catalog
                        setMethods("md", bypassImpl=lambda datasetType, pythonType, location, dataId:
                                   afwImage.readMetadata(os.path.join(location.getStorage().root,
                                                                      location.getLocations()[0]), 2))

                    # Sub-images
                    if subPolicy["storage"] == "FitsStorage":
                        def mapSubClosure(dataId, write=False, mapper=weakref.proxy(self), mapping=mapping):
                            subId = dataId.copy()
                            del subId['bbox']
                            loc = mapping.map(mapper, subId, write)
                            bbox = dataId['bbox']
                            llcX = bbox.getMinX()
                            llcY = bbox.getMinY()
                            width = bbox.getWidth()
                            height = bbox.getHeight()
                            loc.additionalData.set('llcX', llcX)
                            loc.additionalData.set('llcY', llcY)
                            loc.additionalData.set('width', width)
                            loc.additionalData.set('height', height)
                            if 'imageOrigin' in dataId:
                                loc.additionalData.set('imageOrigin',
                                                       dataId['imageOrigin'])
                            return loc
                        def querySubClosure(key, format, dataId, mapping=mapping):
                            subId = dataId.copy()
                            del subId['bbox']
                            return mapping.lookup(format, subId)
                        setMethods("sub", mapImpl=mapSubClosure, queryImpl=querySubClosure)

                    if subPolicy["storage"] == "FitsCatalogStorage":
                        # Length of catalog
                        setMethods("len", bypassImpl=lambda datasetType, pythonType, location, dataId:
                                   afwImage.readMetadata(os.path.join(location.getStorage().root,
                                                                      location.getLocations()[0]),
                                                         2).get("NAXIS2"))

                        # Schema of catalog
                        if not datasetType.endswith("_schema") and datasetType + "_schema" not in datasets:
                            setMethods("schema", bypassImpl=lambda datasetType, pythonType, location, dataId:
                                       afwTable.Schema.readFits(os.path.join(location.getStorage().root,
                                                                             location.getLocations()[0])))

    def _computeCcdExposureId(self, dataId):
        """Compute the 64-bit (long) identifier for a CCD exposure.

        Subclasses must override

        @param dataId (dict) Data identifier with visit, ccd
        """
        raise NotImplementedError()

    def _computeCoaddExposureId(self, dataId, singleFilter):
        """Compute the 64-bit (long) identifier for a coadd.

        Subclasses must override

        @param dataId (dict)       Data identifier with tract and patch.
        @param singleFilter (bool) True means the desired ID is for a single-
                                   filter coadd, in which case dataId
                                   must contain filter.
        """
        raise NotImplementedError()

    def _search(self, path):
        """Search for path in the associated repository's storage.

        Parameters
        ----------
        path : string
            Path that describes an object in the repository associated with
            this mapper.
            Path may contain an HDU indicator, e.g. 'foo.fits[1]'. The
            indicator will be stripped when searching and so will match
            filenames without the HDU indicator, e.g. 'foo.fits'. The path
            returned WILL contain the indicator though, e.g. ['foo.fits[1]'].

        Returns
        -------
        string
            The path for this object in the repository. Will return None if the
            object can't be found. If the input argument path contained an HDU
            indicator, the returned path will also contain the HDU indicator.
        """
        # it would be better if storage was an instance, instead of having to demux the root URI every time.
        return dafPersist.Storage.search(self.root, path)

    def backup(self, datasetType, dataId):
        """Rename any existing object with the given type and dataId.

        The CameraMapper implementation saves objects in a sequence of e.g.:
          foo.fits
          foo.fits~1
          foo.fits~2
        All of the backups will be placed in the output repo, however, and will
        not be removed if they are found elsewhere in the _parent chain.  This
        means that the same file will be stored twice if the previous version was
        found in an input repo.
        """

        # Calling PosixStorage directly is not the long term solution in this
        # function, this is work-in-progress on epic DM-6225. The plan is for
        # parentSearch to be changed to 'search', and search only the storage
        # associated with this mapper. All searching of parents will be handled
        # by traversing the container of repositories in Butler.

        def firstElement(list):
            """Get the first element in the list, or None if that can't be done.
            """
            return list[0] if list is not None and len(list) else None

        n = 0
        newLocation = self.map(datasetType, dataId, write=True)
        newPath = newLocation.getLocations()[0]
        path = dafPersist.PosixStorage.search(self.root, newPath, searchParents=True)
        path = firstElement(path)
        oldPaths = []
        while path is not None:
            n += 1
            oldPaths.append((n, path))
            path = dafPersist.PosixStorage.search(self.root, "%s~%d" % (newPath, n), searchParents=True)
            path = firstElement(path)
        for n, oldPath in reversed(oldPaths):
            self.rootStorage.copyFile(oldPath, "%s~%d" % (newPath, n))

    def keys(self):
        """Return supported keys.
        @return (iterable) List of keys usable in a dataset identifier"""
        return iter(self.keyDict.keys())

    def getKeys(self, datasetType, level):
        """Return supported keys and their value types for a given dataset
        type at a given level of the key hierarchy.

        @param datasetType (str) dataset type or None for all keys
        @param level (str) level or None for all levels
        @return (iterable) Set of keys usable in a dataset identifier"""

        # not sure if this is how we want to do this. what if None was intended?
        if level == '':
            level = self.getDefaultLevel()

        if datasetType is None:
            keyDict = copy.copy(self.keyDict)
        else:
            keyDict = self.mappings[datasetType].keys()
        if level is not None and level in self.levels:
            keyDict = copy.copy(keyDict)
            for l in self.levels[level]:
                if l in keyDict:
                    del keyDict[l]
        return keyDict

    def getDefaultLevel(self):
        return self.defaultLevel

    def getDefaultSubLevel(self, level):
        if level in self.defaultSubLevels:
            return self.defaultSubLevels[level]
        return None

    @classmethod
    def getCameraName(cls):
        """Return the name of the camera that this CameraMapper is for."""
        className = str(cls)
        className = className[className.find('.'):-1]
        m = re.search(r'(\w+)Mapper', className)
        if m is None:
            m = re.search(r"class '[\w.]*?(\w+)'", className)
        name = m.group(1)
        return name[:1].lower() + name[1:] if name else ''

    @classmethod
    def getPackageName(cls):
        """Return the name of the package containing this CameraMapper."""
        if cls.packageName is None:
            raise ValueError('class variable packageName must not be None')
        return cls.packageName

    def map_camera(self, dataId, write=False):
        """Map a camera dataset."""
        if self.camera is None:
            raise RuntimeError("No camera dataset available.")
        actualId = self._transformId(dataId)
        return dafPersist.ButlerLocation(
            pythonType="lsst.afw.cameraGeom.CameraConfig",
            cppType="Config",
            storageName="ConfigStorage",
            locationList=self.cameraDataLocation or "ignored",
            dataId=actualId,
            mapper=self,
            storage=self.rootStorage
        )

    def bypass_camera(self, datasetType, pythonType, butlerLocation, dataId):
        """Return the (preloaded) camera object.
        """
        if self.camera is None:
            raise RuntimeError("No camera dataset available.")
        return self.camera

    def map_defects(self, dataId, write=False):
        """Map defects dataset.

        @return a very minimal ButlerLocation containing just the locationList field
            (just enough information that bypass_defects can use it).
        """
        defectFitsPath = self._defectLookup(dataId=dataId)
        if defectFitsPath is None:
            raise RuntimeError("No defects available for dataId=%s" % (dataId,))

        return dafPersist.ButlerLocation(None, None, None, defectFitsPath,
                                         dataId, self,
                                         storage=self.rootStorage)

    def bypass_defects(self, datasetType, pythonType, butlerLocation, dataId):
        """Return a defect based on the butler location returned by map_defects

        @param[in] butlerLocation: a ButlerLocation with locationList = path to defects FITS file
        @param[in] dataId: the usual data ID; "ccd" must be set

        Note: the name "bypass_XXX" means the butler makes no attempt to convert the ButlerLocation
        into an object, which is what we want for now, since that conversion is a bit tricky.
        """
        detectorName = self._extractDetectorName(dataId)
        defectsFitsPath = butlerLocation.locationList[0]
        with pyfits.open(defectsFitsPath) as hduList:
            for hdu in hduList[1:]:
                if hdu.header["name"] != detectorName:
                    continue

                defectList = []
                for data in hdu.data:
                    bbox = afwGeom.Box2I(
                        afwGeom.Point2I(int(data['x0']), int(data['y0'])),
                        afwGeom.Extent2I(int(data['width']), int(data['height'])),
                    )
                    defectList.append(afwImage.DefectBase(bbox))
                return defectList

        raise RuntimeError("No defects for ccd %s in %s" % (detectorName, defectsFitsPath))

    def map_expIdInfo(self, dataId, write=False):
        return dafPersist.ButlerLocation(
            pythonType="lsst.obs.base.ExposureIdInfo",
            cppType=None,
            storageName="Internal",
            locationList="ignored",
            dataId=dataId,
            mapper=self,
            storage=self.rootStorage
        )

    def bypass_expIdInfo(self, datasetType, pythonType, location, dataId):
        """Hook to retrieve an lsst.obs.base.ExposureIdInfo for an exposure"""
        expId = self.bypass_ccdExposureId(datasetType, pythonType, location, dataId)
        expBits = self.bypass_ccdExposureId_bits(datasetType, pythonType, location, dataId)
        return ExposureIdInfo(expId=expId, expBits=expBits)

    def std_bfKernel(self, item, dataId):
        """Disable standardization for bfKernel

        bfKernel is a calibration product that is numpy array,
        unlike other calibration products that are all images;
        all calibration images are sent through _standardizeExposure
        due to CalibrationMapping, but we don't want that to happen to bfKernel
        """
        return item

    def std_raw(self, item, dataId):
        """Standardize a raw dataset by converting it to an Exposure instead of an Image"""
        exposure = exposureFromImage(item)
        exposureId = self._computeCcdExposureId(dataId)
        md = exposure.getMetadata()
        visitInfo = self.makeRawVisitInfo(md=md, exposureId=exposureId)
        exposure.getInfo().setVisitInfo(visitInfo)
        return self._standardizeExposure(self.exposures['raw'], exposure, dataId,
                                         trimmed=False)

    def map_skypolicy(self, dataId):
        """Map a sky policy."""
        return dafPersist.ButlerLocation("lsst.pex.policy.Policy", "Policy",
                                         "Internal", None, None, self,
                                         storage=self.rootStorage)

    def std_skypolicy(self, item, dataId):
        """Standardize a sky policy by returning the one we use."""
        return self.skypolicy

###############################################################################
#
# Utility functions
#
###############################################################################

    def _getCcdKeyVal(self, dataId):
        """Return CCD key and value used to look a defect in the defect registry

        The default implementation simply returns ("ccd", full detector name)
        """
        return ("ccd", self._extractDetectorName(dataId))

    def _setupRegistry(self, name, path, policy, policyKey, storage, searchParents=True,
                       posixIfNoSql=True):
        """Set up a registry (usually SQLite3), trying a number of possible
        paths.

        Parameters
        ----------
        name : string
            Name of registry.
        path : string
            Path for registry.
        policy : string
            Policy that contains the registry name, used if path is None.
        policyKey : string
            Key in policy for registry path.
        storage : Storage subclass
            Repository Storage to look in.
        searchParents : bool, optional
            True if the search for a registry should follow any Butler v1
            _parent symlinks.
        posixIfNoSql : bool, optional
            If an sqlite registry is not found, will create a posix registry if
            this is True.

        Returns
        -------
        lsst.daf.persistence.Registry
            Registry object
        """
        if path is None and policyKey in policy:
            path = dafPersist.LogicalLocation(policy[policyKey]).locString()
            if os.path.isabs(path):
                raise RuntimeError("Policy should not indicate an absolute path for registry.")
            if not storage.exists(path):
                newPath = storage.instanceSearch(path)

                newPath = newPath[0] if newPath is not None and len(newPath) else None
                if newPath is None:
                    self.log.warn("Unable to locate registry at policy path (also looked in root): %s",
                                  path)
                path = newPath
            else:
                self.log.warn("Unable to locate registry at policy path: %s", path)
                path = None

        # Old Butler API was to indicate the registry WITH the repo folder, New Butler expects the registry to
        # be in the repo folder. To support Old API, check to see if path starts with root, and if so, strip
        # root from path.
        root = storage.root
        if path and (path.startswith(root)):
            path = path[len(root + '/'):]

        # determine if there is an sqlite registry and if not, try the posix registry.
        registry = None

        if path is None:
            path = "%s.sqlite3" % name
            newPath = storage.instanceSearch(path)
            newPath = newPath[0] if newPath is not None and len(newPath) else None
            if newPath is None:
                self.log.info("Unable to locate %s registry in root: %s", name, path)
            path = newPath
        if path is None:
            path = os.path.join(".", "%s.sqlite3" % name)
            newPath = storage.instanceSearch(path)
            newPath = newPath[0] if newPath is not None and len(newPath) else None
            if newPath is None:
                self.log.info("Unable to locate %s registry in current dir: %s", name, path)
            path = newPath
        if path is not None:
            if not storage.exists(path):
                newPath = storage.instanceSearch(path)
                newPath = newPath[0] if newPath is not None and len(newPath) else None
                if newPath is not None:
                    path = newPath
            self.log.debug("Loading %s registry from %s", name, path)
            registry = dafPersist.Registry.create(storage.getLocalFile(path))
        elif not registry and posixIfNoSql:
            self.log.info("Loading Posix registry from %s", storage.root)
            registry = dafPersist.PosixRegistry(storage.root)

        return registry

    def _transformId(self, dataId):
        """Generate a standard ID dict from a camera-specific ID dict.

        Canonical keys include:
        - amp: amplifier name
        - ccd: CCD name (in LSST this is a combination of raft and sensor)
        The default implementation returns a copy of its input.

        @param dataId[in] (dict) Dataset identifier; this must not be modified
        @return (dict) Transformed dataset identifier"""

        return dataId.copy()

    def _mapActualToPath(self, template, actualId):
        """Convert a template path to an actual path, using the actual data
        identifier.  This implementation is usually sufficient but can be
        overridden by the subclass.
        @param template (string) Template path
        @param actualId (dict) Dataset identifier
        @return (string) Pathname"""

        try:
            transformedId = self._transformId(actualId)
            return template % transformedId
        except Exception as e:
            raise RuntimeError("Failed to format %r with data %r: %s" % (template, transformedId, e))

    @staticmethod
    def getShortCcdName(ccdName):
        """Convert a CCD name to a form useful as a filename

        The default implementation converts spaces to underscores.
        """
        return ccdName.replace(" ", "_")

    def _extractDetectorName(self, dataId):
        """Extract the detector (CCD) name from the dataset identifier.

        The name in question is the detector name used by lsst.afw.cameraGeom.

        @param dataId (dict) Dataset identifier
        @return (string) Detector name
        """
        raise NotImplementedError("No _extractDetectorName() function specified")

    def _extractAmpId(self, dataId):
        """Extract the amplifier identifer from a dataset identifier.

        @warning this is deprecated; DO NOT USE IT

        amplifier identifier has two parts: the detector name for the CCD
        containing the amplifier and index of the amplifier in the detector.
        @param dataId (dict) Dataset identifer
        @return (tuple) Amplifier identifier"""

        trDataId = self._transformId(dataId)
        return (trDataId["ccd"], int(trDataId['amp']))

    def _setAmpDetector(self, item, dataId, trimmed=True):
        """Set the detector object in an Exposure for an amplifier.
        Defects are also added to the Exposure based on the detector object.
        @param[in,out] item (lsst.afw.image.Exposure)
        @param dataId (dict) Dataset identifier
        @param trimmed (bool) Should detector be marked as trimmed? (ignored)"""

        return self._setCcdDetector(item=item, dataId=dataId, trimmed=trimmed)

    def _setCcdDetector(self, item, dataId, trimmed=True):
        """Set the detector object in an Exposure for a CCD.
        @param[in,out] item (lsst.afw.image.Exposure)
        @param dataId (dict) Dataset identifier
        @param trimmed (bool) Should detector be marked as trimmed? (ignored)"""

        detectorName = self._extractDetectorName(dataId)
        detector = self.camera[detectorName]
        item.setDetector(detector)

    def _setFilter(self, mapping, item, dataId):
        """Set the filter object in an Exposure.  If the Exposure had a FILTER
        keyword, this was already processed during load.  But if it didn't,
        use the filter from the registry.
        @param mapping (lsst.obs.base.Mapping)
        @param[in,out] item (lsst.afw.image.Exposure)
        @param dataId (dict) Dataset identifier"""

        if not (isinstance(item, afwImage.ExposureU) or isinstance(item, afwImage.ExposureI) or
                isinstance(item, afwImage.ExposureF) or isinstance(item, afwImage.ExposureD)):
            return

        actualId = mapping.need(['filter'], dataId)
        filterName = actualId['filter']
        if self.filters is not None and filterName in self.filters:
            filterName = self.filters[filterName]
        item.setFilter(afwImage.Filter(filterName))

    # Default standardization function for exposures
    def _standardizeExposure(self, mapping, item, dataId, filter=True,
                             trimmed=True):
        """Default standardization function for images.

        This sets the Detector from the camera geometry
        and optionally set the Fiter. In both cases this saves
        having to persist some data in each exposure (or image).

        @param mapping (lsst.obs.base.Mapping)
        @param[in,out] item image-like object; any of lsst.afw.image.Exposure,
                lsst.afw.image.DecoratedImage, lsst.afw.image.Image
                or lsst.afw.image.MaskedImage
        @param dataId (dict) Dataset identifier
        @param filter (bool) Set filter? Ignored if item is already an exposure
        @param trimmed (bool) Should detector be marked as trimmed?
        @return (lsst.afw.image.Exposure) the standardized Exposure"""
        if not hasattr(item, "getMaskedImage"):
            try:
                item = exposureFromImage(item)
            except Exception as e:
                self.log.error("Could not turn item=%r into an exposure: %s" % (repr(item), e))
                raise

        if mapping.level.lower() == "amp":
            self._setAmpDetector(item, dataId, trimmed)
        elif mapping.level.lower() == "ccd":
            self._setCcdDetector(item, dataId, trimmed)

        if filter:
            self._setFilter(mapping, item, dataId)

        return item

    def _defectLookup(self, dataId):
        """Find the defects for a given CCD.
        @param dataId (dict) Dataset identifier
        @return (string) path to the defects file or None if not available"""
        if self.defectRegistry is None:
            return None
        if self.registry is None:
            raise RuntimeError("No registry for defect lookup")

        ccdKey, ccdVal = self._getCcdKeyVal(dataId)

        dataIdForLookup = {'visit': dataId['visit']}
        # .lookup will fail in a posix registry because there is no template to provide.
        rows = self.registry.lookup(('taiObs'), ('raw_visit'), dataIdForLookup)
        if len(rows) == 0:
            return None
        assert len(rows) == 1
        taiObs = rows[0][0]

        # Lookup the defects for this CCD serial number that are valid at the exposure midpoint.
        rows = self.defectRegistry.executeQuery(("path",), ("defect",),
                                                [(ccdKey, "?")],
                                                ("DATETIME(?)", "DATETIME(validStart)", "DATETIME(validEnd)"),
                                                (ccdVal, taiObs))
        if not rows or len(rows) == 0:
            return None
        if len(rows) == 1:
            return os.path.join(self.defectPath, rows[0][0])
        else:
            raise RuntimeError("Querying for defects (%s, %s) returns %d files: %s" %
                               (ccdVal, taiObs, len(rows), ", ".join([_[0] for _ in rows])))

    def _makeCamera(self, policy, repositoryDir):
        """Make a camera (instance of lsst.afw.cameraGeom.Camera) describing the camera geometry

        Also set self.cameraDataLocation, if relevant (else it can be left None).

        This implementation assumes that policy contains an entry "camera" that points to the
        subdirectory in this package of camera data; specifically, that subdirectory must contain:
        - a file named `camera.py` that contains persisted camera config
        - ampInfo table FITS files, as required by lsst.afw.cameraGeom.makeCameraFromPath

        @param policy        (daf_persistence.Policy, or pexPolicy.Policy (only for backward compatibility))
                             Policy with per-camera defaults already merged
        @param repositoryDir (string) Policy repository for the subclassing
                             module (obtained with getRepositoryPath() on the
                             per-camera default dictionary)
        """
        if isinstance(policy, pexPolicy.Policy):
            policy = dafPersist.Policy(pexPolicy=policy)
        if 'camera' not in policy:
            raise RuntimeError("Cannot find 'camera' in policy; cannot construct a camera")
        cameraDataSubdir = policy['camera']
        self.cameraDataLocation = os.path.normpath(
            os.path.join(repositoryDir, cameraDataSubdir, "camera.py"))
        cameraConfig = afwCameraGeom.CameraConfig()
        cameraConfig.load(self.cameraDataLocation)
        ampInfoPath = os.path.dirname(self.cameraDataLocation)
        return afwCameraGeom.makeCameraFromPath(
            cameraConfig=cameraConfig,
            ampInfoPath=ampInfoPath,
            shortNameFunc=self.getShortCcdName
        )

    def getRegistry(self):
        """Get the registry used by this mapper.

        Returns
        -------
        Registry or None
            The registry used by this mapper for this mapper's repository.
        """
        return self.registry

def exposureFromImage(image):
    """Generate an Exposure from an image-like object

    If the image is a DecoratedImage then also set its WCS and metadata
    (Image and MaskedImage are missing the necessary metadata
    and Exposure already has those set)

    @param[in] image  Image-like object (lsst.afw.image.DecoratedImage, Image, MaskedImage or Exposure)
    @return (lsst.afw.image.Exposure) Exposure containing input image
    """
    if hasattr(image, "getVariance"):
        # MaskedImage
        exposure = afwImage.makeExposure(image)
    elif hasattr(image, "getImage"):
        # DecoratedImage
        exposure = afwImage.makeExposure(afwImage.makeMaskedImage(image.getImage()))
        metadata = image.getMetadata()
        wcs = afwImage.makeWcs(metadata, True)
        exposure.setWcs(wcs)
        exposure.setMetadata(metadata)
    elif hasattr(image, "getMaskedImage"):
        # Exposure
        exposure = image
    else:
        # Image
        exposure = afwImage.makeExposure(afwImage.makeMaskedImage(image))

    return exposure

