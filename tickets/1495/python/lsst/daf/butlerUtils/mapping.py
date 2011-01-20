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
from lsst.daf.persistence import ButlerLocation
import lsst.pex.policy as pexPolicy

"""This module defines the Mapping base class."""

class Mapping(object):

    """Mapping is a base class for all mappings.  Mappings are used by
    the Mapper to map (determine a path to some data given some
    identifiers) and standardize (convert data into some standard
    format or type) data, and to query the associated registry to see
    what data is available.

    Public methods: lookup, have, need

    Mappings are specified mainly by policy.  A Mapping policy should
    consist of:

    template (string): a Python string providing the filename for that
    particular data set type based on some data identifiers.  In the
    case of redundancy in the path (e.g., file uniquely specified by
    the exposure number, but filter in the path), the
    redundant/dependent identifiers can be looked up in the registry.

    python (string): the Python type for the data (e.g.,
    lsst.afw.image.ExposureF)

    cpp (string): the C++ type for the data (e.g., ExposureF)

    storage (string): the storage type for the data (e.g., FitsStorage)

    level (string): the level in the camera hierarchy at which the
    data is stored (Amp, Ccd or skyTile)

    In addition, the following optional entries are permitted:

    map (string): the name of a function in the appropriate Mapper
    subclass that will map a data identifier to a path.  The function
    will receive: the Mapper, this Mapping, and the data identifier
    dict.

    standardize (string): the name of a function in the appropriate
    Mapper subclass that will standardize the input data.  The
    function will receive: the Mapper, this Mapping, the item to
    standardize, and the data identifier dict.  The special value
    "None" means no standardization is performed.

    tables (string): a whitespace-delimited list of tables in the
    registry that can be NATURAL JOIN-ed to look up additional
    information.

    query (string): a Python string providing a SQL query to look up
    the registry.  It should typically start as "SELECT *".

    lookup (string): the name of a function in the appropriate Mapper
    subclass that will lookup the desired properties.  The function
    will receive: the Mapper, this Mapping, a list of the desired
    properties, and the data identifier dict.    
    """

    def __init__(self, mapper, policy, datasetType, registry=None, root=None):
        """Constructor for Mapping class.
        @param[in,out] mapper (lsst.daf.persistence.Mapper) Mapper object
        @param policy         (lsst.pex.policy.Policy) Mapping policy
        @param datasetType    (string)
        @param registry       (lsst.daf.butlerUtils.Registry) Registry for metadata lookups
        @param root           (string) Path of root directory"""

        if mapper is None:
            raise RuntimeError, "No mapper provided for mapping"
        if policy is None:
            raise RuntimeError, "No policy provided for mapping"

        self.datasetType = datasetType
        self.registry = registry
        self.root = root

        self.template = policy.getString("template") # Template path
        self.python = policy.getString("python") # Python type
        self.cpp = policy.getString("cpp") # C++ type
        self.storage = policy.getString("storage") # Storage type
        self.level = policy.getString("level") # Level in camera hierarchy
        if not hasattr(mapper, "map_" + datasetType):
            setattr(mapper, "map_" + datasetType,
                    lambda dataId: self._map(mapper, dataId))
        if not hasattr(mapper, "std_" + datasetType):
            setattr(mapper, "std_" + datasetType,
                    lambda item, dataId: mapper._standardize(
                        self, item, dataId))
        if policy.exists("tables"):
            self.tables = policy.getStringArray("tables")
        else:
            self.tables = None
        if policy.exists("match"):
            self.match = policy.getString("match")
        else:
            self.match = None
        self.range = None
        if not hasattr(mapper, "query_" + datasetType):
            setattr(mapper, "query_" + datasetType,
                    lambda key, format, dataId: self.lookup(
                        mapper, format, dataId))

    def _map(self, mapper, dataId):
        """Standard implementation of map function.
        @param mapper (lsst.daf.persistence.Mapper)
        @param dataId (dict) Dataset identifier
        @return (lsst.daf.persistence.ButlerLocation)"""

        properties = re.findall(r'\%\((\w+)\)', self.template)
        actualId = self.need(mapper, properties, dataId)
        path = os.path.join(self.root,
                mapper._mapActualToPath(self.template, actualId))
        return ButlerLocation(self.python, self.cpp, self.storage,
                path, actualId)

    def lookup(self, mapper, properties, dataId):
        """Look up properties for in a metadata registry given a partial
        dataset identifier.
        @param mapper     (lsst.daf.persistence.Mapper) needed for subclasses
        @param properties (list of strings)
        @param dataId     (dict) Dataset identifier
        @return (list of tuples) values of properties"""

        if self.registry is None:
            raise RuntimeError, "No registry for lookup"

        where = {}
        values = []
        if dataId is not None:
            for k, v in dataId.iteritems():
                if k == 'taiObs':
                    continue
                where[k] = '?'
                values.append(v)
            if self.range is not None:
                values.append(dataId['taiObs'])
        return self.registry.executeQuery(properties, self.tables, where,
                self.range, values)

    def have(self, properties, dataId):
        """Returns whether the provided data identifier has all
        the properties in the provided list.
        @param properties (list of strings) Properties required
        @parm dataId      (dict) Dataset identifier
        @return (bool) True if all properties are present"""
        for prop in properties:
            if not dataId.has_key(prop):
                return False
        return True

    def need(self, mapper, properties, dataId, refresh=False, clobber=False):
        """Ensures all properties in the provided list are present in
        the data identifier, looking them up as needed.  This is only
        possible for the case where the data identifies a single
        exposure.
        @param mapper     (lsst.daf.persistence.Mapper)
        @param properties (list of strings) Properties required
        @param dataId     (dict) Partial dataset identifier
        @param refresh    (bool) Refresh values if present?
        @param clobber    (bool) Clobber existing values?
        @return (dict) copy of dataset identifier with enhanced values
        """
        
        if dataId is None:
            newId = dict()
        else:
            newId = dataId.copy()
        if self.have(properties, newId):
            return newId
        if not refresh:
            newProps = []                    # Properties we don't already have
            for prop in properties:
                if not newId.has_key(prop):
                    newProps.append(prop)
            if len(newProps) == 0:
                return newId
        else:
            newProps = properties

        lookups = self.lookup(mapper, newProps, newId)
        if len(lookups) != 1:
            raise RuntimeError, "No unique lookup for %s from %s: %d matches" % (newProps, newId, len(lookups))
        for i, prop in enumerate(newProps):
            if clobber or not dataId.has_key(prop):
                newId[prop] = lookups[0][i]
        return newId


class CalibrationMapping(Mapping):
    """CalibrationMapping is a Mapping subclass for calibration-type products.

    The difference is that data properties in the query or template
    can be looked up using the "raw" Mapping in addition to the calibration's.
    """

    def __init__(self, mapper, policy, datasetType, registry=None, root=None):
        """Constructor for Mapping class.
        @param[in,out] mapper (lsst.daf.persistence.Mapper) Mapper object
        @param policy         (lsst.pex.policy.Policy) Mapping policy
        @param datasetType    (string)
        @param registry       (lsst.daf.butlerUtils.Registry) Registry for metadata lookups
        @param root           (string) Path of root directory"""
        Mapping.__init__(self, mapper, policy, datasetType, registry, root)
        if policy.exists("validRange") and policy.getBool("validRange"):
            self.range = ("DATE(?)", "DATE(validStart)", "DATE(validEnd)")

    def lookup(self, mapper, properties, dataId):
        """Look up properties for in a metadata registry given a partial
        dataset identifier.
        @param properties (list of strings)
        @param dataId     (dict) Dataset identifier
        @return (list of tuples) values of properties"""

        queryProps = ["taiObs"]
        if self.where is not None:
            queryProps.append(self.where)
        if dataId is None or not self.have(queryProps, dataId):
            # Try to get them from the mapping for the raw data
            raw = mapper.getMapping("raw")
            dataId = raw.need(mapper, queryProps, dataId,
                    refresh=False, clobber=False)
        return Mapping.lookup(self, mapper, properties, dataId)

    def need(self, mapper, properties, dataId, refresh=False, clobber=False):
        """Ensures all properties in the provided list are present in
        the data identifier, looking them up as needed.  This is only
        possible for the case where the data identifies a single
        exposure.
        @param mapper     (lsst.daf.persistence.Mapper)
        @param properties (list of strings) Properties required
        @param dataId     (dict) Partial dataset identifier
        @param refresh    (bool) Refresh values if present?
        @param clobber    (bool) Clobber existing values?
        @return (dict) copy of dataset identifier with enhanced values
        """
        # Always want to clobber anything extant, because it's a property of
        # the input data, not the calibration data.
        return Mapping.need(self, mapper, properties, dataId, refresh=False, clobber=True)
