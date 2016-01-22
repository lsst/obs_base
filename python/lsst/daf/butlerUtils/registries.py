#!/usr/bin/env python

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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#

"""This module provides registry classes for maintaining dataset metadata for
use by the Data Butler.  Currently only a SQLite3-based registry is
implemented, but registries based on a text file, a policy file, a MySQL (or
other) relational database, and data gathered from scanning a filesystem are
all anticipated."""

import copy
import fsScanner
import os
import pyfits
import re
try:
    import sqlite3
    haveSqlite3 = True
except ImportError:
    try:
        # try external pysqlite package; deprecated
        import sqlite as sqlite3
        haveSqlite3 = True
    except ImportError:
        haveSqlite3 = False

class Registry(object):
    """The registry base class."""

    def __init__(self):
        pass

    @staticmethod
    def create(location):
        """Create a registry object of an appropriate type.
        @param location (string) Path or URL for registry, or None if
                                 unavailable"""

        # if re.match(r'.*\.registry', location):
        #     return FileRegistry(location)
        # if re.match(r'.*\.paf', location):
        #     return CalibRegistry(location)
        if haveSqlite3 and re.match(r'.*\.sqlite3', location):
            registry = SqliteRegistry(location)
            if registry.conn is None:
                return None
            return registry
        # if re.match(r'mysql:', location):
        #     return DbRegistry(location)
        # return FsRegistry(location)
        raise RuntimeError, \
                "Unable to create registry using location: " + location

class PosixRegistry(Registry):
    """A glob-based filesystem registry"""

    def __init__(self, root):
        Registry.__init__(self)
        self.root = root

    @staticmethod
    def getHduNumber(template, dataId):
        """Looks up the HDU number for a given template+dataId.
        :param template: template with HDU specifier (ends with brackets and an
        identifier that can be populated by a key-value pair in dataId.
        e.g. "%(visit)07d/instcal%(visit)07d.fits.fz[%(ccdnum)d]"
        :param dataId: dictionary that hopefully has a key-value pair whose key
        matches (has the same name) as the key specifier in the template.
        :return: the HDU specified by the template+dataId pair, or None if the
        HDU can not be determined.
        """
        #sanity check that the template at least ends with a brace.
        if not template.endswith(']'):
            return None

        # get the key (with formatting) out of the brances
        hduKey = template[template.rfind('[')+1:template.rfind(']')]
        # extract the key name from the formatting
        hduKey = hduKey[hduKey.rfind('(')+1:hduKey.rfind(')')]

        if hduKey in dataId:
            return dataId[hduKey]
        return None

    class LookupData:
        def __init__(self, lookupProperties, dataId):
            self.dataId = copy.copy(dataId)
            if not hasattr(lookupProperties, '__iter__'):
                lookupProperties = (lookupProperties,)
            self.lookupProperties = copy.copy(lookupProperties)
            self.foundItems = {}
            self.cachedStatus = None
            self.neededKeys = set(lookupProperties).union(dataId.keys())

        def __repr__(self):
            return "LookupData lookupProperties:%s dataId:%s foundItems:%s cachedStatus:%s" % \
                   (self.lookupProperties, self.dataId, self.foundItems, self.cachedStatus)

        def status(self):
            """Query the lookup status

            :return: 'match' if the key+value pairs in dataId have been satisifed and keys in lookupProperties have
            found and their key+value added to resolvedId
            'incomplete' if the found data matches but there are still incomplete data matching in dataId or
            lookupProperties
            'not match' if data in foundId does not match data in dataId
            """
            if self.cachedStatus is not None:
                return self.cachedStatus
            self.cachedStatus = 'match'
            for key in self.dataId:
                if key not in self.foundItems:
                    self.cachedStatus = 'incomplete'
                elif self.dataId[key] != self.foundItems[key]:
                    self.cachedStatus = 'incomplete'
            if self.cachedStatus == 'match':
                for key in self.lookupProperties:
                    if key not in self.foundItems:
                        self.cachedStatus = 'incomplete'
                        break
            return self.cachedStatus

        def setFoundItems(self, items):
            self.cachedStatus = None
            self.foundItems = items

        def addFoundItems(self, items):
            self.cachedStatus = None
            self.foundItems.update(items)

        def getMissingKeys(self):
            return self.neededKeys - set(self.foundItems.keys())

    def lookup(self, lookupProperties, reference, dataId, **kwargs):
        """Perform a lookup in the registry.

        Return values are refined by the values in dataId.
        Returns a list of values that match keys in lookupProperties.
        e.g. if the template is 'raw/raw_v%(visit)d_f%(filter)s.fits.gz', and
        dataId={'visit':1}, and lookupProperties is ['filter'], and the
        filesystem under self.root has exactly one file 'raw/raw_v1_fg.fits.gz'
        then the return value will be [('g',)]

        :param lookupProperties: keys whose values will be returned.
        :param reference: other data types that may be used to search for values.
        :param dataId: must be an iterable. Keys must be string.
        If value is a string then will look for elements in the repository that match value for key.
        If value is a 2-item iterable then will look for elements in the repository are between (inclusive)
        the first and second items in the value.
        :param **kwargs: keys required for the posix registry to search for items. If required keys are not
        provide will return an empty list.
        'template': required. template parameter (typically from a policy) that can be used to look for files
        'storage': optional. Needed to look for metadata in files. Currently supported values: 'FitsStorage'.
        :return: a list of values that match keys in lookupProperties.
        """
        # required kwargs:
        if 'template' in kwargs:
            template = kwargs['template']
        else:
            return []
        # optional kwargs:
        storage = kwargs['storage'] if 'storage' in kwargs else None

        lookupData = PosixRegistry.LookupData(lookupProperties, dataId)
        scanner = fsScanner.FsScanner(template)
        allPaths = scanner.processPath(self.root)
        retItems = [] # one item for each found file that matches
        for path, foundProperties in allPaths.items():
            # check for dataId keys that are not present in found properties
            # search for those keys in metadata of file at path
            # if present, check for matching values
            # if not present, file can not match, do not use it.
            lookupData.setFoundItems(foundProperties)
            if 'incomplete' == lookupData.status():
                PosixRegistry.lookupMetadata(os.path.join(self.root, path), template, lookupData, storage)
            if 'match' == lookupData.status():
                l = tuple(lookupData.foundItems[key] for key in lookupData.lookupProperties)
                retItems.append(l)
        return retItems

    @staticmethod
    def lookupMetadata(filepath, template, lookupData, storage):
        """Dispatcher for looking up metadata in a file of a given storage type
        """
        if storage == 'FitsStorage':
            PosixRegistry.lookupFitsMetadata(filepath, template, lookupData, storage)

    @staticmethod
    def lookupFitsMetadata(filepath, template, lookupData, dataId):
        """Look up metadata in a fits file.
        Will try to discover the correct HDU to look in by testing if the
        template has a value in brackets at the end.
        If the HDU is specified but the metadata key is not discovered in
        that HDU, will look in the primary HDU before giving up.
        :param filepath: path to the file
        :param template: template that was used to discover the file. This can
        be used to look up the correct HDU as needed.
        :param lookupData: an instance if LookupData that contains the
        lookupProperties, the dataId, and the data that has been found so far.
        Will be updated with new information as discovered.
        :param dataId:
        :return:
        """
        try:
            hdulist = pyfits.open(filepath, memmap=True)
        except IOError:
            return;
        hduNumber = PosixRegistry.getHduNumber(template=template, dataId=dataId)
        if hduNumber != None and hduNumber < len(hdulist):
            hdu = hdulist[hduNumber]
        else:
            hdu = None
        if len(hdulist) > 0:
            primaryHdu = hdulist[0]
        else:
            primaryHdu = None

        for property in lookupData.getMissingKeys():
            propertyValue = None
            if hdu is not None and property in hdu.header:
                propertyValue = hdu.header[property]
            # if the value is not in the indicated HDU, try the primary HDU:
            elif primaryHdu is not None and property in primaryHdu.header:
                propertyValue = primaryHdu.header[property]
            lookupData.addFoundItems({property: propertyValue})


class SqliteRegistry(Registry):
    """A SQLite3-based registry."""

    def __init__(self, location):
        """Constructor.
        @param location (string) Path to SQLite3 file"""

        Registry.__init__(self)
        if os.path.exists(location):
            self.conn = sqlite3.connect(location)
            self.conn.text_factory = str
        else:
            self.conn = None

    def lookup(self, lookupProperties, reference, dataId, **kwargs):
        """Perform a lookup in the registry.

        Return values are refined by the values in dataId.
        Returns a list of values that match keys in lookupProperties.
        e.g. if the template is 'raw/raw_v%(visit)d_f%(filter)s.fits.gz', and
        dataId={'visit':1}, and lookupProperties is ['filter'], and the
        filesystem under self.root has exactly one file 'raw/raw_v1_fg.fits.gz'
        then the return value will be [('g',)]

        :param lookupProperties:
        :param dataId: must be an iterable. Keys must be string.
        If key is a string then will look for elements in the repository that match value for key.
        If key is a 2-item iterable then will look for elements in the repository where the value is between
        the values of key[0] and key[1].
        :param reference: other data types that may be used to search for values.
        :param **kwargs: nothing needed for sqlite lookup
        :return: a list of values that match keys in lookupProperties.
        """
        if not self.conn:
            return None

        # input variable sanitization:
        if not hasattr(reference, '__iter__'):
            reference = (reference,)
        if not hasattr(lookupProperties, '__iter__'):
            lookupProperties = (lookupProperties,)

        cmd = "SELECT DISTINCT "
        cmd += ", ".join(lookupProperties)
        cmd += " FROM " + " NATURAL JOIN ".join(reference)
        valueList = []
        if dataId is not None and len(dataId) > 0:
            whereList = []
            for k, v in dataId.items():
                if hasattr(k, '__iter__'):
                    if len(k) is not 2:
                        raise RuntimeError("Wrong number of keys for range:%s" %k)
                    whereList.append("(? BETWEEN %s AND %s)" %(k[0], k[1]))
                    valueList.append(v)
                else:
                    whereList.append("%s = ?" %k)
                    valueList.append(v)
            cmd += " WHERE " + " AND ".join(whereList)
        c = self.conn.execute(cmd, valueList)
        result = []
        for row in c:
            result.append(row)
        return result

    def executeQuery(self, returnFields, joinClause, whereFields, range, values):
        """Extract metadata from the registry.
        @param returnFields (list of strings) Metadata fields to be extracted.
        @param joinClause   (list of strings) Tables in which metadata fields
                            are located.
        @param whereFields  (list of tuples) First tuple element is metadata
                            field to query; second is the value that field
                            must have (often '?').
        @param range        (tuple) Value, lower limit, and upper limit for a
                            range condition on the metadata.  Any of these can
                            be metadata fields.
        @param values       (tuple) Tuple of values to be substituted for '?'
                            characters in the whereFields values or the range
                            values.
        @return (list of tuples) All sets of field values that meet the
                criteria"""
        if not self.conn:
            return None
        cmd = "SELECT DISTINCT "
        cmd += ", ".join(returnFields)
        cmd += " FROM " + " NATURAL JOIN ".join(joinClause)
        whereList = []
        if whereFields: 
            for k, v in whereFields:
                whereList.append("(%s = %s)" % (k, v))
        if range is not None:
            whereList.append("(%s BETWEEN %s AND %s)" % range)
        if len(whereList) > 0:
            cmd += " WHERE " + " AND ".join(whereList)
        c = self.conn.execute(cmd, values)
        result = []
        for row in c:
            result.append(row)
        return result