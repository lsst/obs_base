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

# import glob
import fsScanner
import os
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


    def lookup(self, template, lookupProperties, dataId, storage=None):
        """Looks for files in self.root (file system path) that match the
        described template. Return values are refined by the values in dataId.
        Returns a list of values that match keys in lookupProperties.
        e.g. if the template is 'raw/raw_v%(visit)d_f%(filter)s.fits.gz', and
        dataId={'visit':1}, and lookupProperties is ['filter'], and the
        filesystem under self.root has exactly one file 'raw/raw_v1_fg.fits.gz'
        then the return value will be [('g',)]
        :param template: template parameter (typically from a policy paf) that
        can be used to look for files
        :param lookupProperties: property keys to look up
        :param dataId: property keys and values to match
        :param storage: optional. If specified, partial file name matches will
        look in metadata for more values.
        :return: list of tuples of metadata values that match the list of keys
        in lookupProperties.
        """
        def commonKeysMatch(dict1, dict2):
            """Compare 2 dictionaries. Check all keys in common; if they match
            return true, else return false. Keys that are not common to both
            dictionaries will not be considered."""
            matchingKeys = set(dict1.keys()) & set(dict2.keys())
            for key in matchingKeys:
                if dict1[key] != dict2[key]:
                    return False
            return True

        scanner = fsScanner.FsScanner(template)
        allPaths = scanner.processPath(self.root)
        retItems = [] # one of these for each found file that matches
        foundPropertyList = []
        for path, foundProperties in allPaths.items():
            if commonKeysMatch(dataId, foundProperties):
                if not all(k in foundProperties for k in lookupProperties):
                    mdProperties = self.lookupMetadata(path, template, lookupProperties, dataId, storage)
                    mdProperties.update(foundProperties)
                    foundProperties = mdProperties
                if not all(k in foundProperties for k in lookupProperties):
                    break # incomplete lookup, can't use it.
                for property in lookupProperties:
                    foundPropertyList.append(foundProperties[property])
                retItems.append(tuple(foundPropertyList))
        return retItems


    def lookupMetadata(self, filepath, template, lookupProperties, dataId, storage):
        """
        Dispatcher for looking up metadata in a file of a given storage type
        :param template: template parameter (typically from a policy paf) that
        can be used to look for files
        :param lookupProperties: property keys to look up
        :param dataId: property keys and values to match
        :param storage: optional. If specified, partial file name matches will
        look in metadata for more values.
        :return: list of tuples of metadata values that match the list of keys
        in lookupProperties. e.g.:
        """
        ret = []
        if storage == 'FitsStorage':
            ret = self.lookupFitsMetadata(filepath, template, lookupProperties, dataId)
        return ret


    def lookupFitsMetadata(self, filepath, template, lookupProperties, dataId):
        """Lookup metadata in a fits file.
        Will try to discover the correct HDU to look in by testing if the
        template has a value in brackets at the end.
        If the HDU is specified but the metadata key is not discovered in
        that HDU, will look in the primary HDU before giving up.
        :param filepath: path to the file
        :param template: template that was used to discover the file. This can
        be used to look up the correct HDU as needed.
        :param lookupProperties: a list keys to metadata to be looked up in the
        file.
        :param dataId: key+value pairs to look up the HDU number.
        :return: list of metadata values that matches the list of metadata keys
        in lookupProperties
        """
        from astropy.io import fits
        hdulist = fits.open(filepath, memmap=True)
        hduNumber = PosixRegistry.getHduNumber(template=template, dataId=dataId)
        if hduNumber != None and hduNumber < len(hdulist):
            hdu = hdulist[hduNumber]
        else:
            hdu = None
        if len(hdulist) > 0:
            primaryHdu = hdulist[0]
        else:
            primaryHdu = None
        foundProperties = {}
        for property in lookupProperties:
            propertyValue = None
            if hdu is not None and property in hdu.header:
                propertyValue = hdu.header[property]
            # if the value is not in the indicated HDU, try the primary HDU:
            elif primaryHdu is not None and property in primaryHdu.header:
                propertyValue = primaryHdu.header[property]
            foundProperties[property] = propertyValue
        return foundProperties


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

    def executeQuery(self, returnFields, joinClause, whereFields, range,
            values):
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