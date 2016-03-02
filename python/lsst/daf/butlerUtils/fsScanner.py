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


"""This module provides the FsScanner class."""

import glob
import os
import re
import sys

class FsScanner(object):
    """Class to scan a filesystem location for paths matching a template.

    Decomposes the resulting paths into fields and passes them to a callback
    function.
    """

    def __init__(self, pathTemplate):
        """Constructor.  Takes the path template, which should be in the form
        of a Python string with named format substitution specifications.
        Such a template would be suitable for generating a path given a set of
        fields in a dictionary.  Does not handle hex (%x or %X).
        
        Example:
            %(field)s/%(visit)d/%(exposure)d/raw-%(visit)d-e%(exposure)03d-c%(ccd)03d-a%(amp)03d.fits
            
        Note that fields may appear multiple times; the second and subsequent
        appearances of such fields will have "_{number}" appended to them to
        disambiguate, although it is typically assumed that they will all be
        identical.

        Trailing brackets (and their contents) can be used to indicate which HDU from a file should
        be used. They will not be included in the filename search.
        """

        # Trim any trailing braces off the end of the path template.
        if pathTemplate.endswith(']'):
            pathTemplate = pathTemplate[0:pathTemplate.rfind('[')]

        # Change template into a globbable path specification.
        fmt = re.compile(r'%\((\w+)\).*?([dioueEfFgGcrs])')

        self.globString = fmt.sub('*', pathTemplate)

        # Change template into a regular expression.
        last = 0
        self.fields = {}
        self.reString = ""
        n = 0
        pos = 0
        for m in fmt.finditer(pathTemplate):
            fieldName = m.group(1)
            if fieldName in self.fields:
                fieldName += "_%d" % (n,)
                n += 1

            prefix = pathTemplate[last:m.start(0)]
            last = m.end(0)
            self.reString += prefix
        
            if m.group(2) in 'crs':
                fieldType = str
                self.reString += r'(?P<' + fieldName + '>.+)'
            elif m.group(2) in 'eEfFgG':
                fieldType = float
                self.reString += r'(?P<' + fieldName + '>[\d.eE+-]+)'
            else:
                fieldType = int
                self.reString += r'(?P<' + fieldName + '>[\d+-]+)'

            self.fields[fieldName] = dict(pos=pos, fieldType=fieldType)
            pos += 1

        self.reString += pathTemplate[last:] 

    def getFields(self):
        """Return the list of fields that will be returned from matched
        paths, in order."""

        fieldList = ["" for i in xrange(len(self.fields))]
        for f in self.fields.keys():
            fieldList[self.fields[f]['pos']] = f
        return fieldList

    def isNumeric(self, name):
        """Return true if the given field contains a number."""

        return self.fields[name]['fieldType'] in (float, int)

    def isInt(self, name):
        """Return true if the given field contains an integer."""

        return self.fields[name]['fieldType'] == int

    def isFloat(self, name):
        """Return true if the given field contains an float."""

        return self.fields[name]['fieldType'] == float

    def processPath(self, location):
        """
        Scan a given path location. Return info about paths that conform to the path template:
        :param location:
        :return: Path info: {path: {key:value ...}, ...} e.g.:
            {'0239622/instcal0239622.fits.fz': {'visit_0': 239622, 'visit': 239622}}
        """
        ret = {}
        curdir = os.getcwd()
        os.chdir(location)
        pathList = glob.glob(self.globString)
        for path in pathList:
            m = re.search(self.reString, path)
            if m:
                dataId = m.groupdict()
                for f in self.fields.keys():
                    if self.isInt(f):
                        dataId[f] = int(dataId[f])
                    elif self.isFloat(f):
                        dataId[f] = float(dataId[f])
                ret[path] = dataId
            else:
                print >> sys.stderr, "Warning: unmatched path: %s" % (path,)
        os.chdir(curdir)
        return ret
