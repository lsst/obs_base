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
import os
import re
try:
    import sqlite as sqlite3
    haveSqlite3 = True
except ImportError:
    try:
        import sqlite3
        haveSqlite3 = True
    except ImportError:
        haveSqlite3 = False

try:
    import psycopg2 as pgsql
    havePgSql = True
except ImportError:
    try:
        from pg8000 import DBAPI as pgsql
        havePgSql = True
    except ImportError:
        havePgSql = False

from lsst.pex.config import Config, Field
import lsst.pex.exceptions as pexExcept
from lsst.pex.config import Config, Field

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
        if re.match(r'pgsql:', location):
            return PgSqlRegistry(location)
        raise RuntimeError, \
                "Unable to create registry using location: " + location

    # TODO - flesh out the generic interface as more registry types are
    # defined.

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

#        print "SQL: "+ cmd
#        for s in values:
#            print s
        
        c = self.conn.execute(cmd, values)
        result = []
        for row in c:
            result.append(row)
#            print result
        return result

class PgSqlConfig(Config):
    host = Field(dtype=str, default="localhost", doc="hostname or IP address for database")
    port = Field(dtype=int, default=5432, doc="port for database")
    db = Field(dtype=str, doc="name of database for registry")
    user = Field(dtype=str, doc="username for database")
    password = Field(dtype=str, optional=True, doc="password for database")

class PgSqlRegistry(Registry):
    """A PostgreSQL-based registry."""

    def __init__(self, pgsqlConf):
        """Constructor.
        @param location (string) URL containing the PostgreSQL host, port, and database, of the form
                                 'pgsql://<host>:<port>/<database>'
        """
        if not havePgSql:
            raise RuntimeError("Cannot use PgSqlRegistry: could not import psycopg2 or pg8000")
        Registry.__init__(self)
        self.db = pgsql.connect(host=pgsqlConf.host, port=pgsqlConf.port, database=pgsqlConf.db,
                                user=pgsqlConf.user, password=pgsqlConf.password)
        self.conn = self.db.cursor()

    def __del__(self):
        if self.conn:
            self.conn.close()   # close cursor
        if self.db:
            self.db.close()     # close DB connection

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
        i = 0
        if whereFields: 
            for (k, v), val in zip(whereFields, values):
                whereList.append("(%s = '%s')" % (k, val))
                i += 1
                
        if range is not None:
            range_pgsql = (values[i], 'validStart', 'validEnd')
#            print range[0],range[1],range[2]
            whereList.append("(select cast('%s' as date) BETWEEN cast(%s as date) AND cast(%s as date))" 
                             % range_pgsql)
        if len(whereList) > 0:
            cmd += " WHERE " + " AND ".join(whereList)

#        print "SQL: "+ cmd

#        for (k, v), val in zip(whereFields, values):
#            print "##########"
#            print k, v, val
            
        self.conn.execute(cmd)
        c = self.conn.fetchall()
        result = []
        for row in c:
            new_row=[]
            for val in row:
                new_val=''
                if isinstance(val,unicode):
                    new_val = str(val)
                else:
                    new_val = val

                new_row.append(new_val)
            result.append(new_row)
#            print result
        return result

