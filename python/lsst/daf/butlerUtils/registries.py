#!/usr/bin/env python

# import glob
# import os
import re
try:
    import sqlite3
    haveSqlite3 = True
except:
    haveSqlite3 = False

class Registry(object):
    def __init__(self):
        pass

    @staticmethod
    def create(location):
        # if re.match(r'.*\.registry', location):
        #     return FileRegistry(location)
        # if re.match(r'.*\.paf', location):
        #     return CalibRegistry(location)
        if haveSqlite3 and re.match(r'.*\.sqlite3', location):
            return SqliteRegistry(location)
        # if re.match(r'mysql:', location):
        #     return DbRegistry(location)
        # return FsRegistry(location)
        return None

class SqliteRegistry(Registry):
    def __init__(self, location):
        Registry.__init__(self)
        if os.path.exists(location):
            self.conn = sqlite3.connect(location)
        else:
            self.conn = None

    def __del__(self):
        if self.conn:
            self.conn.close()

    def executeQuery(self, returnFields, joinClause, whereFields, range,
            values):
        if not self.conn:
            return None
        cmd = "SELECT DISTINCT "
        cmd += ", ".join(returnFields)
        cmd += " FROM " + " NATURAL JOIN ".join(joinClause)
        whereList = []
        if whereFields: 
            for k, v in whereFields.iteritems():
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
