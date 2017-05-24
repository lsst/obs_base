from __future__ import print_function, division, absolute_import

import sqlite3
from . import base


class Backend(object):

    def __init__(self, config=None):
        self.config = {"idcol": "id"}
        if config is not None:
            self.config.update(config)

    def getTableName(self, UnitClass):
        return UnitClass.__name__

    def createTable(self, UnitClass):
        raise NotImplementedError()


class SqliteBackend(Backend):

    SQL_TYPES = {
        base.RegionField: "BLOB",
        base.IntField: "INTEGER",
        base.StrField: "TEXT",
        base.DateTimeField: "INTEGER",
        base.ForeignKey: "INTEGER",
    }

    def __init__(self, dbname, config=None):
        Backend.__init__(self, config)
        self.db = sqlite3.connect(dbname)

    def createTable(self, UnitClass):
        self.db.execute(self.getCreateTableString(UnitClass))

    def commit(self):
        self.db.commit()

    def getCreateTableString(self, UnitClass):
        items = ["{} INTEGER PRIMARY KEY".format(self.config["idcol"])]
        for k, v in UnitClass.fields.items():
            t = [k, self.SQL_TYPES[type(v)]]
            if isinstance(v, base.ForeignKey):
                t.append(
                    "REFERENCES {} ({})".format(
                        self.getTableName(v.UnitClass),
                        self.config["idcol"]
                    )
                )
            if not v.optional:
                t.append("NOT NULL")
            items.append(" ".join(t))
        items.append(
            "UNIQUE ({})".format(
                ", ".join(f.name for f in UnitClass.unique)
            )
        )
        return "CREATE TABLE {} (\n    {}\n)".format(
            self.getTableName(UnitClass), ",\n    ".join(items)
        )
