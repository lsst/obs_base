from __future__ import print_function, division, absolute_import

import sqlite3
import datetime
from . import base
from . import graph


class Backend(object):

    def __init__(self, config=None):
        self.config = {"idcol": "id"}
        if config is not None:
            self.config.update(config)

    def getTableName(self, UnitClass):
        return UnitClass.__name__

    def createTable(self, UnitClass):
        raise NotImplementedError()

    def makeGraph(self, UnitClasses, DatasetClasses):
        raise NotImplementedError()


class SqliteBackend(Backend):

    SQL_TYPES = {
        base.RegionField: "BLOB",
        base.IntField: "INTEGER",
        base.StrField: "TEXT",
        base.DateTimeField: "INTEGER",
        base.ForeignKey: "INTEGER",
    }

    CONVERTERS = {
        base.DateTimeField: datetime.datetime.fromtimestamp,
    }

    def __init__(self, dbname, config=None):
        Backend.__init__(self, config)
        self.db = sqlite3.connect(dbname)
        self.db.row_factory = sqlite3.Row

    def createTable(self, UnitClass):
        sql = self._buildCreateTable(UnitClass)
        self.db.execute(sql)
        self.db.commit()

    def makeGraph(self, UnitClasses, DatasetClasses):
        sql = self._buildGraphQuery(UnitClasses)
        cursor = self.db.execute(sql)
        units = self._makeGraphUnits(UnitClasses, cursor)
        # Convert nested dicts (indexed by id) to frozensets: after we're done
        # constructing Units, we have no more need for indexing by ID, and
        # we want the set of all Units in a dict to be immutable.
        units = {k: frozenset(v.values()) for k, v in units.items()}
        datasets = {k: set() for k in DatasetClasses}
        return graph.RepoGraph(units=units, datasets=datasets)

    def _buildCreateTable(self, UnitClass):
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

    def _buildGraphQuery(self, UnitClasses):
        columns = []
        tables = []
        joins = []
        for UnitClass in UnitClasses:
            table = self.getTableName(UnitClass)
            tables.append(table)
            columns.append(
                "{}.{} AS {}_id".format(table, self.config["idcol"], table)
            )
            for f in UnitClass.fields.values():
                columns.append(
                    "{}.{} AS {}_{}".format(table, f.name, table, f.name)
                )
                if isinstance(f, base.ForeignKey):
                    joins.append(
                        "({current}.{column} = {foreign}.{idcol})".format(
                            current=table,
                            column=f.name,
                            foreign=self.getTableName(f.UnitClass),
                            idcol=self.config["idcol"]
                        )
                    )
        sql = "SELECT {} FROM {} WHERE ({})".format(
            ", ".join(columns),
            ", ".join(tables),
            " AND ".join(joins)
        )
        return sql

    def _makeGraphUnits(self, UnitClasses, cursor):
        # Precomputes some SQL-specific per-Field quantities so we don't
        # have to recompute them for every row in the query results.
        helpers = {}

        def noop(value):
            return value

        for UnitClass in UnitClasses:
            table = self.getTableName(UnitClass)
            helpers[UnitClass] = {
                "id": ("{}_{}".format(table, self.config["idcol"]), noop)
            }
            for f in UnitClass.fields.values():
                qname = "{}_{}".format(table, f.name)
                helpers[UnitClass][f.name] = (
                    qname,
                    self.CONVERTERS.get(UnitClass, noop)
                )

        # Iterate over the query results, constructing new Unit instances.
        units = {}
        row = cursor.fetchone()
        while row is not None:
            for UnitClass, helper in helpers.items():
                instances = units.setdefault(UnitClass, {})
                try:
                    id = row[helper["id"][0]]
                except IndexError:
                    print(helper["id"][0], row.keys())
                    raise
                if id in instances:
                    continue
                storage = {}
                for fname, (qname, converter) in helper.iteritems():
                    storage[fname] = converter(row[qname])
                instances[id] = UnitClass(storage)
            row = cursor.fetchone()
        # Iterate over newly created Unit instances, turning integer ForeignKey
        # values into objects and populating ReverseForeignKeys.
        for UnitClass, instances in units.items():
            for instance in instances.values():
                for f in UnitClass.fields.values():
                    f.finalize(instance, units)
        # TODO: convert ReverseForeignKey sets into frozensets to guarantee
        # immutability.
        return units
