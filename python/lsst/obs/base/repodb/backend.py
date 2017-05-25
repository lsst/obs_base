from __future__ import print_function, division, absolute_import

import sqlite3
import datetime
import itertools
from . import base
from . import graph


class Backend(object):

    def __init__(self, config=None):
        self.config = {"idcol": "id"}
        if config is not None:
            self.config.update(config)

    def getUnitTableName(self, UnitClass):
        return UnitClass.__name__

    def getDatasetTableName(self, DatasetClass):
        return DatasetClass.name

    def createUnitTable(self, UnitClass):
        raise NotImplementedError()

    def createDatasetTable(self, DatasetClass):
        raise NotImplementedError()

    def makeGraph(self, UnitClasses, where=None,
                  NeededDatasets=(), FutureDatasets=()):
        raise NotImplementedError()


def noop(value):
    return value


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

    def createUnitTable(self, UnitClass):
        items = ["{} INTEGER PRIMARY KEY".format(self.config["idcol"])]
        for k, v in UnitClass.fields.items():
            t = [k, self.SQL_TYPES[type(v)]]
            if isinstance(v, base.ForeignKey):
                t.append(
                    "REFERENCES {} ({})".format(
                        self.getUnitTableName(v.UnitClass),
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
        sql = "CREATE TABLE {} (\n    {}\n)".format(
            self.getUnitTableName(UnitClass), ",\n    ".join(items)
        )
        self.db.execute(sql)
        self.db.commit()

    def createDatasetTable(self, DatasetClass):
        items = []
        for k, v in DatasetClass.properties.items():
            t = [k, self.SQL_TYPES[base.ForeignKey]]
            t.append(
                "REFERENCES {} ({})".format(
                    self.getUnitTableName(v.UnitClass),
                    self.config["idcol"]
                )
            )
            t.append("NOT NULL")
            items.append(" ".join(t))
        items.append(
            "UNIQUE ({})".format(
                ", ".join(k for k in DatasetClass.properties.keys())
            )
        )
        sql = "CREATE TABLE {} (\n    {}\n)".format(
            self.getDatasetTableName(DatasetClass), ",\n    ".join(items)
        )
        self.db.execute(sql)
        self.db.commit()

    def makeGraph(self, UnitClasses, where=None,
                  NeededDatasets=(), FutureDatasets=()):
        temp = self._makeTemporaryGraphTable(UnitClasses, NeededDatasets,
                                             where=where)
        units = {}
        for UnitClass in UnitClasses:
            units[UnitClass] = self._readPartialUnits(UnitClass, temp)
        datasets = {}
        for DatasetClass in itertools.chain(NeededDatasets, FutureDatasets):
            datasets[DatasetClass] = self._readPartialDatasets(
                DatasetClass,
                temp
            )
        self.db.execute("DROP TABLE {}".format(temp))
        datasets = self._finalizeDatasets(datasets, units)
        units = self._finalizeUnits(units)
        return graph.RepoGraph(units=units, datasets=datasets)

    def _makeTemporaryGraphTable(self, UnitClasses, NeededDatasets,
                                 where=None, temp="graph"):
        columns = []
        tables = []
        if where is None:
            where = []
        else:
            assert where
            where = ["({})".format(where)]
        for UnitClass in UnitClasses:
            table = self.getUnitTableName(UnitClass)
            tables.append(table)
            columns.append(
                "{}.{} AS {}_id".format(table, self.config["idcol"], table)
            )
            for f in UnitClass.fields.values():
                if isinstance(f, base.ForeignKey):
                    where.append(
                        "({current}.{column} = {foreign}.{idcol})".format(
                            current=table,
                            column=f.name,
                            foreign=self.getUnitTableName(f.UnitClass),
                            idcol=self.config["idcol"]
                        )
                    )
        for DatasetClass in NeededDatasets:
            table = self.getDatasetTableName(DatasetClass)
            tables.append(table)
            for p in DatasetClass.properties.values():
                where.append(
                    "({current}.{column} = {foreign}.{idcol})".format(
                        current=table,
                        column=p.name,
                        foreign=self.getUnitTableName(p.UnitClass),
                        idcol=self.config["idcol"]
                    )
                )
        sql = ("CREATE TEMPORARY TABLE {temp} AS "
               "SELECT {columns} FROM {tables} WHERE ({where})").format(
            temp=temp,
            columns=", ".join(columns),
            tables=", ".join(tables),
            where=" AND ".join(where)
        )
        self.db.execute(sql)
        return temp

    def _readPartialUnits(self, UnitClass, temp):
        sql = ("SELECT DISTINCT {idcol}, {columns} FROM {table} "
               "INNER JOIN {temp} ON ({table}.{idcol} = {temp}.{table}_id)"
               ).format(
                    columns=", ".join(UnitClass.fields.keys()),
                    table=self.getUnitTableName(UnitClass),
                    idcol=self.config["idcol"],
                    temp=temp
                )
        converters = [self.CONVERTERS.get(type(f), noop)
                      for f in UnitClass.fields.values()]
        result = {}
        for row in self.db.execute(sql):
            id = row["id"]
            assert id not in result
            storage = {"id": id}
            for f, converter in zip(UnitClass.fields.values(), converters):
                storage[f.name] = converter(row[f.name])
            result[id] = UnitClass(storage)
        return result

    def _finalizeUnits(self, units):
        # Iterate over newly created Unit instances, turning integer ForeignKey
        # values into objects and populating ReverseForeignKeys.
        for UnitClass, instances in units.items():
            for instance in instances.values():
                for f in UnitClass.fields.values():
                    f.finalize(instance, units)
        # Convert from dict to frozenset: want the set of Units to
        # be immutable, and we don't need to look up by ID after
        # dereferencing all the foreign keys.
        return {k: frozenset(v.values()) for k, v in units.items()}

    def _readPartialDatasets(self, DatasetClass, temp):
        columns = []
        joins = []
        table = self.getDatasetTableName(DatasetClass)
        for p in DatasetClass.properties.values():
            columns.append(p.name)
            joins.append(
                "({table}.{field} = {temp}.{unit}_id)".format(
                    table=table,
                    field=p.name,
                    temp=temp,
                    unit=self.getUnitTableName(p.UnitClass)
                )
            )
        sql = ("SELECT DISTINCT {columns} FROM {table}, {temp} "
               "WHERE ({joins})").format(
                    columns=", ".join(DatasetClass.properties.keys()),
                    table=table,
                    temp=temp,
                    joins=" AND ".join(joins)
                )
        result = set()
        for row in self.db.execute(sql):
            storage = {}
            for p in DatasetClass.properties.values():
                storage[p.name] = row[p.name]
            result.add(DatasetClass(storage))
        return result

    def _finalizeDatasets(self, datasets, units):
        # Iterate over newly created Dataset instances, turning integer
        # Unit IDs into Unit instances and adding references from the Units
        # back to these Datasets.
        for DatasetClass, instances in datasets.items():
            for instance in instances:
                for p in DatasetClass.properties.values():
                    p.finalize(instance, units)
        return datasets
