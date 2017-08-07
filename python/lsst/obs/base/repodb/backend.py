from __future__ import print_function, division, absolute_import

import sqlite3
import datetime
import itertools
from .unit import Unit
from . import fields
from . import graph


class Backend(object):
    """An abstract base class that sets the interface between RepoDatabase
    and a SQL database.

    A Backend implementation should not depend on any concrete `Unit` or
    `Dataset` classes, but it is expected to be aware of these base classes
    and the set of descriptor classes used to define their subclasses.

    Design Notes
    ------------
    The `Backend` interface is almost certainly incomplete, given the current
    limitations of `RepoDatabase`.  It's also probably not quite the right
    interface for the insulation layer it tries to provide, as its design has
    been so far by a single implementation (`SqliteBackend`).

    When its design is more complete, `Backend` *should* be the only interface
    that needs to be reimplemented to enable the `RepoDatabase` system to work
    with the monolithic production database.

    Parameters
    ----------
    config : `dict`
        A dict of configuration values shared by all backends.  At present the
        only entry is `idcol`, which provides the string name used for integer
        primary key fields in the database.

    """

    def __init__(self, config=None):
        self.config = {"idcol": "id"}
        if config is not None:
            self.config.update(config)

    def getUnitTableName(self, UnitClass):
        """Return the name to use for a SQL table or view containing instances
        for the given `Unit` type.
        """
        name = UnitClass.__name__
        while UnitClass != Unit:
            name = UnitClass.__name__
            UnitClass = UnitClass.__bases__[0]
        return name

    def getDatasetTableName(self, DatasetClass):
        """Return the name to use for a SQL table or view containing instances
        for the given `Dataset` type.
        """
        return DatasetClass.name

    def createUnitTable(self, UnitClass):
        """Create a table or view that can hold instances of the given `Unit`
        type.

        If the given unit type does not inherit from `Unit`, any intermediate
        base classes with one or more `Fields` must already have associated
        tables.
        """
        raise NotImplementedError()

    def createDatasetTable(self, DatasetClass):
        """Create a table or view that can hold instances of the given
        `Dataset` type.

        Must be a silent no-op if the table already exists.
        """
        raise NotImplementedError()

    def insertUnit(self, unit):
        """Insert a `Unit` instance into the appropriate table.
        """
        raise NotImplementedError()

    def insertDataset(self, dataset):
        """Insert a `Dataset` instance into the appropriate table.
        """
        raise NotImplementedError()

    def makeGraph(self, UnitClasses=(), where=None,
                  NeededDatasets=(), FutureDatasets=()):
        """Execute a multi-table SQL query and return the results
        as a `RepoGraph`.

        See `RepoDatabase.makeGraph` for more information.
        """
        raise NotImplementedError()


def noop(value):
    """A simple pass-through function, used as a do-nothing converter
    for trivial types when mapping Python types to SQL types or vice-versa.
    """
    return value


def expandUnitClasses(UnitClasses):
    """Expand a Unit class sequence by recursively adding any
    classes used as ForeignKeys.
    """
    incomplete = set(UnitClasses)
    complete = set()
    everything = set(UnitClasses)
    while incomplete:
        UnitClass = incomplete.pop()
        for f in UnitClass.fields.values():
            if (isinstance(f, fields.ForeignKey) and
                    f.UnitClass not in everything):
                incomplete.add(f.UnitClass)
                everything.add(f.UnitClass)
        complete.add(UnitClass)
    return everything


class SqliteBackend(Backend):
    """A `Backend` implementation that utilizes a single SQLite database.

    Because SqliteBackend assumes a single database for what could be a set
    of multiple chained data repositories, it probably isn't usable long-term
    as-is.

    Parameters
    ----------
    dbname : `str`
        Name of a SQLite database file, or `:memory:` for an in-memory database
        (for e.g. testing purposes).
    config : `dict`
        Dictionary of configuration values passed to `Backend.__init__`.
    """

    # Mapping from Unit Field types to SQLite types.
    SQL_TYPES = {
        fields.RegionField: "BLOB",
        fields.IntField: "INTEGER",
        fields.StrField: "TEXT",
        fields.DateTimeField: "INTEGER",
        fields.ForeignKey: "INTEGER",
    }

    # Mapping from Unit Field type to a function that converts SQL types
    # to Python types.
    FROM_SQL = {
        fields.DateTimeField: datetime.datetime.fromtimestamp,
    }

    # Mapping from Unit Field type to a function that converts Python types
    # to SQL types.
    TO_SQL = {
        fields.DateTimeField:
            lambda dt: (dt - datetime.datetime(1970, 1, 1)).total_seconds(),
        fields.ForeignKey:
            lambda x: x.id
    }

    def __init__(self, dbname, config=None):
        Backend.__init__(self, config)
        self.dbname = dbname
        self.db = sqlite3.connect(dbname)
        self.db.row_factory = sqlite3.Row

    def __getstate__(self):
        return (self.dbname, self.config)

    def __setstate__(self, state):
        dbname, config = state
        self.__init__(dbname, config=config)

    def createUnitTable(self, UnitClass):
        """Create a table or view that can hold instances of the given `Unit`
        type.

        If the given unit type does not inherit from `Unit`, any intermediate
        base classes with one or more `Fields` must already have associated
        tables.
        """
        if len(UnitClass.__bases__) > 1:
            raise ValueError("Unit classes with multiple base classes are not supported.")
        if len(UnitClass.fields) == 0:
            # Do nothing for derived classes that don't add any fields.
            return
        idColDef = "{} INTEGER PRIMARY KEY".format(self.config["idcol"])
        BaseClass = UnitClass.__bases__[0]
        if BaseClass != Unit:
            idColDef += " REFERENCES {} ({})".format(
                self.getUnitTableName(BaseClass),
                self.config["idcol"]
            )
        items = [idColDef]
        for k, v in UnitClass.fields.items():
            t = [k, self.SQL_TYPES[type(v)]]
            if isinstance(v, fields.ForeignKey):
                t.append(
                    "REFERENCES {} ({})".format(
                        self.getUnitTableName(v.UnitClass),
                        self.config["idcol"]
                    )
                )
            if not v.optional:
                t.append("NOT NULL")
            items.append(" ".join(t))
        if UnitClass.unique is not BaseClass.unique and UnitClass.unique is not None:
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
        """Create a table or view that can hold instances of the given
        `Dataset` type.
        """
        items = []
        for k, v in DatasetClass.properties.items():
            t = [k, self.SQL_TYPES[fields.ForeignKey]]
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
        sql = "CREATE TABLE IF NOT EXISTS {} (\n    {}\n)".format(
            self.getDatasetTableName(DatasetClass), ",\n    ".join(items)
        )
        self.db.execute(sql)
        self.db.commit()

    def insertUnit(self, unit):
        """Insert a `Unit` instance into the appropriate table.
        """
        UnitClass = type(unit)
        sql, columns, converters = self._buildUnitInsertQuery(UnitClass, id=unit.id)
        values = tuple(convert(getattr(unit, k))
                       for k, convert in zip(columns, converters))
        cursor = self.db.cursor()
        cursor.execute(sql, values)
        if unit.id is None:
            unit._storage["id"] = cursor.lastrowid
        self.db.commit()

    def _buildUnitInsertQuery(self, UnitClass, id=None):
        """Build a SQL query string that inserts `Unit` instances into the
        appropriate table.

        In addition to the SQL query string, also returns a list of column
        names and a list of converter functions that should be applied
        to column values.
        """
        columns = [self.config['idcol']]
        converters = [noop]
        for f in UnitClass.fields.values():
            columns.append(f.name)
            converters.append(self.TO_SQL.get(type(f), noop))
        sql = "INSERT INTO {table} ({columns}) VALUES ({placeholders})".format(
            table=self.getUnitTableName(UnitClass),
            columns=", ".join(columns),
            placeholders=", ".join(["?"] * len(columns))
        )
        return sql, columns, converters

    def insertDataset(self, dataset):
        """Insert a `Dataset` instance into the appropriate table.
        """
        sql, columns = self._buildDatasetInsertQuery(type(dataset))
        values = tuple(getattr(dataset, k).id for k in columns)
        self.db.execute(sql, values)
        self.db.commit()

    def _buildDatasetInsertQuery(self, DatasetClass):
        """Build a SQL query string that inserts `Dataset` instances into the
        appropriate table.

        In addition to the SQL query string, also returns a list of column
        names.
        """
        columns = []
        for f in DatasetClass.properties.values():
            columns.append(f.name)
        sql = "INSERT INTO {table} ({columns}) VALUES ({placeholders})".format(
            table=self.getDatasetTableName(DatasetClass),
            columns=", ".join(columns),
            placeholders=", ".join(["?"] * len(columns))
        )
        return sql, columns

    def makeGraph(self, UnitClasses=(), where=None,
                  NeededDatasets=(), FutureDatasets=()):
        """Execute a multi-table SQL query and return the results
        as a `RepoGraph`.

        See `RepoDatabase.makeGraph` for more information.
        """
        UnitClasses = set(UnitClasses)
        for DatasetClass in itertools.chain(NeededDatasets, FutureDatasets):
            for p in DatasetClass.properties.values():
                UnitClasses.add(p.UnitClass)
        UnitClasses = expandUnitClasses(UnitClasses)
        temp = self._makeTemporaryGraphTable(UnitClasses, NeededDatasets,
                                             where=where)
        units = {}
        for UnitClass in UnitClasses:
            units[UnitClass] = self._readPartialUnits(UnitClass, temp)
        datasets = {}
        for DatasetClass in FutureDatasets:
            # This is probably a rather expensive way to ensure the table
            # exists; might want to consider an in-memory set of Dataset
            # classes that have tables in the DB.
            self.createDatasetTable(DatasetClass)
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
        """Create a temporary table containing the results of a SELECT
        query that returns integer primary key values for all `Unit` classes
        constrained by all joins and the user-provided WHERE clause.
        """
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
                if isinstance(f, fields.ForeignKey):
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
               "SELECT {columns} FROM {tables}").format(
            temp=temp,
            columns=", ".join(columns),
            tables=", ".join(tables),
        )
        if where:
            sql += " WHERE {}".format(" AND ".join(where))
        self.db.execute(sql)
        return temp

    def _readPartialUnits(self, UnitClass, temp):
        """Read 'partial' Units (with integer IDs for all ForeignKey fields)
        by executing a query joined with the temporary table created by
        _makeTemporaryGraphTable.
        """
        sql = ("SELECT DISTINCT {idcol}, {columns} FROM {table} "
               "INNER JOIN {temp} ON ({table}.{idcol} = {temp}.{table}_id)"
               ).format(
                    columns=", ".join(UnitClass.fields.keys()),
                    table=self.getUnitTableName(UnitClass),
                    idcol=self.config["idcol"],
                    temp=temp
                )
        converters = [self.FROM_SQL.get(type(f), noop)
                      for f in UnitClass.fields.values()]
        result = {}
        for row in self.db.execute(sql):
            id = row["id"]
            assert id not in result
            storage = {"id": id}
            for f, converter in zip(UnitClass.fields.values(), converters):
                storage[f.name] = converter(row[f.name])
            result[id] = UnitClass(**storage)
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
        """Read 'partial' Datasets (with integer IDs for all UnitPropertys)
        by executing a query joined with the temporary table created by
        _makeTemporaryGraphTable.
        """
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
            result.add(DatasetClass(**storage))
        return result

    def _finalizeDatasets(self, datasets, units):
        """ Iterate over newly created Dataset instances, turning integer
        Unit IDs into Unit instances and adding references from the Units
        back to these Datasets.
        """
        for DatasetClass, instances in datasets.items():
            for instance in instances:
                for p in DatasetClass.properties.values():
                    p.finalize(instance, units)
        return datasets
