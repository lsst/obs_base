from __future__ import print_function, division, absolute_import

import itertools
import datetime

from .backend import SqliteBackend
from . import common


__all__ = ("EXAMPLE_DATA", "makeBackend")


def makeDateTime(hours):
    dt = datetime.datetime(2016, 5, 10, hours)
    return (dt - datetime.datetime(1970, 1, 1)).total_seconds()

EXAMPLE_DATA = [(
        common.TractUnit,
        ("id", "number", "skymap"),
        [(1, 8766, "RINGS_120"),
         (2, 8767, "RINGS_120")]
    ), (
        common.PatchUnit,
        ("id", "tract", "x", "y"),
        [(n, tract, x, y) for n, (tract, x, y) in enumerate(
            itertools.product((1, 2), range(8), range(8))
        )]
    ), (
        common.FilterUnit,
        ("id", "name", "camera"),
        [(n, f, "HSC") for n, f in enumerate("gri")]
    ), (
        common.VisitUnit,
        ("id", "number", "filter", "camera", "dateobs"),
        [(1, 1001, 0, "HSC", makeDateTime(10)),
         (2, 1002, 0, "HSC", makeDateTime(11)),
         (3, 1003, 1, "HSC", makeDateTime(12)),
         (4, 1004, 2, "HSC", makeDateTime(13))],
    ), (
        common.SensorUnit,
        ("id", "number", "visit"),
        [(n, number, visit) for n, (number, visit) in enumerate(
            itertools.product(range(45, 55), (1, 2, 3, 4))
         )]
    )]


def makeBackend(filename=":memory:"):
    backend = SqliteBackend(filename)
    for UnitClass in common.COMMON_UNITS:
        backend.createTable(UnitClass)
    for UnitClass, names, values in EXAMPLE_DATA:
        tableName = backend.getTableName(UnitClass)
        sql = "INSERT INTO {} ({}) VALUES ({})".format(
            tableName,
            ", ".join(names),
            ", ".join(["?"]*len(names))
        )
        backend.db.executemany(sql, values)
    backend.db.commit()
    return backend
