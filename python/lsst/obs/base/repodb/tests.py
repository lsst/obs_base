from __future__ import print_function, division, absolute_import

from lsst.skymap import DiscreteSkyMap

from .backend import SqliteBackend
from .repodb import RepoDatabase, CameraDataSpec
from .datasets import Dataset
from . import common

__all__ = ("makeRepoDatabase", "Coadd")


HSC = CameraDataSpec("HSC", filters="grizy")

DISCRETE_2 = DiscreteSkyMap(
    config=DiscreteSkyMap.ConfigClass(
        raList=[40.0, 15.0],
        decList=[60.0, 26.0],
        radiusList=[1.0, 1.0]
    )
)

Coadd = Dataset.subclass(
    "Coadd",
    tract=common.TractUnit,
    patch=common.PatchUnit,
    filter=common.FilterUnit
)


def makeRepoDatabase(filename=":memory:"):
    backend = SqliteBackend(filename)
    db = RepoDatabase(backend)
    db.create()
    db.addCamera(HSC)
    db.addSkyMap(DISCRETE_2, "DISCRETE_2")
    db.addTracts("DISCRETE_2")
    db.registerDatasetType(Coadd)
    graph = db.makeGraph(db.UNIT_CLASSES)
    for filterUnit in graph.units[common.FilterUnit]:
        if filterUnit.name != "r":
            continue
        for tractUnit in graph.units[common.TractUnit]:
            for patchUnit in graph.units[common.PatchUnit]:
                db.addDataset(Coadd(filter=filterUnit, tract=tractUnit,
                                    patch=patchUnit))
    return db
