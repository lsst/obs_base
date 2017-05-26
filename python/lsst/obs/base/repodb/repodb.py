from __future__ import print_function, division, absolute_import

from . import common
from . import graph


class Camera(object):

    def __init__(self, name, filters):
        self.name = name
        self.filters = filters


class RepoDatabase(object):

    UNIT_CLASSES = (common.CameraUnit, common.SkyMapUnit,
                    common.TractUnit, common.PatchUnit,
                    common.FilterUnit,)

    def __init__(self, backend):
        self.backend = backend
        self._cameras = {}
        self._skymaps = {}

    def create(self):
        for UnitClass in self.UNIT_CLASSES:
            self.backend.createUnitTable(UnitClass)

    def addCamera(self, camera):
        cameraUnit = common.CameraUnit(name=camera.name)
        self.backend.insertUnit(cameraUnit)
        self._cameras[camera.name] = (camera, cameraUnit)
        for f in camera.filters:
            filterUnit = common.FilterUnit(name=f, camera=cameraUnit)
            self.backend.insertUnit(filterUnit)
        # TODO: add table for raw Dataset type

    def addSkyMap(self, skyMap, name):
        # n.b. name has to uniquely identify SkyMap + configuration, so it
        # isn't just mapped to the name of the class.
        skyMapUnit = common.SkyMapUnit(name=name)
        self.backend.insertUnit(skyMapUnit)
        self._skymaps[name] = (skyMap, skyMapUnit)
        allPatches = set()
        for tract in skyMap:
            tractUnit = common.TractUnit(number=tract.getId(),
                                         skymap=skyMapUnit)
            self.backend.insertUnit(tractUnit)
            for patch in tract:
                x, y = patch.getIndex()
                allPatches.add((x, y))
        for x, y in allPatches:
            patchUnit = common.PatchUnit(x=x, y=y, skymap=skyMapUnit)
            self.backend.insertUnit(patchUnit)
        # TODO: tract-patch join table

    def registerDatasetType(self, DatasetClass):
        self.backend.createDatasetTable(DatasetClass)

    def addDataset(self, dataset):
        self.backend.insertDataset(dataset)

    def makeGraph(self, UnitClasses, where=None,
                  NeededDatasets=(), FutureDatasets=()):
        return self.backend.makeGraph(
            UnitClasses, where=where,
            NeededDatasets=NeededDatasets,
            FutureDatasets=FutureDatasets
        )
