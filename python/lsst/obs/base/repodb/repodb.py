from __future__ import print_function, division, absolute_import

from . import base
from . import common

__all__ = ("RepoDatabase",)


class RepoDatabase(object):
    """An interface to metadata in a repository, as represented by `Units`
    and `Datasets`.

    Design Notes
    ------------
    `RepoDatabase` is a concrete class that is aware of the concrete set of
    Units that are used when processing optical/NIR astronomical imaging data.
    It probably implicitly assumes that it's backed by a SQL database (at least
    it would probably be hard to implement without one), but it is unaware of
    the details of the schema or the DBMS; all of that is hidden behind the
    `Backend` interface.

    `RepoDatabase` plays a role somewhat similar to the current
    `obs.base.CameraMapper` class, which is the lowest layer in the current
    `Butler`/`Mapper` system that is aware of astronomical concepts.  Unlike
    '`CameraMapper`, however, it puts camera-based data units (visits,
    sensors) and SkyMap-based units (tracts, patches) on an even footing, and
    via `Units` it demands the per-camera specializations be mapped to a
    fairly rigid common data model that can be used by camera-generic
    algorithmic code.  A single `RepoDatabase` also explicitly supports
    multiple cameras (via `CameraDataSpec`) and multiple skymaps (via
    `lsst.skymap` objects).

    The ideal relationship between `RepoDatabase` and `Mapper` is as yet
    unclear. At least at first, `RepoDatabase` can work with the existing
    `Mappers` without requiring them to be changed.  In that mode, a
    `RepoDatabase` would sit "next to" a `Mapper` (or, more generally, the set
    of `Mapper`s in that repo's parent chain), and generate `Dataset`s that can
    be translated to `datasetType` strings and dictionary-style data IDs for
    use by the `Mapper`.  But as `RepoDatabase` becomes more capable, we
    probably want to retire the existing registry databases, and it would
    probably make sense to have `Mapper` use `RepoDatabase`'s in that role
    instead: dictionary-style data IDs would be expanded into fully-qualified
    `Datasets` by queries against a `RepoDatabase`, and these would be used
    directly to map to a location and storage for actual retrieval.  This
    *could* ultimately eliminate the need to specialize `Mapper` for each
    camera; all camera-specific content could move to `CameraDataSpec`.  I
    think it's more likely that we'll want to keep camera-specific `Mapper`s
    but trim them down to just YAML configurations that override the locations
    for a few datasets, with all camera-specific *code* moved to
    `CameraDataSpec`.

    At present, `RepoDatabase` just uses `lsst.skymap` objects directly, and
    considers that package's `TractInfo` and `PatchInfo` classes to be a
    parallel description of what's in it's own `TractUnit` and `PatchUnit`
    classes.  That redundancy is more confusing than helpful, though, and once
    the rest of the `RepoDatabase` design settles down, it'd make sense to
    unify the tract and patch classes.  A particular `SkyMap` class would then
    just be responsible for generating a set of `TractUnit` and `PatchUnit`
    classes that support the full functionality of `TractInfo` and
    `PatchInfo`, which would then be stored directly in the `RepoDatabase`
    itself (including their WCSs and detailed bounding boxes) instead of as a
    separate `Butler`-accessed dataset.

    Because `RepoDatabase` contains Python objects (types, in particular)
    as well as a SQL database connection, we need a way to perist and unpersist
    it to a special file stored in a repository.  At present we use pickle,
    but it probably should be stored as part of the YAML repository
    configuration used by `Butler`.

    The other major remaining design challenges for `RepoDatabase` are:

     - How do we split storage across multiple chained repositories?  The
       current design represents the entire content of a repository (including
       its parents) via a single database backend.  This *should* be just a
       `Backend`/`Butler` problem: a `RepoDatabase` shouldn't care how it is
       stored, and it should (in the future) be the responsibility of a
       `Butler` to construct a `Backend` for a particular endpoint repository
       and a `RepoDatabase` from that.  But it's quite possible there's
       something in the current interface between `RepoDatabase` and `Backend`
       that will need to change to keep that separation of concerns.

     - How do we let cameras specialize `Unit`s, and in particular add their
       own labels for camera-generic concepts like `Visit`?  The general plan
       is to add per-camera tables for `Unit`s that are identified as
       belonging to a camera; these would be joined with the camera-generic
       tables for those `Unit`s in the queries run by `makeGraph`, allowing
       camera-specific labels to be used in the where clause that represents
       user data ID expression.  The details of that need to be worked out.

     - How do we supports date-based joins between calibration `Unit`s and
       visit/sensor `Unit`s?

     - How do we support spatial joins between visit/sensors `Unit`s, skymap
       `Unit`s, reference-catalog shard `Unit`s, and any other partition of
       the sky?

     - Do we need, and if so, how do we support the addition of new `Unit`
       instances to a `RepoDatabase` by a `SuperTask` supervisory framework?
       Normally, new `Unit`s are added by "ingest-like" steps (this includes
       both ingesting raw data, ingesting raw calibrations, and adding a new
       SkyMap), which define new `Unit`s and only add `Dataset`s that
       represent data products that already exist.  That model may not work
       for master calibration data products, which are identified by date-like
       `Unit`s that it doesn't make sense to ask the user to "ingest" in
       advance of actually running the `Pipeline` that produces them.

     - Do we need, and if so, how do we support the addition of new `Unit`
       *types* without modifying the `RepoDatabase` implementation itself?
       In discussions of an early version of this design, some concern about
       losing the flexibility to define new data ID keys was expressed, and
       while a clear use case for user-defined `Unit` classes has not been
       identified, it might be prudent to find a way to support them, even
       if the new `Unit`s are not treated exactly the same way as those
       a `RepoDatabase` is intrinsically aware of.
    """

    DEFAULT_UNIT_CLASSES = (common.CameraUnit, common.SkyMapUnit,
                            common.TractUnit, common.PatchUnit,
                            common.AbstractFilterUnit, common.PhysicalFilterUnit,
                            common.VisitUnit, common.SensorUnit)

    def __init__(self, backend):
        self.backend = backend
        self._cameras = {}
        self._skyMaps = {}
        self.UnitClasses = set(self.DEFAULT_UNIT_CLASSES)
        self.DatasetClasses = set()

    def registerUnitClass(self, UnitClass):
        RootUnit, hasTable = base.categorizeUnit(UnitClass)
        if RootUnit and RootUnit not in self.UnitClasses:
            self.backend.createUnitTable(RootUnit)
            self.UnitClasses.add(RootUnit)
        if hasTable and UnitClass not in self.UnitClasses:
            self.backend.createUnitTable(UnitClass)
            self.UnitClasses.add(UnitClass)

    def registerDatasetClass(self, DatasetClass):
        if DatasetClass not in self.DatasetClasses:
            self.backend.createUnitTable(DatasetClass)
            self.DatasetClasses.add(DatasetClass)

    def create(self):
        """Create all `Unit` tables required by the `RepoDatabase`.

        This should only be once when a `RepoDatabase` is first constructed
        (not merely unpersisted).
        """
        for UnitClass in self.DEFAULT_UNIT_CLASSES:
            self.backend.createUnitTable(UnitClass)

    def addCamera(self, camera):
        """Add a `CameraUnit` and its associated `FilterUnits` and
        `SensorUnits`.
        """
        self.registerUnitClass(type(camera))
        self.insertUnit(camera)
        self._cameras[camera.name] = camera
        camera.register(self)

    def insertUnit(self, unit):
        self.backend.insertUnit(unit)

    def addSkyMap(self, skyMap, name):
        """Add `Unit`s to the `RepoDatabase defined by a `SkyMap`.

        This adds a `SkyMapUnit` to the databse, enabling the user
        to call `addTracts` to actually add `TractUnit` and `PatchUnit`
        instances to the database.

        Parameters
        ----------
        skyMap : subclass of `lsst.skymap.BaseSkyMap`
            An object that describes a set of tracts and patches that tile
            the sky.
        name : `str`
            A unique name for this skymap.  This needs to uniquely identify
            the skymap *instance* (i.e. including configuration), not just
            its type.
        """
        skyMapUnit = common.SkyMapUnit(name=name)
        self.insertUnit(skyMapUnit)
        self._skyMaps[name] = (skyMap, skyMapUnit)

    def addTracts(self, skyMapName, only=None):
        """Add `TractUnit` and `PatchUnit` instances to the database.

        Parameters
        ----------
        skyMapName : `str`
            Name the skymap that generates these tracts was registered with
            in the call to `addSkyMap`.
        only : sequence of `int`
            A list of `lsst.skymap` tract IDs (i.e. `TractUnit.number`) values
            to limit which tracts to add.  `None` (default) adds all tracts.
        """
        skyMap, skyMapUnit = self._skyMaps[skyMapName]
        allPatches = set()
        if only is None:
            iterable = skyMap
        else:
            iterable = (skyMap[t] for t in only)
        for tract in iterable:
            tractUnit = common.TractUnit(number=tract.getId(),
                                         skymap=skyMapUnit)
            self.insertUnit(tractUnit)
            for patch in tract:
                x, y = patch.getIndex()
                allPatches.add((x, y))
        for x, y in allPatches:
            patchUnit = common.PatchUnit(x=x, y=y, skymap=skyMapUnit)
            self.insertUnit(patchUnit)
        # TODO: tract-patch join table

    def registerDatasetType(self, DatasetClass):
        """Add a table for a new `Dataset` type to the database.

        This is a no-op if the dataset already exists.
        """
        self.backend.createDatasetTable(DatasetClass)

    def addDataset(self, dataset):
        """Add an instance of a `Dataset` to the database.

        This should only be called when the corresponding data product actually
        exists in the repository, i.e. during ingest or after a successful call
        to `Butler.put`.
        """
        self.backend.insertDataset(dataset)

    def makeGraph(self, UnitClasses=(), where=None,
                  NeededDatasets=(), FutureDatasets=()):
        """Create a `RepoGraph` that represents a possibly-restricted view
        into the database.

        Parameters
        ----------
        UnitClasses : sequence of type objects that inherit from `Unit`
            Include at least these unit types in the graph, which naturally
            restricts the graph to the intersection (across the predefined
            relationships between these units) of that is in the database for
            all of these units.  This sequence is expanded to include any unit
            type related to the unit types in the sequence and any units
            related to any `Dataset` types in `NeededDatasets` or
            `FutureDatasets`, and hence can frequently be an empty sequence.
        where : `str`
            An optional SQL where clause operating on the tables for the
            `Unit`s and `NeededDataset`s that restricts the graph.
        NeededDatasets : sequence of type objects that inherit from `Dataset`.
            Include these `Dataset` types in the graph, and restrict the
            graph to the intersection of the instances of these `Datasets`
            that already exist in the database.  Typically this should be the
            set of pure input datasets needed by a `Pipeline`.
        FutureDatasets : sequence of type objects that inherit from `Dataset`.
            Include these `Dataset` types in the graph, but do not restrict
            the graph based on whether they already exist in the database.
            Typically this should be the set of datasets produced by a
            `Pipeline`.

        Design Notes
        ------------
        The `where` argument currently requires the user to know about the
        actual database schema.  We need to abstract this user-provided
        expression somehow and have the `Backend` turn it into SQL.
        """
        return self.backend.makeGraph(
            UnitClasses=UnitClasses, where=where,
            NeededDatasets=NeededDatasets,
            FutureDatasets=FutureDatasets
        )
