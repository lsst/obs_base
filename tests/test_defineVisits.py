# This file is part of obs_base.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import pickle
import shutil
import tempfile
import unittest
import warnings
from collections import defaultdict

import lsst.daf.butler.tests as butlerTests
from lsst.daf.butler import DataCoordinate, DimensionRecord, SerializedDimensionRecord
from lsst.daf.butler.registry import ConflictingDefinitionError
from lsst.obs.base import DefineVisitsConfig, DefineVisitsTask
from lsst.obs.base.instrument_tests import DummyCam
from lsst.utils.iteration import ensure_iterable

TESTDIR = os.path.dirname(__file__)
DATADIR = os.path.join(TESTDIR, "data", "visits")


class DefineVisitsBase:
    """General set up that can be shared."""

    use_data_ids = True
    """Use data IDs when calling defineVisits, else use dimension records."""

    def setUpExposures(self):
        """Create a new butler for each test since we are changing dimension
        records.
        """
        self.root = tempfile.mkdtemp(dir=TESTDIR)
        self.creatorButler = butlerTests.makeTestRepo(self.root, {})
        self.enterContext(self.creatorButler)
        self.butler = butlerTests.makeTestCollection(self.creatorButler, uniqueId=self.id())
        self.enterContext(self.butler)

        self.config = self.get_config()
        self.task = DefineVisitsTask(config=self.config, butler=self.butler)

        # Need to register the instrument.
        DummyCam().register(self.butler.registry)

        # Choose serializations based on universe.
        universe = self.butler.dimensions
        uversion = universe.version
        # Not all universe changes result in visible changes.
        match uversion:
            case uversion if uversion < 2:
                raise unittest.SkipTest(f"Universe {uversion} is not compatible with these test files.")
            case 2 | 3 | 4 | 5:
                # has_simulated, azimuth, seq_start, seq_end.
                v = 2
            case 6:
                # group not group_name, group_id dropped.
                v = 6
            case 7:
                # can_see_sky.
                v = 7
            case _:
                # Might work.
                warnings.warn(f"Universe {uversion} has not been validated.")
                v = 7

        # Read the exposure records.
        self.records: dict[int, DimensionRecord] = {}
        for i in (347, 348, 349):
            with open(os.path.join(DATADIR, f"exp_v{v}_{i}.json")) as fh:
                simple = SerializedDimensionRecord.model_validate_json(fh.read())
            self.records[i] = DimensionRecord.from_simple(simple, registry=self.butler.registry)

    def define_visits(
        self,
        exposures: list[DimensionRecord | list[DimensionRecord]],
        incremental: bool,
    ) -> None:
        for records in exposures:
            records = list(ensure_iterable(records))
            if "group" in self.butler.dimensions["exposure"].implied:
                # This is a group + day_obs universe.
                for rec in records:
                    self.butler.registry.syncDimensionData(
                        "group", dict(instrument=rec.instrument, name=rec.group)
                    )
                    self.butler.registry.syncDimensionData(
                        "day_obs", dict(instrument=rec.instrument, id=rec.day_obs)
                    )

            self.butler.registry.insertDimensionData("exposure", *records)
            # Include all records so far in definition.
            if self.use_data_ids:
                dataIds = list(self.butler.registry.queryDataIds("exposure", instrument="DummyCam"))
            else:
                dataIds = records
            self.task.run(dataIds, incremental=incremental)


class DefineVisitsTestCase(unittest.TestCase, DefineVisitsBase):
    """Test visit definition."""

    def setUp(self):
        self.setUpExposures()

    def tearDown(self):
        if self.root is not None:
            shutil.rmtree(self.root, ignore_errors=True)

    def get_config(self) -> DefineVisitsConfig:
        config = DefineVisitsTask.ConfigClass()
        config.groupExposures.name = "one-to-one-and-by-counter"
        return config

    def assertVisits(self):
        """Check that the visits were registered as expected."""
        visits = list(self.butler.registry.queryDimensionRecords("visit"))
        self.assertEqual(len(visits), 4)
        self.assertEqual(
            {visit.id for visit in visits}, {2022040500347, 2022040500348, 2022040500349, 92022040500348}
        )

        # Ensure that the definitions are correct (ignoring order).
        defmap = defaultdict(set)
        definitions = list(self.butler.registry.queryDimensionRecords("visit_definition"))
        for defn in definitions:
            defmap[defn.visit].add(defn.exposure)

        self.assertEqual(
            dict(defmap),
            {
                92022040500348: {2022040500348},
                2022040500347: {2022040500347},
                2022040500348: {2022040500348, 2022040500349},
                2022040500349: {2022040500349},
            },
        )

    def test_defineVisits(self):
        # Test visit definition with all the records.
        self.define_visits([list(self.records.values())], incremental=False)  # list inside a list
        self.assertVisits()

    def test_incremental_cumulative(self):
        # Define the visits after each exposure.
        self.define_visits(list(self.records.values()), incremental=True)
        self.assertVisits()

    def test_incremental_cumulative_reverse(self):
        # In reverse order we should still eventually end up with the right
        # answer.
        with self.assertLogs("lsst.defineVisits.groupExposures", level="WARNING") as cm:
            self.define_visits(list(reversed(self.records.values())), incremental=True)
        self.assertIn("Skipping the multi-snap definition", "\n".join(cm.output))
        self.assertVisits()

    def define_visits_incrementally(self, exposure: DimensionRecord) -> None:
        if "group" in self.butler.dimensions["exposure"].implied:
            self.butler.registry.syncDimensionData(
                "group", dict(instrument=exposure.instrument, name=exposure.group)
            )
            self.butler.registry.syncDimensionData(
                "day_obs",
                dict(
                    instrument=exposure.instrument,
                    id=exposure.day_obs,
                ),
            )
        self.butler.registry.insertDimensionData("exposure", exposure)
        dataIds = [
            DataCoordinate.standardize(
                instrument="DummyCam", exposure=exposure.id, universe=self.butler.dimensions
            )
        ]
        self.task.run(dataIds, incremental=True)

    def test_incremental(self):
        for record in self.records.values():
            self.define_visits_incrementally(record)
        self.assertVisits()

    def test_incremental_reverse(self):
        for record in reversed(self.records.values()):
            self.define_visits_incrementally(record)
        self.assertVisits()

    def testPickleTask(self):
        stream = pickle.dumps(self.task)
        copy = pickle.loads(stream)
        self.enterContext(copy.butler)
        self.assertEqual(self.task.getFullName(), copy.getFullName())
        self.assertEqual(self.task.log.name, copy.log.name)
        self.assertEqual(self.task.config, copy.config)
        self.assertEqual(self.task.butler._config, copy.butler._config)
        self.assertEqual(list(self.task.butler.collections.defaults), list(copy.butler.collections.defaults))
        self.assertEqual(self.task.butler.run, copy.butler.run)
        self.assertEqual(self.task.universe, copy.universe)


class DefineVisitsRecordsTestCase(DefineVisitsTestCase):
    """Define visits using only dimension records."""

    use_data_ids = False


class DefineVisitsGroupingTestCase(unittest.TestCase, DefineVisitsBase):
    """Test visit grouping by group metadata."""

    def setUp(self):
        self.setUpExposures()

    def tearDown(self):
        if self.root is not None:
            shutil.rmtree(self.root, ignore_errors=True)

    def get_config(self) -> DefineVisitsConfig:
        config = DefineVisitsTask.ConfigClass()
        config.groupExposures.name = "by-group-metadata"
        return config

    def test_defineVisits(self):
        # Test visit definition with all the records.
        self.define_visits([list(self.records.values())], incremental=False)  # list inside a list
        self.assertVisits()

    def assertVisits(self):
        """Check that the visits were registered as expected."""
        visits = list(self.butler.registry.queryDimensionRecords("visit"))
        self.assertEqual(len(visits), 2)

        # The visit ID itself depends on which universe we are using.
        # It is either calculated or comes from the JSON record.
        if "group" in self.butler.dimensions["exposure"].implied:
            visit_ids = [20220406025653255, 20220406025807181]
        else:
            visit_ids = [2291434132550000, 2291434871810000]
        self.assertEqual({visit.id for visit in visits}, set(visit_ids))

        # Ensure that the definitions are correct (ignoring order).
        defmap = defaultdict(set)
        definitions = list(self.butler.registry.queryDimensionRecords("visit_definition"))
        for defn in definitions:
            defmap[defn.visit].add(defn.exposure)

        self.assertEqual(
            dict(defmap),
            {
                visit_ids[0]: {2022040500347},
                visit_ids[1]: {2022040500348, 2022040500349},
            },
        )


class DefineVisitsGroupingRecordsTestCase(DefineVisitsGroupingTestCase):
    """Test using dimension records instead of Data IDs."""

    use_data_ids = False


class DefineVisitsOneToOneTestCase(unittest.TestCase, DefineVisitsBase):
    """Test visit grouping by group metadata."""

    def setUp(self):
        self.setUpExposures()

    def tearDown(self):
        if self.root is not None:
            shutil.rmtree(self.root, ignore_errors=True)

    def get_config(self) -> DefineVisitsConfig:
        config = DefineVisitsTask.ConfigClass()
        config.groupExposures.name = "one-to-one"
        return config

    def test_defineVisits(self):
        # Test visit definition with all the records.
        self.define_visits([list(self.records.values())], incremental=False)  # list inside a list
        self.assertVisits()

    def assertVisits(self):
        """Check that the visits were registered as expected."""
        visits = list(self.butler.registry.queryDimensionRecords("visit"))
        self.assertEqual(len(visits), 3)

        # For one-to-one the visit ID is the exposure ID.
        visit_ids = [rec.id for rec in self.records.values()]
        self.assertEqual({visit.id for visit in visits}, set(visit_ids))

        # Ensure that the definitions map an exposure ID to an identical visit
        # ID.
        definitions = list(self.butler.registry.queryDimensionRecords("visit_definition"))
        for defn in definitions:
            self.assertEqual(defn.visit, defn.exposure)

    def test_update_records(self):
        self.define_visits([list(self.records.values())], incremental=False)  # list inside a list
        self.assertVisits()

        # Modify one of the records.
        records = self.records
        simple = records[348].to_simple()
        simple.record["target_name"] = "new target"
        records[348] = DimensionRecord.from_simple(simple, universe=self.butler.dimensions)
        self.butler.registry.syncDimensionData("exposure", records[348], update=True)

        # Re-run without updates or skipping should fail.
        with self.assertRaises(ConflictingDefinitionError):
            self.task.run(records.values())

        result = self.task.run(records.values(), skip_conflicting=True)
        self.assertEqual(result.n_skipped, 3, str(result))

        # Check that the visit definition did not change.
        visit_348 = self.butler.query_dimension_records("visit", where="visit.seq_num = 348")[0]
        self.assertEqual(visit_348.target_name, "LATISS_E6A_00000040", visit_348)

        # Run with forced update.
        result = self.task.run(records.values(), skip_conflicting=True, update_records=True)

        # Every record reports it was updated if we are updating, even if
        # a record was not really changed.
        self.assertEqual(result.n_skipped, 0, str(result))
        self.assertEqual(result.n_fully_updated, 3, str(result))
        visit_348 = self.butler.query_dimension_records("visit", where="visit.seq_num = 348")[0]
        self.assertEqual(visit_348.target_name, "new target", visit_348)


if __name__ == "__main__":
    unittest.main()
