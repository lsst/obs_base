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
from collections import defaultdict

import lsst.daf.butler.tests as butlerTests
from lsst.daf.butler import DataCoordinate, DimensionRecord, SerializedDimensionRecord
from lsst.obs.base import DefineVisitsTask
from lsst.obs.base.instrument_tests import DummyCam
from lsst.utils.iteration import ensure_iterable

TESTDIR = os.path.dirname(__file__)
DATADIR = os.path.join(TESTDIR, "data", "visits")


class DefineVisitsTestCase(unittest.TestCase):
    """Test visit definition."""

    def setUp(self):
        """Create a new butler for each test since we are changing dimension
        records.
        """
        self.root = tempfile.mkdtemp(dir=TESTDIR)
        self.creatorButler = butlerTests.makeTestRepo(self.root, [])
        self.butler = butlerTests.makeTestCollection(self.creatorButler, uniqueId=self.id())

        self.config = DefineVisitsTask.ConfigClass()
        self.task = DefineVisitsTask(config=self.config, butler=self.butler)

        # Need to register the instrument.
        DummyCam().register(self.butler.registry)

        # Read the exposure records.
        self.records: dict[int, DimensionRecord] = {}
        for i in (347, 348, 349):
            with open(os.path.join(DATADIR, f"exp_{i}.json")) as fh:
                simple = SerializedDimensionRecord.model_validate_json(fh.read())
            self.records[i] = DimensionRecord.from_simple(simple, registry=self.butler.registry)

    def tearDown(self):
        if self.root is not None:
            shutil.rmtree(self.root, ignore_errors=True)

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

    def define_visits(
        self, exposures: list[DimensionRecord | list[DimensionRecord]], incremental: bool
    ) -> None:
        for records in exposures:
            self.butler.registry.insertDimensionData("exposure", *ensure_iterable(records))
            # Include all records so far in definition.
            dataIds = list(self.butler.registry.queryDataIds("exposure", instrument="DummyCam"))
            self.task.run(dataIds, incremental=incremental)

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
        self.assertEqual(self.task.getFullName(), copy.getFullName())
        self.assertEqual(self.task.log.name, copy.log.name)
        self.assertEqual(self.task.config, copy.config)
        self.assertEqual(self.task.butler._config, copy.butler._config)
        self.assertEqual(self.task.butler.collections, copy.butler.collections)
        self.assertEqual(self.task.butler.run, copy.butler.run)
        self.assertEqual(self.task.universe, copy.universe)


if __name__ == "__main__":
    unittest.main()
