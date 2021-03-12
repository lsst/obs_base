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

import lsst.log
import lsst.daf.butler.tests as butlerTests
from lsst.daf.butler import DatasetType, Butler, DataCoordinate, Config
from lsst.daf.butler.registry import ConflictingDefinitionError
from lsst.daf.butler.core.utils import getFullTypeName

from lsst.obs.base.ingest_tests import IngestTestBase
from lsst.obs.base.instrument_tests import DummyCam
from lsst.obs.base import RawIngestTask


TESTDIR = os.path.abspath(os.path.dirname(__file__))
INGESTDIR = os.path.join(TESTDIR, "data", "ingest")


class DummyCamRawIngestTask(RawIngestTask):
    """For DummyCam we ingest a different dataset type that can return
    a non-Exposure."""

    def getDatasetType(self):
        """Return the DatasetType of the datasets ingested by this Task.
        """
        return DatasetType("raw_dict", ("instrument", "detector", "exposure"), "StructuredDataDict",
                           universe=self.butler.registry.dimensions)


class RawIngestTestCase(IngestTestBase, unittest.TestCase):
    """Test ingest using JSON sidecar files."""

    ingestDatasetTypeName = "raw_dict"
    rawIngestTask = getFullTypeName(DummyCamRawIngestTask)
    curatedCalibrationDatasetTypes = ()
    ingestDir = TESTDIR
    instrumentClassName = "lsst.obs.base.instrument_tests.DummyCam"
    file = os.path.join(INGESTDIR, "sidecar_data", "dataset_1.yaml")
    dataIds = [dict(instrument="DummyCam", exposure=100, detector=0)]

    @property
    def visits(self):
        butler = Butler(self.root, collections=[self.outputRun])
        return {
            DataCoordinate.standardize(
                instrument="DummyCam",
                visit=100,
                universe=butler.registry.dimensions
            ): [
                DataCoordinate.standardize(
                    instrument="DummyCam",
                    exposure=100,
                    universe=butler.registry.dimensions
                )
            ]
        }

    def testWriteCuratedCalibrations(self):
        """There are no curated calibrations in this test instrument"""
        pass


class RawIngestImpliedIndexTestCase(RawIngestTestCase):
    """Test ingest using JSON index files."""
    file = os.path.join(INGESTDIR, "indexed_data", "dataset_1.yaml")


class RawIngestEdgeCaseTestCase(unittest.TestCase):
    """Test ingest using non-standard approaches including failures."""

    @classmethod
    def setUpClass(cls):
        butlerConfig = """
datastore:
  # Want to ingest real files so can't use in-memory datastore
  cls: lsst.daf.butler.datastores.fileDatastore.FileDatastore
"""
        cls.root = tempfile.mkdtemp(dir=TESTDIR)
        cls.creatorButler = butlerTests.makeTestRepo(cls.root, {}, config=Config.fromYaml(butlerConfig))
        DummyCam().register(cls.creatorButler.registry)

    @classmethod
    def tearDownClass(cls):
        if cls.root is not None:
            shutil.rmtree(cls.root, ignore_errors=True)

    def setUp(self):
        self.butler = butlerTests.makeTestCollection(self.creatorButler)
        self.outputRun = self.butler.run

        config = RawIngestTask.ConfigClass()
        self.task = DummyCamRawIngestTask(config=config, butler=self.butler)

        # Different test files.
        self.bad_metadata_file = os.path.join(TESTDIR, "data", "small.fits")
        self.good_file = os.path.join(INGESTDIR, "sidecar_data", "dataset_2.yaml")
        self.bad_instrument_file = os.path.join(TESTDIR, "data", "calexp.fits")

    def testSimpleIngest(self):
        # Use the default per-instrument run for this.
        self.task.run([self.good_file])
        datasets = list(self.butler.registry.queryDatasets("raw_dict", collections="DummyCam/raw/all"))
        self.assertEqual(len(datasets), 1)

        # Now parallelized.
        files = [self.good_file,
                 os.path.join(INGESTDIR, "sidecar_data", "dataset_1.yaml")]
        self.task.run(files, processes=2, run=self.outputRun)
        datasets = list(self.butler.registry.queryDatasets("raw_dict", collections=self.outputRun))
        self.assertEqual(len(datasets), 2)

    def testExplicitIndex(self):
        files = [os.path.join(INGESTDIR, "indexed_data", "_index.json")]
        self.task.run(files, run=self.outputRun)

        datasets = list(self.butler.registry.queryDatasets("raw_dict", collections=self.outputRun))
        self.assertEqual(len(datasets), 2)

        # Try again with an explicit index and a file that is in that index.
        files.append(os.path.join(INGESTDIR, "indexed_data", "dataset_2.yaml"))
        new_run = self.outputRun + "b"
        self.task.run(files, run=new_run)

        datasets = list(self.butler.registry.queryDatasets("raw_dict", collections=self.outputRun))
        self.assertEqual(len(datasets), 2)

        # Now with two index files that point to the same files.
        # Look for the warning from duplication.
        files = [os.path.join(INGESTDIR, "indexed_data", "_index.json"),
                 os.path.join(INGESTDIR, "indexed_data", "translated_subdir", "_index.json")]
        new_run = self.outputRun + "c"

        with self.assertLogs(level="WARNING") as cm:
            with lsst.log.UsePythonLogging():
                self.task.run(files, run=new_run)
        self.assertIn("already specified in an index file, ignoring content", cm.output[0])

        datasets = list(self.butler.registry.queryDatasets("raw_dict", collections=self.outputRun))
        self.assertEqual(len(datasets), 2)

        # Again with an index file of metadata and one of translated.
        # Translated should win.
        # Put the metadata one first to test that order is preserved.
        files = [os.path.join(INGESTDIR, "indexed_data", "metadata_subdir", "_index.json"),
                 os.path.join(INGESTDIR, "indexed_data", "_index.json")]
        new_run = self.outputRun + "d"
        with self.assertLogs(level="WARNING") as cm:
            with lsst.log.UsePythonLogging():
                self.task.run(files, run=new_run)
        self.assertIn("already specified in an index file but overriding", cm.output[0])

        # Reversing the order should change the warning.
        # Again with an index file of metadata and one of translated.
        # Translated should win.
        # Put the metadata one first to test that order is preserved.
        files = [os.path.join(INGESTDIR, "indexed_data", "_index.json"),
                 os.path.join(INGESTDIR, "indexed_data", "metadata_subdir", "_index.json")]

        new_run = self.outputRun + "e"
        with self.assertLogs(level="WARNING") as cm:
            with lsst.log.UsePythonLogging():
                self.task.run(files, run=new_run)
        self.assertIn("already specified in an index file, ignoring", cm.output[0])

        # Bad index file.
        files = [os.path.join(INGESTDIR, "indexed_data", "bad_index", "_index.json")]
        with self.assertRaises(RuntimeError):
            self.task.run(files, run=self.outputRun)

        # Bad index file due to bad instrument.
        files = [os.path.join(INGESTDIR, "indexed_data", "bad_instrument", "_index.json")]
        with self.assertLogs(level="WARNING") as cm:
            with lsst.log.UsePythonLogging():
                with self.assertRaises(RuntimeError):
                    self.task.run(files, run=self.outputRun)
        self.assertIn("Instrument HSC for file", cm.output[0])

    def testBadExposure(self):
        """Test that bad exposures trigger the correct failure modes.

        This is the only test that uses the bad definition of dataset 4
        because exposure definitions are defined globally in a butler registry.
        """

        # Ingest 3 files. 2 of them will implicitly find an index and one
        # will use a sidecar.
        files = [os.path.join(INGESTDIR, "indexed_data", f"dataset_{n}.yaml") for n in (1, 2, 3)]
        new_run = self.outputRun
        self.task.run(files, run=new_run)

        datasets = list(self.butler.registry.queryDatasets("raw_dict", collections=new_run))
        self.assertEqual(len(datasets), 3)

        # Test fail fast.
        self.task.config.failFast = True

        # Ingest files with conflicting exposure definitions.
        # Ingest 3 files. One of them will implicitly find an index and one
        # will use a sidecar. The 3rd will fail due to exposure conflict.
        files = [os.path.join(INGESTDIR, "indexed_data", f"dataset_{n}.yaml") for n in (1, 3, 4)]
        new_run = self.outputRun + "_bad_exposure"
        with self.assertRaises(ConflictingDefinitionError):
            self.task.run(files, run=new_run)

    def testBadFile(self):
        """Try to ingest a bad file."""
        files = [self.bad_metadata_file]

        with self.assertRaises(RuntimeError) as cm:
            # Default is to raise an error at the end.
            self.task.run(files, run=self.outputRun)
        self.assertIn("Some failures", str(cm.exception))

        # Including a good file will result in ingest working but still
        # raises (we might want to move this to solely happen in the
        # command line invocation).
        files.append(self.good_file)

        # Also include a file with unknown instrument.
        files.append(self.bad_instrument_file)

        with self.assertRaises(RuntimeError):
            self.task.run(files, run=self.outputRun)
        datasets = list(self.butler.registry.queryDatasets("raw_dict", collections=self.outputRun))
        self.assertEqual(len(datasets), 1)

        # Fail fast will trigger a run time error with different text.
        # Use a different output run to be sure we are not failing because
        # of the attempt to ingest twice.
        self.task.config.failFast = True
        new_run = self.outputRun + "b"
        with self.assertRaises(RuntimeError) as cm:
            self.task.run([self.bad_metadata_file, self.good_file], run=new_run)
        self.assertIn("Problem extracting metadata", str(cm.exception))

        # Attempt to ingest good file again -- this will fail for a different
        # reason than failed metadata extraction.
        with self.assertRaises(ConflictingDefinitionError):
            self.task.run([self.good_file], run=self.outputRun)

        # Ingest a file with good metadata but unknown instrument.
        with self.assertRaises(RuntimeError) as cm:
            self.task.run([self.bad_instrument_file], run=self.outputRun)
        self.assertIn("Instrument HSC", str(cm.exception))

        # Ingest of a metadata index file that will fail translation.
        with self.assertRaises(RuntimeError) as cm:
            self.task.run([os.path.join(INGESTDIR, "indexed_data", "metadata_subdir", "_index.json")])
        self.assertIn("Problem extracting metadata", str(cm.exception))

        # Ingest of a bad index file.
        with self.assertRaises(RuntimeError) as cm:
            self.task.run([os.path.join(INGESTDIR, "indexed_data", "bad_index", "_index.json")])
        self.assertIn("Problem reading index file", str(cm.exception))

        # Ingest of an implied bad index file.
        with self.assertRaises(RuntimeError) as cm:
            self.task.run([os.path.join(INGESTDIR, "indexed_data", "bad_implied", "dataset_2.yaml")])

    def testCallbacks(self):
        """Test the callbacks for failures."""

        # Define the callbacks.
        metadata_failures = []
        successes = []
        ingest_failures = []

        def on_metadata_failure(filename, exc):
            metadata_failures.append(filename)

        def on_success(datasets):
            successes.append(datasets)

        def on_ingest_failure(exposure, exc):
            ingest_failures.append(exposure)

        # Need our own task instance
        config = RawIngestTask.ConfigClass()
        self.task = DummyCamRawIngestTask(config=config, butler=self.butler,
                                          on_metadata_failure=on_metadata_failure,
                                          on_success=on_success,
                                          on_ingest_failure=on_ingest_failure)

        files = [self.good_file, self.bad_metadata_file, self.bad_instrument_file]

        with self.assertRaises(RuntimeError):
            self.task.run(files, run=self.outputRun)

        self.assertEqual(len(successes), 1)
        self.assertEqual(len(metadata_failures), 2)
        self.assertEqual(len(ingest_failures), 0)

        # Try the good one a second time.
        with self.assertRaises(RuntimeError):
            self.task.run([self.good_file], run=self.outputRun)

        self.assertEqual(len(successes), 1)
        self.assertEqual(len(ingest_failures), 1)

        # An index file with metadata that won't translate.
        metadata_failures[:] = []
        files = [os.path.join(INGESTDIR, "indexed_data", "metadata_subdir", "_index.json")]
        with self.assertRaises(RuntimeError):
            self.task.run(files, run=self.outputRun)
        self.assertEqual(len(metadata_failures), 2)

        # Bad index file.
        metadata_failures[:] = []
        files = [os.path.join(INGESTDIR, "indexed_data", "bad_index", "_index.json")]
        with self.assertRaises(RuntimeError):
            self.task.run(files, run=self.outputRun)
        self.assertEqual(len(metadata_failures), 1)

        # Ingest two files that have conflicting exposure metadata.
        ingest_failures[:] = []
        successes[:] = []
        # Ingest 4 files. 2 of them will implicitly find an index and one
        # will use a sidecar. The 4th will fail due to exposure conflict.
        files = [os.path.join(INGESTDIR, "indexed_data", f"dataset_{n}.yaml") for n in (1, 2, 3, 4)]
        new_run = self.outputRun + "_fail"
        with self.assertRaises(RuntimeError):
            self.task.run(files, run=new_run)
        self.assertEqual(len(ingest_failures), 1)
        self.assertEqual(len(successes), 3)


class TestRawIngestTaskPickle(unittest.TestCase):
    """Test that pickling of the RawIngestTask works properly."""

    @classmethod
    def setUpClass(cls):
        cls.root = tempfile.mkdtemp(dir=TESTDIR)
        cls.creatorButler = butlerTests.makeTestRepo(cls.root, {})

    @classmethod
    def tearDownClass(cls):
        if cls.root is not None:
            shutil.rmtree(cls.root, ignore_errors=True)

    def setUp(self):
        self.butler = butlerTests.makeTestCollection(self.creatorButler)

        self.config = RawIngestTask.ConfigClass()
        self.config.transfer = "copy"  # safe non-default value
        self.task = RawIngestTask(config=self.config, butler=self.butler)

    def testPickleTask(self):
        stream = pickle.dumps(self.task)
        copy = pickle.loads(stream)
        self.assertEqual(self.task.getFullName(), copy.getFullName())
        self.assertEqual(self.task.log.getName(), copy.log.getName())
        self.assertEqual(self.task.config, copy.config)
        self.assertEqual(self.task.butler._config, copy.butler._config)
        self.assertEqual(self.task.butler.collections, copy.butler.collections)
        self.assertEqual(self.task.butler.run, copy.butler.run)
        self.assertEqual(self.task.universe, copy.universe)
        self.assertEqual(self.task.datasetType, copy.datasetType)


if __name__ == "__main__":
    unittest.main()
