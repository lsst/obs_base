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
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Base class for writing Gen3 raw data ingest tests.
"""

__all__ = ("IngestTestBase",)

import abc
import tempfile
import unittest
import os
import shutil

from lsst.daf.butler import Butler
from lsst.daf.butler.script import createRepo
import lsst.obs.base
from .utils import getInstrument
from .script import ingestRaws, registerInstrument, writeCuratedCalibrations


class IngestTestBase(metaclass=abc.ABCMeta):
    """Base class for tests of gen3 ingest. Subclass from this, then
    `unittest.TestCase` to get a working test suite.
    """

    ingestDir = ""
    """Root path to ingest files into. Typically `obs_package/tests/`; the
    actual directory will be a tempdir under this one.
    """

    instrument = None
    """The instrument to be registered and tested."""

    dataIds = []
    """list of butler data IDs of files that should have been ingested."""

    file = ""
    """Full path to a file to ingest in tests."""

    rawIngestTask = "lsst.obs.base.RawIngestTask"
    """The task to use in the Ingest test."""

    curatedCalibrationDatasetTypes = None
    """List or tuple of Datasets types that should be present after calling
    writeCuratedCalibrations. If `None` writeCuratedCalibrations will
    not be called and the test will be skipped."""

    DefineVisitsTask = lsst.obs.base.DefineVisitsTask
    """The task to use to define visits from groups of exposures.

    This is ignored if ``visits`` is `None`.
    """

    visits = {}
    """A dictionary mapping visit data IDs the lists of exposure data IDs that
    are associated with them.

    If this is empty (but not `None`), visit definition will be run but no
    visits will be expected (e.g. because no exposures are on-sky
    observations).
    """

    instrument = ""
    """The fully qualified name of the instrument.
    """

    instrumentName = ""
    """The name of the instrument.
    """

    outputRun = "raw"
    """The name of the output run to use in tests.
    """

    def setUp(self):
        # Use a temporary working directory
        self.root = tempfile.mkdtemp(dir=self.ingestDir)
        createRepo(self.root)

        # Register the instrument and its static metadata
        registerInstrument(self.root, self.instrument)

    def tearDown(self):
        if os.path.exists(self.root):
            shutil.rmtree(self.root, ignore_errors=True)

    def verifyIngest(self, files=None, cli=False):
        """
        Test that RawIngestTask ingested the expected files.

        Parameters
        ----------
        files : `list` [`str`], or None
            List of files to be ingested, or None to use ``self.file``
        """
        butler = Butler(self.root, run=self.outputRun)
        datasets = butler.registry.queryDatasets(self.outputRun, collections=...)
        self.assertEqual(len(list(datasets)), len(self.dataIds))
        for dataId in self.dataIds:
            exposure = butler.get(self.outputRun, dataId)
            metadata = butler.get("raw.metadata", dataId)
            self.assertEqual(metadata.toDict(), exposure.getMetadata().toDict())

            # Since components follow a different code path we check that
            # WCS match and also we check that at least the shape
            # of the image is the same (rather than doing per-pixel equality)
            wcs = butler.get("raw.wcs", dataId)
            self.assertEqual(wcs, exposure.getWcs())

            rawImage = butler.get("raw.image", dataId)
            self.assertEqual(rawImage.getBBox(), exposure.getBBox())

            self.checkRepo(files=files)

    def checkRepo(self, files=None):
        """Check the state of the repository after ingest.

        This is an optional hook provided for subclasses; by default it does
        nothing.

        Parameters
        ----------
        files : `list` [`str`], or None
            List of files to be ingested, or None to use ``self.file``
        """
        pass

    def testLink(self):
        ingestRaws(self.root, self.outputRun, file=self.file, transfer="link",
                   ingest_task=self.rawIngestTask)
        self.verifyIngest()

    def testSymLink(self):
        ingestRaws(self.root, self.outputRun, file=self.file, transfer="symlink",
                   ingest_task=self.rawIngestTask)
        self.verifyIngest()

    def testCopy(self):
        ingestRaws(self.root, self.outputRun, file=self.file, transfer="copy",
                   ingest_task=self.rawIngestTask)
        self.verifyIngest()

    def testHardLink(self):
        try:
            ingestRaws(self.root, self.outputRun, file=self.file, transfer="hardlink",
                       ingest_task=self.rawIngestTask)
            self.verifyIngest()
        except PermissionError as err:
            raise unittest.SkipTest("Skipping hard-link test because input data"
                                    " is on a different filesystem.") from err

    def testInPlace(self):
        """Test that files already in the directory can be added to the
        registry in-place.
        """
        # symlink into repo root manually
        butler = Butler(self.root, run=self.outputRun)
        newPath = os.path.join(butler.datastore.root, os.path.basename(self.file))
        os.symlink(os.path.abspath(self.file), newPath)
        ingestRaws(self.root, self.outputRun, file=newPath, transfer=None, ingest_task=self.rawIngestTask)
        self.verifyIngest()

    def testFailOnConflict(self):
        """Re-ingesting the same data into the repository should fail.
        """
        ingestRaws(self.root, self.outputRun, file=self.file, transfer="symlink",
                   ingest_task=self.rawIngestTask)
        with self.assertRaises(Exception):
            ingestRaws(self.root, self.outputRun, file=self.file, transfer="symlink",
                       ingest_task=self.rawIngestTask)

    def testWriteCuratedCalibrations(self):
        """Test that we can ingest the curated calibrations"""
        if self.curatedCalibrationDatasetTypes is None:
            raise unittest.SkipTest("Class requests disabling of writeCuratedCalibrations test")

        writeCuratedCalibrations(self.root, self.instrumentName, self.outputRun)

        dataId = {"instrument": self.instrumentName}
        butler = Butler(self.root, run=self.outputRun)
        for datasetTypeName in self.curatedCalibrationDatasetTypes:
            with self.subTest(dtype=datasetTypeName, dataId=dataId):
                datasets = list(butler.registry.queryDatasets(datasetTypeName, collections=...,
                                                              dataId=dataId))
                self.assertGreater(len(datasets), 0, f"Checking {datasetTypeName}")

    def testDefineVisits(self):
        if self.visits is None:
            self.skipTest("Expected visits were not defined.")
        ingestRaws(self.root, self.outputRun, file=self.file, transfer="link",
                   ingest_task=self.rawIngestTask)

        config = self.DefineVisitsTask.ConfigClass()
        butler = Butler(self.root, run=self.outputRun)
        instrument = getInstrument(self.instrumentName, butler.registry)
        instrument.applyConfigOverrides(self.DefineVisitsTask._DefaultName, config)
        task = self.DefineVisitsTask(config=config, butler=butler)
        task.run(self.dataIds)

        # Test that we got the visits we expected.
        visits = set(butler.registry.queryDimensions(["visit"], expand=True))
        self.assertCountEqual(visits, self.visits.keys())
        camera = instrument.getCamera()
        for foundVisit, (expectedVisit, expectedExposures) in zip(visits, self.visits.items()):
            # Test that this visit is associated with the expected exposures.
            foundExposures = set(butler.registry.queryDimensions(["exposure"], dataId=expectedVisit,
                                                                 expand=True))
            self.assertCountEqual(foundExposures, expectedExposures)
            # Test that we have a visit region, and that it contains all of the
            # detector+visit regions.
            self.assertIsNotNone(foundVisit.region)
            detectorVisitDataIds = set(butler.registry.queryDimensions(["visit", "detector"],
                                                                       dataId=expectedVisit,
                                                                       expand=True))
            self.assertEqual(len(detectorVisitDataIds), len(camera))
            for dataId in detectorVisitDataIds:
                self.assertTrue(foundVisit.region.contains(dataId.region))
