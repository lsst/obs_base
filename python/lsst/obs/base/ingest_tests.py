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

import lsst.afw.cameraGeom
from lsst.daf.butler import Butler
from lsst.daf.butler.cli.butler import cli as butlerCli
from lsst.daf.butler.cli.utils import LogCliRunner
import lsst.obs.base
from lsst.utils import doImport
from .utils import getInstrument
from . import script


class IngestTestBase(metaclass=abc.ABCMeta):
    """Base class for tests of gen3 ingest. Subclass from this, then
    `unittest.TestCase` to get a working test suite.
    """

    ingestDir = ""
    """Root path to ingest files into. Typically `obs_package/tests/`; the
    actual directory will be a tempdir under this one.
    """

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

    defineVisitsTask = lsst.obs.base.DefineVisitsTask
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

    outputRun = "raw"
    """The name of the output run to use in tests.
    """

    @property
    @abc.abstractmethod
    def instrumentClassName(self):
        """The fully qualified instrument class name.

        Returns
        -------
        `str`
            The fully qualified instrument class name.
        """
        pass

    @property
    def instrumentClass(self):
        """The instrument class."""
        return doImport(self.instrumentClassName)

    @property
    def instrumentName(self):
        """The name of the instrument.

        Returns
        -------
        `str`
            The name of the instrument.
        """
        return self.instrumentClass.getName()

    def setUp(self):
        # Use a temporary working directory
        self.root = tempfile.mkdtemp(dir=self.ingestDir)
        self._createRepo()

        # Register the instrument and its static metadata
        self._registerInstrument()

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

    def _createRepo(self):
        """Use the Click `testing` module to call the butler command line api
        to create a repository."""
        runner = LogCliRunner()
        result = runner.invoke(butlerCli, ["create", self.root])
        self.assertEqual(result.exit_code, 0, f"output: {result.output} exception: {result.exception}")

    def _ingestRaws(self, transfer):
        """Use the Click `testing` module to call the butler command line api
        to ingest raws.

        Parameters
        ----------
        transfer : `str`
            The external data transfer type.
        """
        runner = LogCliRunner()
        result = runner.invoke(butlerCli, ["ingest-raws", self.root, self.file,
                                           "--output-run", self.outputRun,
                                           "--transfer", transfer,
                                           "--ingest-task", self.rawIngestTask])
        self.assertEqual(result.exit_code, 0, f"output: {result.output} exception: {result.exception}")

    def _registerInstrument(self):
        """Use the Click `testing` module to call the butler command line api
        to register the instrument."""
        runner = LogCliRunner()
        result = runner.invoke(butlerCli, ["register-instrument", self.root, self.instrumentClassName])
        self.assertEqual(result.exit_code, 0, f"output: {result.output} exception: {result.exception}")

    def _writeCuratedCalibrations(self):
        """Use the Click `testing` module to call the butler command line api
        to write curated calibrations."""
        runner = LogCliRunner()
        result = runner.invoke(butlerCli, ["write-curated-calibrations", self.root, self.instrumentName])
        self.assertEqual(result.exit_code, 0, f"output: {result.output} exception: {result.exception}")

    def testLink(self):
        self._ingestRaws(transfer="link")
        self.verifyIngest()

    def testSymLink(self):
        self._ingestRaws(transfer="symlink")
        self.verifyIngest()

    def testCopy(self):
        self._ingestRaws(transfer="copy")
        self.verifyIngest()

    def testHardLink(self):
        try:
            self._ingestRaws(transfer="hardlink")
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
        newPath = butler.datastore.root.join(os.path.basename(self.file))
        os.symlink(os.path.abspath(self.file), newPath.ospath)
        self._ingestRaws(transfer=None)
        self.verifyIngest()

    def testFailOnConflict(self):
        """Re-ingesting the same data into the repository should fail.
        """
        self._ingestRaws(transfer="symlink")
        with self.assertRaises(Exception):
            self._ingestRaws(transfer="symlink")

    def testWriteCuratedCalibrations(self):
        """Test that we can ingest the curated calibrations, and read them
        with `loadCamera` both before and after.
        """
        if self.curatedCalibrationDatasetTypes is None:
            raise unittest.SkipTest("Class requests disabling of writeCuratedCalibrations test")

        butler = Butler(self.root, writeable=False)
        collection = self.instrumentClass.makeCalibrationCollectionName()

        # Trying to load a camera with a data ID not known to the registry
        # is an error, because we can't get any temporal information.
        with self.assertRaises(LookupError):
            lsst.obs.base.loadCamera(butler, self.dataIds[0], collections=collection)

        # Ingest raws in order to get some exposure records.
        self._ingestRaws(transfer="auto")

        # Load camera should returned an unversioned camera because there's
        # nothing in the repo.
        camera, isVersioned = lsst.obs.base.loadCamera(butler, self.dataIds[0], collections=collection)
        self.assertFalse(isVersioned)
        self.assertIsInstance(camera, lsst.afw.cameraGeom.Camera)

        self._writeCuratedCalibrations()

        # Make a new butler instance to make sure we don't have any stale
        # caches (e.g. of DatasetTypes).  Note that we didn't give
        # _writeCuratedCalibrations the butler instance we had, because it's
        # trying to test the CLI interface anyway.
        butler = Butler(self.root, writeable=False)

        for datasetTypeName in self.curatedCalibrationDatasetTypes:
            with self.subTest(dtype=datasetTypeName):
                found = list(
                    butler.registry.queryDatasetAssociations(
                        datasetTypeName,
                        collections=collection,
                    )
                )
                self.assertGreater(len(found), 0, f"Checking {datasetTypeName}")

        # Load camera should returned the versioned camera from the repo.
        camera, isVersioned = lsst.obs.base.loadCamera(butler, self.dataIds[0], collections=collection)
        self.assertTrue(isVersioned)
        self.assertIsInstance(camera, lsst.afw.cameraGeom.Camera)

    def testDefineVisits(self):
        if self.visits is None:
            self.skipTest("Expected visits were not defined.")
        self._ingestRaws(transfer="link")

        # Calling defineVisits tests the implementation of the butler command line interface "define-visits"
        # subcommand. Functions in the script folder are generally considered protected and should not be used
        # as public api.
        script.defineVisits(self.root, config_file=None, collections=self.outputRun,
                            instrument=self.instrumentName)

        # Test that we got the visits we expected.
        butler = Butler(self.root, run=self.outputRun)
        visits = butler.registry.queryDataIds(["visit"]).expanded().toSet()
        self.assertCountEqual(visits, self.visits.keys())
        instr = getInstrument(self.instrumentName, butler.registry)
        camera = instr.getCamera()
        for foundVisit, (expectedVisit, expectedExposures) in zip(visits, self.visits.items()):
            # Test that this visit is associated with the expected exposures.
            foundExposures = butler.registry.queryDataIds(["exposure"], dataId=expectedVisit
                                                          ).expanded().toSet()
            self.assertCountEqual(foundExposures, expectedExposures)
            # Test that we have a visit region, and that it contains all of the
            # detector+visit regions.
            self.assertIsNotNone(foundVisit.region)
            detectorVisitDataIds = butler.registry.queryDataIds(["visit", "detector"], dataId=expectedVisit
                                                                ).expanded().toSet()
            self.assertEqual(len(detectorVisitDataIds), len(camera))
            for dataId in detectorVisitDataIds:
                self.assertTrue(foundVisit.region.contains(dataId.region))
