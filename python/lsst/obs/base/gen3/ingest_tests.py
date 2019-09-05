# This file is part of obs_subaru.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
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

"""Base class for writing Gen3 raw data ingest tests.
"""

__all__ = ("IngestTestBase",)

import abc
import tempfile
import unittest
import os
import shutil

from lsst.daf.butler import Butler
from lsst.obs.base.gen3 import RawIngestTask


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

    dataId = {}
    """Butler data ID of a file to ingest when testing."""

    file = ""
    """Full path to a file to ingest in tests."""

    def setUp(self):
        # Use a temporary working directory
        self.root = tempfile.mkdtemp(dir=self.ingestDir)
        Butler.makeRepo(self.root)
        self.butler = Butler(self.root, run="raw")

        # Register the instrument and its static metadata
        self.instrument.register(self.butler.registry)

        # Make a default config for test methods to play with
        self.config = RawIngestTask.ConfigClass()
        self.config.onError = "break"

    def tearDown(self):
        if os.path.exists(self.root):
            shutil.rmtree(self.root, ignore_errors=True)

    def runIngest(self, files=None):
        """
        Initialize and run RawIngestTask on a list of files.

        Parameters
        ----------
        files : `list`, [`str`], or None
            List of files to be ingested, or None to use ``self.file``
        """
        if files is None:
            files = [self.file]
        task = RawIngestTask(config=self.config, butler=self.butler)
        task.log.setLevel(task.log.FATAL)  # silence logs, since we expect a lot of warnings
        task.run(files)

    def runIngestTest(self, files=None):
        """
        Test that RawIngestTask ingested the expected files.

        Parameters
        ----------
        files : `list`, [`str`], or None
            List of files to be ingested, or None to use ``self.file``
        """
        self.runIngest(files)
        exposure = self.butler.get("raw", self.dataId)
        metadata = self.butler.get("raw.metadata", self.dataId)
        image = self.butler.get("raw.image", self.dataId)
        self.assertImagesEqual(exposure.image, image)
        self.assertEqual(metadata.toDict(), exposure.getMetadata().toDict())

    def testSymLink(self):
        self.config.transfer = "symlink"
        self.runIngestTest()

    def testCopy(self):
        self.config.transfer = "copy"
        self.runIngestTest()

    def testHardLink(self):
        self.config.transfer = "hardlink"
        try:
            self.runIngestTest()
        except PermissionError as err:
            raise unittest.SkipTest("Skipping hard-link test because input data"
                                    " is on a different filesystem.") from err

    def testInPlace(self):
        """Test that files already in the directory can be added to the
        registry in-place.
        """
        # copy into repo root manually
        newPath = os.path.join(self.butler.datastore.root, os.path.basename(self.file))
        shutil.copyfile(self.file, newPath)
        self.config.transfer = None
        self.runIngestTest([newPath])

    def testOnConflictFail(self):
        """Re-ingesting the same data into the repository should fail, if
        configured to do so.
        """
        self.config.transfer = "symlink"
        self.config.conflict = "fail"
        self.runIngest()
        with self.assertRaises(Exception):
            self.runIngest()

    def testOnConflictIgnore(self):
        """Re-ingesting the same data into the repository does not fail, if
        configured to ignore conflict errors.
        """
        self.config.transfer = "symlink"
        self.config.conflict = "ignore"
        self.runIngest()   # this one should succeed
        n1, = self.butler.registry.query("SELECT COUNT(*) FROM Dataset")
        self.runIngest()   # this one should silently fail
        n2, = self.butler.registry.query("SELECT COUNT(*) FROM Dataset")
        self.assertEqual(n1, n2)

    def testOnConflictStash(self):
        """Re-ingesting the same data will be put into a different collection,
        if configured to do so.
        """
        self.config.transfer = "symlink"
        self.config.conflict = "ignore"
        self.config.stash = "stash"
        self.runIngest()   # this one should write to 'raw'
        self.runIngest()   # this one should write to 'stash'
        dt = self.butler.registry.getDatasetType("raw.metadata")
        ref1 = self.butler.registry.find(self.butler.collection, dt, self.dataId)
        ref2 = self.butler.registry.find("stash", dt, self.dataId)
        self.assertNotEqual(ref1.id, ref2.id)
        self.assertEqual(self.butler.get(ref1).toDict(), self.butler.getDirect(ref2).toDict())

    def testOnErrorBreak(self):
        """Test that errors do not roll back success, when configured to do so.

        Failing to ingest a nonexistent file after ingesting the valid one should
        leave the valid one in the registry, despite raising an exception.
        """
        self.config.transfer = "symlink"
        self.config.onError = "break"
        with self.assertRaises(Exception):
            self.runIngest(files=[self.file, "nonexistent.fits"])
        dt = self.butler.registry.getDatasetType("raw.metadata")
        self.assertIsNotNone(self.butler.registry.find(self.butler.collection, dt, self.dataId))

    def testOnErrorContinue(self):
        """Failing to ingest nonexistent files before and after ingesting the
        valid one should leave the valid one in the registry and not raise
        an exception.
        """
        self.config.transfer = "symlink"
        self.config.onError = "continue"
        self.runIngest(files=["nonexistent.fits", self.file, "still-not-here.fits"])
        dt = self.butler.registry.getDatasetType("raw.metadata")
        self.assertIsNotNone(self.butler.registry.find(self.butler.collection, dt, self.dataId))

    def testOnErrorRollback(self):
        """Failing to ingest nonexistent files after ingesting the
        valid one should leave the registry unchanged.
        """
        self.config.transfer = "symlink"
        self.config.onError = "rollback"
        with self.assertRaises(Exception):
            self.runIngest(file=[self.file, "nonexistent.fits"])
        try:
            dt = self.butler.registry.getDatasetType("raw.metadata")
        except KeyError:
            # If we also rollback registering the DatasetType, that's fine,
            # but not required.
            pass
        else:
            self.assertIsNotNone(self.butler.registry.find(self.butler.collection, dt, self.dataId))
