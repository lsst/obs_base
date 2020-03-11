# This file is part of obs_base.
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

"""Unit test base class for the gen2 to gen3 converter.
"""

import itertools
import shutil
import subprocess
import tempfile
import unittest


import lsst.afw.image
import lsst.afw.table
import lsst.daf.persistence
import lsst.daf.butler
import lsst.meas.algorithms
import lsst.utils.tests


class ConvertGen2To3TestCase:
    """Test the `convert_gen2_repo_to_gen3.py` script.

    Subclass this, and then `lsst.utils.tests.TestCase` and set the below
    attributes.
    """
    gen2root = ""
    """Root path to the gen2 repo to be converted."""

    gen2calib = ""
    """Path to the gen2 calib repo to be converted."""

    instrumentName = None
    """Name of the instrument for the gen3 registry, e.g. "DECam"."""

    instrumentClass = None
    """Full path to the `Instrument` class of the data to be converted, e.g.
    ``lsst.obs.decam.DarkEnergyCamera``."""

    config = None
    """Full path to a config override for ConvertRepoTask, to be applied after
    the Instrument overrides when running `convert_gen2_repo_to_gen3.py`."""

    biases = []
    """List dataIds to use to load gen3 biases to test that they exist."""

    biasName = "bias"
    """Name of the dataset that the biases are loaded into."""

    flats = []
    """List dataIds to use to load gen3 flats to test that they exist."""

    flatName = "flat"
    """Name of the dataset that the flats are loaded into."""

    darks = []
    """List dataIds to use to load gen3 darks to test that they exist."""

    darkName = "dark"
    """Name of the dataset that the darks are loaded into."""

    args = None
    """Other arguments to pass directly to the converter script, as a tuple."""

    refcats = []
    """Names of the reference catalogs to query for the existence of in the
    converted gen3 repo."""

    collections = set()
    """Additional collections that should appear in the gen3 repo, beyond the
    instrument name and any refcats, if ``refcats`` is non-empty above.
    Typically the only additional one necessary would be "skymaps"."""

    detectorKey = "ccd"
    """Key to use in a gen2 dataId to refer to a detector."""

    exposureKey = "visit"
    """Key to use in a gen2 dataId to refer to a visit or exposure."""

    def setUp(self):
        self.gen3root = tempfile.mkdtemp()
        self.gen2Butler = lsst.daf.persistence.Butler(root=self.gen2root, calibRoot=self.gen2calib)
        # This command is in obs_base, and we use the one that has been setup and scons'ed.
        self.cmd = "convert_gen2_repo_to_gen3.py"

        # if the collections set is empty we do not add to it, we create
        # a new instance version. Without this each subclass would add
        # to the same set.
        if not self.collections:
            self.collections = set()
        self.collections.add(self.instrumentName)
        if len(self.refcats) > 0:
            self.collections.add("refcats")

    def tearDown(self):
        shutil.rmtree(self.gen3root, ignore_errors=True)

    def _run_convert(self):
        """Convert a gen2 repo to gen3 for testing.
        """
        cmd = [self.cmd, self.instrumentClass,
               "--gen2root", self.gen2root,
               "--gen3root", self.gen3root,
               "--calibs", self.gen2calib
               ]
        if self.config is not None:
            cmd.extend(("--config", self.config))
        if self.args is not None:
            cmd.extend(self.args)
        print(f"Running command: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)

    def check_raw(self, gen3Butler, exposure, detector):
        """Check that a raw was converted correctly.

        Parameters
        ----------
        gen3Butler : `lsst.daf.butler.Butler`
            The Butler to be tested.
        exposure : `int`
            The exposure/vist identifier ``get`` from both butlers.
        detector : `int`
            The detector identifier to ``get`` from both butlers.
        """
        dataIdGen2 = {self.detectorKey: detector, self.exposureKey: exposure}
        try:
            gen2Exposure = self.gen2Butler.get("raw", dataId=dataIdGen2)
        except lsst.daf.persistence.butlerExceptions.NoResults:
            # ignore datasets that don't actually exist in the gen2 butler.
            return
        dataIdGen3 = dict(detector=detector, exposure=exposure, instrument=self.instrumentName)
        gen3Exposure = gen3Butler.get("raw", dataId=dataIdGen3)
        # Check that we got an Exposure, but not what type; there is
        # inconsistency between different obs packages.
        self.assertIsInstance(gen3Exposure, lsst.afw.image.Exposure)
        self.assertEqual(gen3Exposure.getInfo().getDetector().getId(), detector)
        self.assertMaskedImagesEqual(gen2Exposure.maskedImage, gen3Exposure.maskedImage)

    def check_calibs(self, calibName, calibIds, gen3Butler):
        """Test that we can get converted bias/dark/flat from the gen3 repo.

        Note: because there is no clear way to get calibrations from a gen2
        repo, we just test that the thing we got is an ExposureF here, and
        assume that formatter testing is handled properly elsewhere.

        Parameters
        ----------
        calibName : `str`
            The name of the calibration to attempt to get ("bias", "flat").
        calibIds : `list` of `dict`
            The list of calibration dataIds to get.
        gen3Butler : `lsst.daf.butler.Butler`
            The Butler to use to get the data.
        """
        for dataId in calibIds:
            datasets = list(gen3Butler.registry.queryDatasets(calibName, collections=..., dataId=dataId))
            gen3Exposure = gen3Butler.getDirect(datasets[0])
            self.assertIsInstance(gen3Exposure, lsst.afw.image.ExposureF)

    def check_defects(self, gen3Butler, detectors):
        """Test that we can get converted defects from the gen3 repo.

        Parameters
        ----------
        gen3Butler : `lsst.daf.butler.Butler`
            The Butler to be tested.
        detector : `int`
            The detector identifiers to ``get`` from the gen3 butler.
        """
        for detector in detectors:
            dataId = dict(detector=detector, instrument=self.instrumentName)
            # Fill out the missing parts of the dataId, as we don't a-priori
            # know e.g. the "calibration_label". Use the first element of the
            # result because we only need to check one.
            datasets = list(gen3Butler.registry.queryDatasets("defects", collections=..., dataId=dataId))
            if datasets:
                gen3Defects = gen3Butler.getDirect(datasets[0])
                self.assertIsInstance(gen3Defects, lsst.meas.algorithms.Defects)

    def check_refcat(self, gen3Butler):
        """Test that each expected refcat is in the gen3 repo.

        Parameters
        ----------
        gen3Butler : `lsst.daf.butler.Butler`
            The Butler to be tested.
        """
        if len(self.refcats) > 0:
            for refcat in self.refcats:
                query = gen3Butler.registry.queryDatasets(refcat, collections=["refcats"])
                self.assertGreater(len(list(query)), 0,
                                   msg=f"refcat={refcat} has no entries in collection 'refcats'.")

    def check_collections(self, gen3Butler):
        """Test that the correct set of collections is in the gen3 repo.

        Parameters
        ----------
        gen3Butler : `lsst.daf.butler.Butler`
            The Butler to be tested.
        """
        self.assertEqual(self.collections, gen3Butler.registry.getAllCollections())

    def test_convert(self):
        """Test that raws are converted correctly.
        """
        self._run_convert()
        gen3Butler = lsst.daf.butler.Butler(self.gen3root, run=self.instrumentName)
        self.check_collections(gen3Butler)

        # check every raw detector that the gen2 butler knows about
        detectors = self.gen2Butler.queryMetadata("raw", self.detectorKey)
        exposures = self.gen2Butler.queryMetadata("raw", self.exposureKey)
        for exposure, detector in itertools.product(exposures, detectors):
            self.check_raw(gen3Butler, exposure, detector)

        self.check_refcat(gen3Butler)
        self.check_defects(gen3Butler, detectors)
        self.check_calibs(self.biasName, self.biases, gen3Butler)
        self.check_calibs(self.flatName, self.flats, gen3Butler)
        self.check_calibs(self.darkName, self.darks, gen3Butler)


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
