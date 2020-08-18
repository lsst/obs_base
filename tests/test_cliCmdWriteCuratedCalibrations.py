# This file is part of daf_butler.
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

"""Unit tests for daf_butler CLI define-visits command.
"""

import unittest

from lsst.daf.butler.tests import CliCmdTestBase
from lsst.obs.base.cli.cmd import write_curated_calibrations


class WriteCuratedCalibrationsTest(CliCmdTestBase, unittest.TestCase):

    @staticmethod
    def defaultExpected():
        return dict()

    @staticmethod
    def command():
        return write_curated_calibrations

    def test_repoBasic(self):
        """Test the most basic required arguments."""
        self.run_test(["write-curated-calibrations", "here",
                       "--instrument", "a.b.c",
                       "--output-run", "foo"],
                      self.makeExpected(repo="here",
                                        instrument="a.b.c",
                                        output_run="foo"))

    def test_missing(self):
        """test a missing argument"""
        self.run_missing(["write-curated-calibrations"], 'Missing argument "REPO"')
        self.run_missing(["write-curated-calibrations", "here"], 'Missing option "-i" / "--instrument"')


if __name__ == "__main__":
    unittest.main()
