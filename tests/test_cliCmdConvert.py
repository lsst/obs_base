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

"""Unit tests for daf_butler CLI convert command.
"""

import unittest

from lsst.daf.butler.tests import CliCmdTestBase
from lsst.obs.base.cli.cmd import convert


class ConvertTestCase(CliCmdTestBase):

    defaultExpected = dict(skymap_name=None,
                           skymap_config=None,
                           calibs=None,
                           reruns=[],
                           transfer="auto",
                           config_file=None)

    command = convert

    def test_repoInstrGen2root(self):
        """Test the most basic required arguments."""
        self.run_test(["convert", "here",
                       "--gen2root", "from",
                       "--instrument", "a.b.c"],
                      self.makeExpected(repo="here",
                                        gen2root="from",
                                        instrument="a.b.c"))

    def test_all(self):
        """Test all the arguments."""
        self.run_test(["convert", "here",
                       "--gen2root", "from",
                       "--instrument", "a.b.c",
                       "--skymap-name", "sky",
                       "--skymap-config", "/path/to/config",
                       "--calibs", "path/to/calib/repo",
                       "--reruns", "one,two",
                       "--reruns", "three",
                       "--transfer", "symlink",
                       "--config-file", "/path/to/config"],
                      self.makeExpected(repo="here",
                                        gen2root="from",
                                        instrument="a.b.c",
                                        skymap_name="sky",
                                        skymap_config="/path/to/config",
                                        calibs="path/to/calib/repo",
                                        reruns=["one", "two", "three"],
                                        transfer="symlink",
                                        config_file="/path/to/config"))

    def test_missing(self):
        """test a missing argument"""
        self.run_missing(["convert"], 'Missing argument "REPO"')
        self.run_missing(["convert", "here", "--gen2root", "from"], 'Missing option "-i" / "--instrument"')
        self.run_missing(["convert", "here", "--gen2root", "from"], 'Missing option "-i" / "--instrument"')
        self.run_missing(["convert", "here", "--instrument", "instr"], 'Missing option "--gen2root"')

    def test_help(self):
        self.help_test()


if __name__ == "__main__":
    unittest.main()
