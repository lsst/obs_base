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


class ConvertTestCase(CliCmdTestBase, unittest.TestCase):

    mockFunc = "lsst.obs.base.cli.cmd.commands.script.convert"

    @staticmethod
    def defaultExpected():
        return dict(skymap_name=None,
                    skymap_config=None,
                    calibs=None,
                    reruns=(),
                    transfer="auto",
                    processes=1,
                    config_file=None)

    @staticmethod
    def command():
        return convert

    def test_repoInstrGen2root(self):
        """Test the most basic required arguments."""
        self.run_test(["convert", "here",
                       "--gen2root", "from"],
                      self.makeExpected(repo="here",
                                        gen2root="from"))

    def test_all(self):
        """Test all the arguments."""
        self.run_test(["convert", "here",
                       "--gen2root", "from",
                       "--skymap-name", "sky",
                       "--skymap-config", "/path/to/config",
                       "--calibs", "path/to/calib/repo",
                       "--reruns", "one,two",
                       "--reruns", "three",
                       "--transfer", "symlink",
                       "--processes", 1,
                       "--config-file", "/path/to/config"],
                      self.makeExpected(repo="here",
                                        gen2root="from",
                                        skymap_name="sky",
                                        skymap_config="/path/to/config",
                                        calibs="path/to/calib/repo",
                                        reruns=("one", "two", "three"),
                                        transfer="symlink",
                                        processes=1,
                                        config_file="/path/to/config"))

    def test_missing(self):
        """test a missing argument"""
        self.run_missing(["convert"], "Missing argument ['\"]REPO['\"]")
        self.run_missing(["convert", "here"], "Missing option ['\"]--gen2root['\"]")


if __name__ == "__main__":
    unittest.main()
