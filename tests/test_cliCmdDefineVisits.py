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

"""Unit tests for daf_butler CLI define-visits command."""

import unittest

from lsst.daf.butler.tests import CliCmdTestBase
from lsst.obs.base.cli.cmd import define_visits


class DefineVisitsTest(CliCmdTestBase, unittest.TestCase):
    """Test the define-visits command-line tool."""

    mockFuncName = "lsst.obs.base.cli.cmd.commands.script.defineVisits"

    @staticmethod
    def defaultExpected():
        return {"config_file": None, "collections": ()}

    @staticmethod
    def command():
        return define_visits

    def test_repoBasic(self):
        """Test the most basic required arguments."""
        self.run_test(
            ["define-visits", "here", "a.b.c"],
            self.makeExpected(
                repo="here", instrument="a.b.c", where="", update_records=False, incremental=False
            ),
        )

    def test_all(self):
        """Test all the arguments."""
        self.run_test(
            [
                "define-visits",
                "here",
                "a.b.c",
                "--collections",
                "foo/bar,baz",
                "--config-file",
                "/path/to/config",
                "--collections",
                "boz",
                "--where",
                "x=y",
                "--incremental",
                "--update-records",
            ],
            self.makeExpected(
                repo="here",
                instrument="a.b.c",
                config_file="/path/to/config",
                # The list of collections must be in
                # exactly the same order as it is
                # passed in the list of arguments to
                # run_test.
                collections=("foo/bar", "baz", "boz"),
                where="x=y",
                update_records=True,
                incremental=True,
            ),
        )

    def test_missing(self):
        """Test a missing argument."""
        self.run_missing(["define-visits"], "Missing argument ['\"]REPO['\"]")
        self.run_missing(["define-visits", "here"], "Missing argument ['\"]INSTRUMENT['\"]")


if __name__ == "__main__":
    unittest.main()
