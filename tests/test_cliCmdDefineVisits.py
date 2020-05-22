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

import click
import click.testing
import unittest

from lsst.daf.butler.cli import butler
from lsst.daf.butler.cli.utils import Mocker, mockEnvVar


def makeExpectedKwargs(**kwargs):
    expected = dict(config_file=None,
                    collections=[])
    expected.update(kwargs)
    return expected


class Case(unittest.TestCase):

    def run_test(self, inputs, expectedKwargs):
        """Test command line interaction with the defineVisits command function.

        Parameters
        ----------
        inputs : [`str`]
            A list of the arguments to the butler command, starting with
            `define-visits`
        expectedKwargs : `dict` [`str`, `str`]
            The expected arguments to the define-visits command function, keys are
            the argument name and values are the argument value.
        """
        runner = click.testing.CliRunner(env=mockEnvVar)
        result = runner.invoke(butler.cli, inputs)
        self.assertEqual(result.exit_code, 0, f"output: {result.output} exception: {result.exception}")
        Mocker.mock.assert_called_with(**expectedKwargs)

    def test_repoBasic(self):
        """Test the most basic required arguments."""
        expected = makeExpectedKwargs(repo="here", instrument="a.b.c")
        self.run_test(["define-visits", "here",
                       "--instrument", "a.b.c"], expected)

    def test_all(self):
        """Test all the arguments."""
        expected = dict(repo="here", instrument="a.b.c", config_file="/path/to/config",
                        # The list of collections must be in exactly the same order as it is passed in the
                        # list of arguments to run_test.
                        collections=["foo/bar", "baz", "boz"])
        self.run_test(["define-visits", "here",
                       "--instrument", "a.b.c",
                       "--collections", "foo/bar,baz",
                       "--config-file", "/path/to/config",
                       "--collections", "boz"], expected)

    def test_missing(self):
        """test a missing argument"""
        with self.assertRaises(TypeError):
            self.run_test(["define-visits", "from"])


if __name__ == "__main__":
    unittest.main()
