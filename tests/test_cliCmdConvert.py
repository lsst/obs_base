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

import click
import click.testing
import unittest

from lsst.daf.butler.cli import butler
from lsst.daf.butler.cli.utils import Mocker, mockEnvVar


def makeExpectedKwargs(**kwargs):
    expected = dict(skymap_name=None,
                    skymap_config=None,
                    calibs=None,
                    reruns=[],
                    transfer="auto",
                    config_file=None)
    expected.update(kwargs)
    return expected


class Case(unittest.TestCase):

    def run_test(self, inputs, expectedKwargs):
        """Test command line interaction with convert command function.

        Parameters
        ----------
        inputs : [`str`]
            A list of the arguments to the butler command, starting with
            `ingest-raws`
        expectedKwargs : `dict` [`str`, `str`]
            The expected arguments to the ingestRaws command function, keys are
            the argument name and values are the argument value.
        """
        runner = click.testing.CliRunner(env=mockEnvVar)
        result = runner.invoke(butler.cli, inputs)
        self.assertEqual(result.exit_code, 0, f"output: {result.output} exception: {result.exception}")
        Mocker.mock.assert_called_with(**expectedKwargs)

    def test_repoInstrGen2root(self):
        """Test the most basic required arguments."""
        expected = makeExpectedKwargs(repo="here", gen2root="from", instrument="a.b.c")
        self.run_test(["convert", "here",
                       "--gen2root", "from",
                       "--instrument", "a.b.c"], expected)

    def test_all(self):
        """Test all the arguments."""
        expected = dict(repo="here", gen2root="from", instrument="a.b.c", skymap_name="sky",
                        skymap_config="/path/to/config", calibs="path/to/calib/repo",
                        reruns=["one", "two", "three"], transfer="symlink", config_file="/path/to/config")
        self.run_test(["convert", "here",
                       "--gen2root", "from",
                       "--instrument", "a.b.c",
                       "--skymap-name", "sky",
                       "--skymap-config", "/path/to/config",
                       "--calibs", "path/to/calib/repo",
                       "--reruns", "one,two",
                       "--reruns", "three",
                       "--transfer", "symlink",
                       "--config-file", "/path/to/config"], expected)

    def test_missing(self):
        """test a missing argument"""
        with self.assertRaises(TypeError):
            self.run_test(["convert", "--gen2root", "from"])


if __name__ == "__main__":
    unittest.main()
