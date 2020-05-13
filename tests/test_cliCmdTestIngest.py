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

"""Unit tests for daf_butler CLI config-dump command.
"""

import click
import click.testing
import unittest

from lsst.daf.butler.cli import butler
from lsst.daf.butler.cli.utils import Mocker, mockEnvVar


def makeExpectedKwargs(**kwargs):
    expected = dict(directory=None,
                    file=None,
                    transfer="auto",
                    ingest_task="lsst.obs.base.RawIngestTask",
                    config={},
                    config_file=None)
    expected.update(kwargs)
    return expected


class Suite(unittest.TestCase):

    def run_test(self, inputs, expectedKwargs):
        """Test command line interaction with ingest_raws command function.

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

    def test_repoAndOutput(self):
        """Test the most basic required arguments, repo and output run"""
        expected = makeExpectedKwargs(repo="here", output_run="out")
        self.run_test(["ingest-raws", "here",
                       "--output-run", "out"], expected)

    def test_configMulti(self):
        """Test config overrides"""
        expected = makeExpectedKwargs(repo="here", output_run="out", config=dict(foo="1", bar="2", baz="3"))
        self.run_test(["ingest-raws", "here",
                       "--output-run", "out",
                       "-c", "foo=1",
                       "--config", "bar=2,baz=3"], expected)

    def test_configFile(self):
        """Test config file override"""
        expected = makeExpectedKwargs(repo="here", output_run="out", config_file="path/to/file.txt")
        self.run_test(["ingest-raws", "here",
                       "--output-run", "out",
                       "--config-file", "path/to/file.txt"], expected)

    def test_transfer(self):
        """Test the transfer argument"""
        expected = makeExpectedKwargs(repo="here", output_run="out", transfer="foo")
        self.run_test(["ingest-raws", "here",
                       "--output-run", "out",
                       "--transfer", "foo"], expected)

    def test_ingestTask(self):
        """Test the ingest task argument"""
        expected = makeExpectedKwargs(repo="here", output_run="out", ingest_task="foo.bar.baz")
        self.run_test(["ingest-raws", "here",
                       "--output-run", "out",
                       "--ingest-task", "foo.bar.baz"], expected)


if __name__ == "__main__":
    unittest.main()
