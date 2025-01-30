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

"""Unit tests for daf_butler CLI ingest-raws command."""

import unittest

import lsst.obs.base
from lsst.daf.butler import Butler
from lsst.daf.butler.cli.butler import cli as butlerCli
from lsst.daf.butler.cli.utils import LogCliRunner, clickResultMsg
from lsst.daf.butler.tests import CliCmdTestBase
from lsst.obs.base.cli.cmd import ingest_raws
from lsst.obs.base.cli.cmd.commands import fits_re
from lsst.obs.base.ingest import RawIngestConfig


class IngestRawsTestCase(CliCmdTestBase, unittest.TestCase):
    """Test the ingest-raws command-line tool."""

    mockFuncName = "lsst.obs.base.cli.cmd.commands.script.ingestRaws"

    @staticmethod
    def defaultExpected():
        return dict(
            config={},
            config_file=None,
            ingest_task="lsst.obs.base.RawIngestTask",
            locations=(),
            output_run=None,
            processes=1,
            regex=fits_re,
            transfer="auto",
            track_file_attrs=True,
            update_records=False,
            fail_fast=False,
        )

    @staticmethod
    def command():
        return ingest_raws

    def test_repoAndOutput(self):
        """Test the most basic required arguments, repo and output run."""
        self.run_test(
            ["ingest-raws", "repo", "resources", "--output-run", "out"],
            self.makeExpected(repo="repo", locations=("resources",), output_run="out"),
        )

    def test_configMulti(self):
        """Test config overrides."""
        self.run_test(
            [
                "ingest-raws",
                "repo",
                "resources",
                "--output-run",
                "out",
                "-c",
                "foo=1",
                "--config",
                "bar=2",
                "--config",
                "baz=3",
            ],
            self.makeExpected(
                repo="repo",
                locations=("resources",),
                output_run="out",
                config={"foo": "1", "bar": "2", "baz": "3"},
            ),
        )

    def test_configFile(self):
        """Test config file override."""
        configFile = "path/to/file.txt"
        self.run_test(
            ["ingest-raws", "repo", "resources", "--output-run", "out", "--config-file", configFile],
            self.makeExpected(
                repo="repo", locations=("resources",), output_run="out", config_file=configFile
            ),
            withTempFile=configFile,
        )

    def test_transfer(self):
        """Test the transfer argument."""
        self.run_test(
            ["ingest-raws", "repo", "resources", "--output-run", "out", "--transfer", "symlink"],
            self.makeExpected(repo="repo", locations=("resources",), output_run="out", transfer="symlink"),
        )

    def test_ingestTask(self):
        """Test the ingest task argument."""
        self.run_test(
            ["ingest-raws", "repo", "resources", "--output-run", "out", "--ingest-task", "foo.bar.baz"],
            self.makeExpected(
                repo="repo", locations=("resources",), output_run="out", ingest_task="foo.bar.baz"
            ),
        )

    def test_locations(self):
        """Test that the locations argument accepts multiple inputs and splits
        commas.
        """
        self.run_test(
            ["ingest-raws", "repo", "in/directory/,in/another/dir/", "other/file.fits"],
            self.makeExpected(repo="repo", locations=("in/directory/", "in/another/dir/", "other/file.fits")),
        )


class PatchRawIngestTask(lsst.obs.base.RawIngestTask):
    """Ingest task with run() method disabled."""

    init_args = []

    def __init__(self, *args, **kwargs):
        self.init_args.append((args, kwargs))
        super().__init__(*args, **kwargs)

    def run(self, *args, **kwargs):
        pass


class RawIngestMockTest(unittest.TestCase):
    """Run ingest tests with mock."""

    def setUp(self):
        self.runner = LogCliRunner()

    def test(self):
        """Verify config gets applied properly."""
        with self.runner.isolated_filesystem():
            result = self.runner.invoke(butlerCli, ["create", "repo"])
            self.assertEqual(result.exit_code, 0, clickResultMsg(result))
            with unittest.mock.patch("lsst.obs.base.RawIngestTask", new=PatchRawIngestTask) as mock:
                # Call, override the name parameter of the config and set
                # fail-fast.
                result = self.runner.invoke(
                    butlerCli,
                    ["ingest-raws", "repo", "resources", "--config", "transfer=hardlink", "--fail-fast"],
                )
                self.assertEqual(result.exit_code, 0, clickResultMsg(result))
                # Verify the mock class was initialized exactly once:
                self.assertEqual(len(mock.init_args), 1)
                # Verify that the task was initialized with a 'butler' kwarg
                # that received a butler instance:
                self.assertIsInstance(mock.init_args[0][1]["butler"], Butler)
                expectedConfig = RawIngestConfig()
                # Verify that the task was initialized with a 'config' kwarg
                # that received an expected config:
                expectedConfig.update(transfer="hardlink")
                # Verfiy that --failfast caused the config's failFast
                # parameter to be set to True.
                expectedConfig.update(failFast=True)
                self.assertEqual(mock.init_args[0][1]["config"], expectedConfig)


if __name__ == "__main__":
    unittest.main()
