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

"""Unit tests for daf_butler CLI ingest-raws command.
"""

import unittest

from lsst.daf.butler.tests import CliCmdTestBase
from lsst.obs.base.cli.cmd import ingest_raws


class IngestRawsTestCase(CliCmdTestBase, unittest.TestCase):

    @staticmethod
    def defaultExpected():
        return dict(directory=None,
                    file=None,
                    transfer="auto",
                    ingest_task="lsst.obs.base.RawIngestTask",
                    config={},
                    config_file=None)

    @staticmethod
    def command():
        return ingest_raws

    def test_repoAndOutput(self):
        """Test the most basic required arguments, repo and output run"""
        self.run_test(["ingest-raws", "here",
                       "--output-run", "out"],
                      self.makeExpected(repo="here", output_run="out"))

    def test_configMulti(self):
        """Test config overrides"""
        self.run_test(["ingest-raws", "here",
                       "--output-run", "out",
                       "-c", "foo=1",
                       "--config", "bar=2,baz=3"],
                      self.makeExpected(repo="here",
                                        output_run="out",
                                        config=dict(foo="1", bar="2", baz="3")))

    def test_configFile(self):
        """Test config file override"""
        configFile = "path/to/file.txt"
        self.run_test(["ingest-raws", "here",
                       "--output-run", "out",
                       "--config-file", configFile],
                      self.makeExpected(repo="here",
                                        output_run="out",
                                        config_file=configFile),
                      withTempFile=configFile)

    def test_transfer(self):
        """Test the transfer argument"""
        self.run_test(["ingest-raws", "here",
                       "--output-run", "out",
                       "--transfer", "symlink"],
                      self.makeExpected(repo="here",
                                        output_run="out",
                                        transfer="symlink"))

    def test_ingestTask(self):
        """Test the ingest task argument"""
        self.run_test(["ingest-raws", "here",
                       "--output-run", "out",
                       "--ingest-task", "foo.bar.baz"],
                      self.makeExpected(repo="here", output_run="out", ingest_task="foo.bar.baz"))


if __name__ == "__main__":
    unittest.main()
