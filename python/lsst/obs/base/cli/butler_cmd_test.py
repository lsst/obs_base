# This file is part of obs_base.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
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
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Base class for writing CLI butler command tests.
"""

__all__ = ("ButlerCmdTestBase",)

import abc
import click
import click.testing

from lsst.daf.butler.cli import butler


class ButlerCmdTestBase(metaclass=abc.ABCMeta):
    """Base class for tests of butler command line interface subcommands.
    Subclass from this, then `unittest.TestCase` to get a working test suite.
    """

    instrument_class = None
    """The fully qualified instrument class.
    """

    instrument_name = None
    """The instrument name."""

    def test_cli(self):
        runner = click.testing.CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(butler.cli, ["create", "here"])
            self.assertEqual(result.exit_code, 0, result.output)
            result = runner.invoke(butler.cli, ["register-instrument",
                                                "here",
                                                "-i", self.instrument_class])
            self.assertEqual(result.exit_code, 0, result.output)
            result = runner.invoke(butler.cli, ["write-curated-calibrations",
                                                "here",
                                                "-i", self.instrument_name,
                                                "--output-run", "calib/hsc"])
            self.assertEqual(result.exit_code, 0, result.output)
