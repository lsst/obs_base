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

from lsst.daf.butler.cli import butler
from lsst.daf.butler.cli.utils import LogCliRunner
from lsst.utils import doImport


class ButlerCmdTestBase(metaclass=abc.ABCMeta):
    """Base class for tests of butler command line interface subcommands.
    Subclass from this, then `unittest.TestCase` to get a working test suite.
    """

    @property
    @abc.abstractmethod
    def instrumentClassName(self):
        """The fully qualified instrument class name.

        Returns
        -------
        `str`
            The fully qualified instrument class name.
        """
        pass

    @property
    def secondInstrumentClassName(self):
        """Optional; if provided the register-instrument test will try to
        register two instruments.

        Returns
        -------
        `str` or `None`
            The fully qualified instrument class name.
        """

    @property
    def instrument(self):
        """The instrument class."""
        return doImport(self.instrumentClassName)

    @property
    def instrumentName(self):
        """The name of the instrument.

        Returns
        -------
        `str`
            The name of the instrument.
        """
        return self.instrument.getName()

    def test_cli(self):
        runner = LogCliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(butler.cli, ["create", "here"])
            self.assertEqual(result.exit_code, 0, f"output: {result.output} exception: {result.exception}")
            registerInstrumentArgs = ["register-instrument", "here", self.instrumentClassName]
            if self.secondInstrumentClassName is not None:
                registerInstrumentArgs.append(self.secondInstrumentClassName)
            result = runner.invoke(butler.cli, registerInstrumentArgs)
            self.assertEqual(result.exit_code, 0, f"output: {result.output} exception: {result.exception}")
            result = runner.invoke(butler.cli, ["write-curated-calibrations",
                                                "here",
                                                "--instrument", self.instrumentName,
                                                "--collection", "collection"])
            self.assertEqual(result.exit_code, 0, f"output: {result.output} exception: {result.exception}")
