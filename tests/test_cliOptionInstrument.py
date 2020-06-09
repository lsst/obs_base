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

"""Unit tests for the butler instrument CLI option.
"""

import click
import click.testing
import unittest

from lsst.daf.butler.cli.utils import ParameterType
from lsst.daf.butler.tests.cliOptionTest import OptionTestBase
from lsst.obs.base.cli.opt import instrument_parameter


class InstrumentTestCase(OptionTestBase):

    optionClass = instrument_parameter

    def verify(self, result, verifyArgs):
        self.assertEqual(result.exit_code, 0, f"output: {result.output} exception: {result.exception}")
        self.assertIn(verifyArgs, result.stdout)

    def verifyMissing(self, result, verifyArgs):
        self.assertNotEqual(result.exit_code, 0, f"output: {result.output} exception: {result.exception}")
        self.assertIn(verifyArgs, result.stdout)

    def test_argument(self):
        """test argument"""
        @click.command()
        @instrument_parameter(parameterType=ParameterType.ARGUMENT)
        def cli(instrument):
            if instrument is None:
                instrument = "None"
            print(instrument)

        self.run_test(cli, ["foo*"], self.verify, "foo*")
        self.run_test(cli, [], self.verify, "None")

    def test_argument_required(self):
        """test with argument required"""
        @click.command()
        @instrument_parameter(parameterType=ParameterType.ARGUMENT, required=True)
        def cli(instrument):
            print(instrument)

        self.run_test(cli, ["foo*"], self.verify, "foo*")
        self.run_test(cli, [], self.verifyMissing, 'Error: Missing argument "INSTRUMENT"')

    def test_option(self):
        """test option"""
        @click.command()
        @instrument_parameter()
        def cli(instrument):
            if instrument is None:
                instrument = "None"
            print(instrument)

        self.run_test(cli, ["--instrument", "foo*"], self.verify, "foo*")
        self.run_test(cli, [], self.verify, "None")

    def test_option_required(self):
        """test with argument required"""
        @click.command()
        @instrument_parameter(parameterType=ParameterType.ARGUMENT, required=True)
        def cli(instrument):
            print(instrument)

        self.run_test(cli, ["foo"], self.verify, "foo")
        self.run_test(cli, [], self.verifyMissing, 'Error: Missing argument "INSTRUMENT"')

    def test_argument_multiple(self):
        """test with multiple argument values"""
        @click.command()
        @instrument_parameter(parameterType=ParameterType.ARGUMENT, multiple=True)
        def cli(instrument):
            print(instrument)

        self.run_test(cli, ["foo*", "bar", "b?z"], self.verify, "('foo*', 'bar', 'b?z')")

    def test_option_multiple(self):
        """test with multiple option values"""
        @click.command()
        @instrument_parameter(multiple=True)
        def cli(instrument):
            print(instrument)

        self.run_test(cli, ["--instrument", "foo*", "--instrument", "bar", "--instrument", "b?z"],
                      self.verify, "('foo*', 'bar', 'b?z')")

    def test_help(self):
        self.help_test()
        self.custom_help_test()


if __name__ == "__main__":
    unittest.main()
