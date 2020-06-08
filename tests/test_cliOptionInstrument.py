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

from lsst.daf.butler.tests.cliOptionTest import OptionTestBase
from lsst.obs.base.cli.opt import instrument_option


class InstrumentOptionTestCase(OptionTestBase):

    optionClass = instrument_option

    def test_set(self):

        @click.command()
        @instrument_option()
        def cli(instrument):
            click.echo(instrument, nl=False)

        def verify(result, verifyArgs):
            self.assertEqual(result.exit_code, 0, f"output: {result.output} exception: {result.exception}")
            self.assertEqual(result.stdout, "myInstr")

        self.run_test(cli, ["--instrument", "myInstr"], verify, None)

    def test_required(self):

        @click.command()
        @instrument_option(required=True)
        def cli(instrument):
            click.echo(instrument, nl=False)

        def verify(result, verifyArgs):
            self.assertNotEqual(result.exit_code, 0, f"output: {result.output} exception: {result.exception}")
            self.assertIn('Missing option "-i" / "--instrument"', result.output)

        self.run_test(cli, [], verify, None)

    def test_notRequired(self):

        @click.command()
        @instrument_option()
        def cli(instrument):
            click.echo(instrument, nl=False)

        def verify(result, verifyArgs):
            self.assertEqual(result.exit_code, 0, f"output: {result.output} exception: {result.exception}")
            self.assertEqual(result.stdout, "")

        self.run_test(cli, [], verify, None)

    def test_help(self):
        self.custom_help_test()


if __name__ == "__main__":
    unittest.main()
