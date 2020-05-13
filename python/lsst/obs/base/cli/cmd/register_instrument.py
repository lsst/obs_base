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

import click

from lsst.daf.butler.cli.opt import repo_argument
from lsst.daf.butler.cli.utils import cli_handle_exception
from ..opt import instrument_option
from ...script import registerInstrument


@click.command()
@repo_argument(required=True)
@instrument_option(required=True, help="The fully-qualified name of an Instrument subclass.")
def register_instrument(*args, **kwargs):
    """Add an instrument to the data repository.
    """
    cli_handle_exception(registerInstrument, *args, **kwargs)
