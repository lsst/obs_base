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
import logging

from lsst.daf.butler.cli.opt.repo import repo_option
from lsst.daf.butler import Butler
from ..opt import instrument_option
from ...utils import getInstrument

log = logging.getLogger(__name__)


@click.command()
@instrument_option(required=True, helpMsg="The fully-qualified name of an Instrument subclass.")
@repo_option(required=True)
@click.pass_context
def register_instrument(context, repo, instrument):
    """Add an instrument to the data repository.
    """
    butler = Butler(repo, writeable=True)
    try:
        instr = getInstrument(instrument, butler.registry)
    except RuntimeError as err:
        log.critical("Failed getting instrument %s with exception %s", instrument, err)
        raise click.ClickException("Could not import instrument.")
    except TypeError:
        raise click.ClickException(f"{instrument} is not a subclass of obs.base.Instrument")
    instr.register(butler.registry)
