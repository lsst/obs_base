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

from lsst.daf.butler.cli.opt import repo_argument, run_option
from lsst.daf.butler import Butler
from ..opt import instrument_option
from ...utils import getInstrument

log = logging.getLogger(__name__)


@click.command()
@repo_argument(required=True)
@instrument_option(required=True)
@run_option(required=True)
@click.pass_context
def write_curated_calibrations(context, repo, instrument, output_run):
    """Add an instrument's curated calibrations to the data repository.
    """
    butler = Butler(repo, writeable=True, run=output_run)
    try:
        instr = getInstrument(instrument, butler.registry)
    except RuntimeError as err:
        log.critical(f"Exception getting instrument: {err}")
        raise click.BadParameter(f"Failed getting instrument {instrument} from repo {repo}")
    except TypeError as err:
        log.critical(f"{instrument} is not an Instrument subclass. {err}")
        raise click.BadParameter(f"{instrument} is not an Instrument subclass.")
    instr.writeCuratedCalibrations(butler)
