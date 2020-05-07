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
import os

from lsst.daf.butler.cli.opt import repo_argument, config_option, config_file_option, run_option
from lsst.daf.butler.cli.utils import cli_handle_exception
from lsst.daf.butler import Butler
from lsst.pipe.base.configOverrides import ConfigOverrides
from ..opt import instrument_option
from ... import RawIngestTask, RawIngestConfig
from ...utils import getInstrument

log = logging.getLogger(__name__)


@click.command()
@repo_argument(required=True)
@config_option()
@config_file_option()
@run_option(required=True)
@click.option("-d", "--dir", help="The path to the directory containing the raws to ingest.")
@click.option("-t", "--transfer", help="The external data transfer type.", default="auto")
def ingest_raws(repo, config, config_file, output_run, dir, transfer):
    """Ingests raw frames into the butler registry
    /f

    Parameters
    ----------
    repo : `str`
        URI to the repository.
    config : `dict` [`str`, `str`]
        Key-vaule pairs to apply as overrides to the ingest config.
    config_file : `str`
        Path to a config file that contains overrides to the ingest config.
    output_run : `str`
        The path to the location, the run, where datasets should be put.
    dir : `str`
        Path to the directory containing the raws to ingest.
    transfer : `str`
        The external data transfer type.
    """
    butler = Butler(repo, run=output_run)
    ingestConfig = RawIngestConfig()
    ingestConfig.transfer = transfer
    configOverrides = ConfigOverrides()
    if config_file is not None:
        configOverrides.addFileOverride(config_file)
    for name, value in config:
        configOverrides.addValueOverride(name, value)
    cli_handle_exception(configOverrides.applyTo, ingestConfig)
    ingester = RawIngestTask(config=ingestConfig, butler=butler)
    files = [os.path.join(dir, f) for f in os.listdir(dir) if f.endswith("fits") or f.endswith("FITS")]
    ingester.run(files)
