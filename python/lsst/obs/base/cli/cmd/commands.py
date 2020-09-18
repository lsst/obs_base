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

from lsst.daf.butler.cli.opt import (repo_argument, config_option, config_file_option, locations_argument,
                                     file_override_options, regex_option, run_option, transfer_option)
from lsst.daf.butler.cli.utils import (cli_handle_exception, split_commas, typeStrAcceptsMultiple)
from ..opt import instrument_argument, instrument_option
from ... import script


# regular expression that can be used to find supported fits file extensions.
fits_re = r"\.fit[s]?\b"


@click.command(short_help="Convert a gen2 repo to gen3.")
@repo_argument(required=True,
               help="REPO is the URI or path to the gen3 repository. Will be created if it does not already "
               "exist")
@click.option("--gen2root", required=True,
              help="Root path of the gen 2 repo to be converted.")
@click.option("--skymap-name",
              help="Name of the new gen3 skymap (e.g. 'discrete/ci_hsc').")
@click.option("--skymap-config",
              help="Path to skymap config file defining the new gen3 skymap.")
@click.option("--calibs",
              help="Path to the gen 2 calibration repo. It can be absolute or relative to gen2root.")
@click.option("--reruns", multiple=True, callback=split_commas, metavar=typeStrAcceptsMultiple,
              help="List of gen 2 reruns to convert.")
@transfer_option(help="Mode to use to transfer files into the new repository.")
@config_file_option(help="Path to a `ConvertRepoConfig` override to be included after the Instrument config "
                    "overrides are applied.")
@file_override_options()
def convert(*args, **kwargs):
    """Convert a Butler gen 2 repository into a gen 3 repository."""
    cli_handle_exception(script.convert, *args, **kwargs)


@click.command(short_help="Define visits from exposures.")
@repo_argument(required=True)
@config_file_option(help="Path to a pex_config override to be included after the Instrument config overrides "
                         "are applied.")
@click.option("--collections",
              help="The collections to be searched (in order) when reading datasets.",
              multiple=True,
              callback=split_commas,
              metavar=typeStrAcceptsMultiple)
@instrument_option(required=True)
@file_override_options()
def define_visits(*args, **kwargs):
    """Define visits from exposures in the butler registry."""
    cli_handle_exception(script.defineVisits, *args, **kwargs)


@click.command(short_help="Ingest raw frames.")
@repo_argument(required=True)
@locations_argument(help="LOCATIONS specifies files to ingest and/or locations to search for files.",
                    required=True)
@regex_option(default=fits_re,
              help="Regex string used to find files in directories listed in LOCATIONS. "
                   "Searches for fits files by default.")
@config_option(metavar="TEXT=TEXT", multiple=True)
@config_file_option(type=click.Path(exists=True, writable=False, file_okay=True, dir_okay=False))
@run_option(required=False)
@transfer_option()
@click.option("--ingest-task", default="lsst.obs.base.RawIngestTask", help="The fully qualified class name "
              "of the ingest task to use.")
@file_override_options()
def ingest_raws(*args, **kwargs):
    """Ingest raw frames into from a directory into the butler registry"""
    cli_handle_exception(script.ingestRaws, *args, **kwargs)


@click.command(short_help="Add an instrument to the repository")
@repo_argument(required=True)
@instrument_argument(required=True, nargs=-1, help="The fully-qualified name of an Instrument subclass.")
def register_instrument(*args, **kwargs):
    """Add an instrument to the data repository.
    """
    cli_handle_exception(script.registerInstrument, *args, **kwargs)


@click.command(short_help="Add an instrument's curated calibrations.")
@repo_argument(required=True)
@instrument_option(required=True)
@run_option(required=False)
@file_override_options()
def write_curated_calibrations(*args, **kwargs):
    """Add an instrument's curated calibrations to the data repository.
    """
    cli_handle_exception(script.writeCuratedCalibrations, *args, **kwargs)
