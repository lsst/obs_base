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

from collections.abc import Sequence
from typing import Any

import click

from lsst.daf.butler import Butler
from lsst.daf.butler.cli.opt import (
    collections_argument,
    config_file_option,
    config_option,
    locations_argument,
    options_file_option,
    processes_option,
    regex_option,
    repo_argument,
    run_option,
    transfer_option,
    where_option,
)
from lsst.daf.butler.cli.utils import ButlerCommand, split_commas, typeStrAcceptsMultiple
from lsst.pipe.base.cli.opt import instrument_argument

from ... import script
from ..opt import failfast_option, labels_argument

# regular expression that can be used to find supported fits file extensions.
fits_re = r"\.fit[s]?\b"


@click.command(short_help="Define visits from exposures.", cls=ButlerCommand)
@repo_argument(required=True)
@instrument_argument(required=True)
@config_file_option(
    help="Path to a pex_config override to be included after the Instrument config overrides are applied."
)
@click.option(
    "--collections",
    help="The collections to be searched (in order) for camera geometry.",
    multiple=True,
    callback=split_commas,
    metavar=typeStrAcceptsMultiple,
)
@where_option()
@click.option(
    "--update-records/--no-update-records",
    default=False,
    help="Use this option to force updates to the visit definition record. "
    "Should only be used if you know that there has been a change to the "
    "exposure records, such as a change to the metadata translator.",
)
@click.option(
    "--incremental/--no-incremental",
    default=False,
    help="Use this option to force updates to the visit definition record "
    "when multi-snap visits are being ingested incrementally and so you "
    "might encounter partial visits.  Implies --update-records.",
)
@options_file_option()
def define_visits(*args: Any, **kwargs: Any) -> None:
    """Define visits from exposures in the butler registry.

    The calibration collection containing the camera geometry can not
    be specified.
    """
    script.defineVisits(*args, **kwargs)


@click.command(short_help="Ingest raw frames.", cls=ButlerCommand)
@repo_argument(required=True)
@locations_argument(
    help="LOCATIONS specifies files to ingest and/or locations to search for files.", required=True
)
@regex_option(
    default=fits_re,
    help="Regex string used to find files in directories listed in LOCATIONS. "
    "Searches for fits files by default.",
)
@config_option(metavar="TEXT=TEXT", multiple=True)
@config_file_option(type=click.Path(exists=True, writable=False, file_okay=True, dir_okay=False))
@run_option(required=False)
@transfer_option()
@processes_option(help="Number of workers to use.")
@click.option(
    "--ingest-task",
    default="lsst.obs.base.RawIngestTask",
    help="The fully qualified class name of the ingest task to use.",
)
@click.option(
    "--track-file-attrs/--no-track-file-attrs",
    default=True,
    help="Indicate to the datastore whether file attributes such as file size"
    " or checksum should be tracked or not. Whether this parameter is honored"
    " depends on the specific datastore implentation.",
)
@failfast_option()
@click.option(
    "--update-records/--no-update-records",
    default=False,
    help="Use this option to force updates to the exposure records. "
    "Should only be used if you know that there has been a change to the "
    "exposure records, such as a change to the metadata translator.",
)
@click.option(
    "--skip-existing/--no-skip-existing",
    default=True,
    help="Use this option to control whether a raw file already existing "
    "in the Butler repository is an error or not.",
)
@click.option(
    "--search-indexes/--no-search-indexes",
    default=True,
    help="By default ingest searches for index files in each directory being "
    "searched. Turning this off can improve performance if you know for a fact "
    "you have no indexes and you are using an object store.",
)
@options_file_option()
def ingest_raws(*args: Any, **kwargs: Any) -> None:
    """Ingest raw frames into from a directory into the butler registry."""
    script.ingestRaws(*args, **kwargs)


@click.command(short_help="Add an instrument's curated calibrations.", cls=ButlerCommand)
@repo_argument(required=True)
@instrument_argument(required=True)
@labels_argument()
@click.option(
    "--collection",
    required=False,
    help="Name of the calibration collection that associates datasets with validity ranges.",
)
@click.option(
    "--label",
    "labels",
    multiple=True,
    help=(
        "Extra strings to include (with automatic delimiters) in all RUN collection names, "
        "as well as the calibration collection name if it is not provided via --collection. "
        "May be provided as a positional argument instead of or in addition to this option, "
        "as long as at least one label is provided.  Positional-argument labels appear before "
        "those provided by this option."
    ),
)
@click.option(
    "--prefix",
    required=False,
    help=(
        "Prefix for the collection name.  Default is the instrument name. "
        "This is ignored if --collection is passed."
    ),
)
@options_file_option()
def write_curated_calibrations(
    *,
    repo: str,
    instrument: str,
    collection: str,
    labels: list[str],
    labels_arg: list[str],
    prefix: str,
) -> None:
    """Add an instrument's curated calibrations to the data repository."""
    script.writeCuratedCalibrations(
        repo=repo, instrument=instrument, collection=collection, labels=labels_arg + labels, prefix=prefix
    )


@click.command()
@repo_argument(required=True)
@click.argument("instrument")
@collections_argument(help="Collections to search for visit geometry datasets.")
@click.option("--dataset-type", default="visit_geometry", help="Name of the visit geometry dataset type.")
@click.option("--batch-size", default=11, help="Number of visits to process in each transaction.")
@where_option()
def update_dimension_regions(
    *,
    repo: str,
    instrument: str,
    collections: Sequence[str],
    dataset_type: str,
    batch_size: int,
    where: str,
) -> None:
    """Update dimension record regions from visit geometry datasets."""
    from ...visit_geometry import VisitGeometry

    butler = Butler.from_config(repo, writeable=True, collections=collections)
    VisitGeometry.update_dimension_records(
        butler, instrument, where, dataset_type=dataset_type, batch_size=batch_size
    )
