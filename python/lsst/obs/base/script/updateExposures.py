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

import itertools
import logging

from lsst.daf.butler import Butler, DimensionRecord
from lsst.daf.butler.script.queryDatasets import QueryDatasets
from lsst.obs.base import Instrument
from lsst.obs.base.ingest import RawIngestTask
from lsst.pipe.base.configOverrides import ConfigOverrides
from lsst.utils import doImportType

_LOG = logging.getLogger(__name__)


def updateExposures(
    repo: str,
    instrument: str,
    where: str,
    collections: list[str] | None,
    fail_fast: bool = False,
    config: dict[str, str] | None = None,
    config_file: str | None = None,
    processes: int = 1,
    ingest_task: str = "lsst.obs.base.RawIngestTask",
) -> None:
    """Ingest raw frames into the butler registry.

    Parameters
    ----------
    repo : `str`
        URI to the repository.
    instrument : `str`
        Butler name of instrument. This is used to determine the default
        collection and the raw dataset type.
    where : `str`
        The query to use to discover the datasets. Does not need to include
        the instrument name. Must be specified to reduce the number of
        datasets being processed.
    collections : `list` [`str`] or `None`
        Override the default raw collection.
    fail_fast : `bool`, optional
        If True, stop ingest as soon as any problem is encountered with any
        file. Otherwise problem files will be skipped and logged and a report
        issued at completion.
    config : `dict` [`str`, `str`] or `None`, optional
        Key-value pairs to apply as overrides to the ingest config.
    config_file : `str` or `None`, optional
        Path to a config file that contains overrides to the ingest config.
    processes : `int`, optional
        Number of workers to use for ingest.
    ingest_task : `str`, optional
        The fully qualified class name of the ingest task to use by default
        lsst.obs.base.RawIngestTask.

    Raises
    ------
    Exception
        Raised if operations on configuration object fail.
    """
    if not where.strip():
        raise ValueError("A WHERE query string must be defined.")

    TaskClass = doImportType(ingest_task)
    assert issubclass(TaskClass, RawIngestTask)
    ingestConfig = TaskClass.ConfigClass()
    configOverrides = ConfigOverrides()
    if config_file is not None:
        configOverrides.addFileOverride(config_file)
    if config is not None:
        for name, value in config.items():
            configOverrides.addValueOverride(name, value)
    if fail_fast:
        configOverrides.addValueOverride("failFast", True)
    configOverrides.applyTo(ingestConfig)

    records_updated: list[DimensionRecord] = []

    def on_metadata_updates(record: DimensionRecord) -> None:
        records_updated.append(record)

    with Butler.from_config(repo, writeable=True) as butler:
        ingester = TaskClass(
            config=ingestConfig,  # type: ignore[arg-type]
            butler=butler,
            on_exposure_record=on_metadata_updates,
        )
        instr_ = Instrument.from_string(instrument, registry=butler.registry)
        if not collections:
            collections = [instr_.makeDefaultRawIngestRunName()]
        raw_datasetType = ingester.get_raw_datasetType(instr_)

        # Constrain the query by instrument. The QueryDatasets class does
        # not allow kwargs overrides.
        where = where + f" AND instrument = '{instr_.getName()}'"

        # Query for datasets.
        query = QueryDatasets(
            butler=butler,
            glob=[raw_datasetType.name],
            collections=collections,
            where=where,
            find_first=False,
            with_dimension_records=False,
            show_uri=False,
        )
        refs = set(itertools.chain(*query.getDatasets()))

        # The ingest code wants the raw file locations and we want the zips
        # without fragments.
        uris = butler.get_many_uris(refs)
        locations = {
            uri.replace(fragment="")
            for uri in itertools.chain.from_iterable(res.iter_all() for res in uris.values())
        }

        n_refs = len(refs)
        s_refs = "s" if n_refs != 1 else ""
        n_files = len(locations)
        s_files = "s" if n_files != 1 else ""
        _LOG.info(
            "Recalculating exposure records for %d dataset%s from %d file%s.",
            n_refs,
            s_refs,
            n_files,
            s_files,
        )

        ingester.run(
            locations,
            run=None,  # Not writing to a RUN collection.
            num_workers=processes,
            update_exposure_records=True,
            skip_ingest=True,
        )

        if records_updated:
            n_updated = len(records_updated)
            s_updated = "s" if n_updated != 1 else ""
            _LOG.info("Changed %d exposure record%s", n_updated, s_updated)
        else:
            _LOG.info("No exposure records were changed.")
