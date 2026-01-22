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

from lsst.daf.butler import Butler
from lsst.pipe.base.configOverrides import ConfigOverrides
from lsst.utils import doImportType


def ingestRaws(
    repo,
    locations,
    regex,
    output_run,
    fail_fast=False,
    config=None,
    config_file=None,
    transfer="auto",
    processes=1,
    ingest_task="lsst.obs.base.RawIngestTask",
    track_file_attrs=True,
    update_records=False,
    skip_existing=False,
    search_indexes=True,
):
    """Ingest raw frames into the butler registry.

    Parameters
    ----------
    repo : `str`
        URI to the repository.
    locations : `list` [`str`]
        Files to ingest and directories to search for files that match
        ``regex`` to ingest.
    regex : `str`
        Regex string used to find files in directories listed in locations.
    output_run : `str`
        The path to the location, the run, where datasets should be put.
    fail_fast : `bool`
        If True, stop ingest as soon as any problem is encountered with any
        file. Otherwise problem files will be skipped and logged and a report
        issued at completion.
    config : `dict` [`str`, `str`] or `None`
        Key-value pairs to apply as overrides to the ingest config.
    config_file : `str` or `None`
        Path to a config file that contains overrides to the ingest config.
    transfer : `str` or None
        The external data transfer type, by default "auto".
    processes : `int`
        Number of processes to use for ingest.
    ingest_task : `str`
        The fully qualified class name of the ingest task to use by default
        lsst.obs.base.RawIngestTask.
    track_file_attrs : `bool`, optional
        Control whether file attributes such as the size or checksum should
        be tracked by the datastore. Whether this parameter is honored
        depends on the specific datastore implementation.
    update_records : `bool`, optional
        Control whether recalculated exposure definitions will be accepted or
        not.
    skip_existing : `bool`, optional
        Control whether raws that are already in the Butler repo will be
        skipped without error.
    search_indexes : `bool`, optional
        Control whether raw ingest will search for per-directory index files
        or not. Disabling this can improve performance if you know that there
        are no indexes.

    Raises
    ------
    Exception
        Raised if operations on configuration object fail.
    """
    TaskClass = doImportType(ingest_task)
    ingestConfig = TaskClass.ConfigClass()
    ingestConfig.transfer = transfer
    configOverrides = ConfigOverrides()
    if config_file is not None:
        configOverrides.addFileOverride(config_file)
    if config is not None:
        for name, value in config.items():
            configOverrides.addValueOverride(name, value)
    if fail_fast:
        configOverrides.addValueOverride("failFast", True)
    configOverrides.applyTo(ingestConfig)
    with Butler.from_config(repo, writeable=True) as butler:
        ingester = TaskClass(config=ingestConfig, butler=butler)
        ingester.run(
            locations,
            run=output_run,
            processes=processes,
            file_filter=regex,
            track_file_attrs=track_file_attrs,
            update_exposure_records=update_records,
            skip_existing_exposures=skip_existing,
            search_indexes=search_indexes,
        )
