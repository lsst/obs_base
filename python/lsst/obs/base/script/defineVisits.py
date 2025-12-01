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
import logging

from lsst.daf.butler import Butler
from lsst.obs.base import DefineVisitsConfig, DefineVisitsTask
from lsst.pipe.base import Instrument

log = logging.getLogger("lsst.obs.base.defineVisits")


def defineVisits(
    repo,
    config_file,
    collections,
    instrument,
    where=None,
    update_records=False,
    incremental=False,
):
    """Implement the command line interface `butler define-visits` subcommand,
    should only be called by command line tools and unit test code that tests
    this function.

    Defines visits from exposures in the butler registry

    Parameters
    ----------
    repo : `str`
        URI to the location to create the repo.
    config_file : `str` or `None`
        Path to a config file that contains overrides to the ingest config.
    collections : `list` [`str`]
        An expression specifying the collections to be searched (in order) when
        reading datasets, and optionally dataset type restrictions on them.
        If empty it will be passed as `None` to Butler.
    instrument : `str`
        The name or fully-qualified class name of an instrument.
    where : `str`, optional
        Query clause to use when querying for dataIds. Can be used to limit
        the relevant exposures.
    update_records : `bool`, optional
        Control whether recalculated visit definitions will be accepted or
        not.
    incremental : `bool`, optional
        Declare that the visit definitions are being run in a situation
        where data from multi-snap visits are being ingested incrementally
        and so the visit definition could change as new data arrive.

    Notes
    -----
    Camera geometry is not currently found in registry but instead a default
    camera will be used for the relevant instrument.
    """
    if not collections:
        collections = None
    with Butler.from_config(repo, collections=collections, writeable=True) as butler:
        instr = Instrument.from_string(instrument, butler.registry)
        config = DefineVisitsConfig()
        instr.applyConfigOverrides(DefineVisitsTask._DefaultName, config)

        # If this is old schema but is using modern visit grouping algorithm,
        # (which is the default for new code) revert to one-to-one (which
        # was the old default).
        exposure_dimension = butler.dimensions["exposure"]
        modern = "one-to-one-and-by-counter"
        if "seq_end" not in exposure_dimension.metadata and config.groupExposures.name == modern:
            legacy = "one-to-one"
            log.warning(
                "Request to use %s grouping algorithm but registry schema is too old. Using %s instead.",
                modern,
                legacy,
            )
            config.groupExposures.name = legacy

        if not where:
            where = ""

        if config_file is not None:
            config.load(config_file)
        task = DefineVisitsTask(config=config, butler=butler)

        with butler.query() as query:
            query = query.join_dimensions(["exposure"]).where(where, instrument=instr.getName())
            data_ids = list(query.data_ids(["exposure"]).with_dimension_records())

        task.run(
            data_ids,
            collections=collections,
            update_records=update_records,
            incremental=incremental,
        )
