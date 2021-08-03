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
from lsst.obs.base import DefineVisitsTask, DefineVisitsConfig
from ..utils import getInstrument

log = logging.getLogger("defineVisits")


def defineVisits(repo, config_file, collections, instrument, processes=1):
    """Implements the command line interface `butler define-visits` subcommand,
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

    Notes
    -----
    Camera geometry is not currently found in registry but instead a default
    camera will be used for the relevant instrument.
    """
    if not collections:
        collections = None
    butler = Butler(repo, collections=collections, writeable=True)
    instr = getInstrument(instrument, butler.registry)
    config = DefineVisitsConfig()
    instr.applyConfigOverrides(DefineVisitsTask._DefaultName, config)

    if collections is None:
        # Default to the raw collection for this instrument
        collections = instr.makeDefaultRawIngestRunName()
        log.info("Defaulting to searching for raw exposures in collection %s", collections)

    if config_file is not None:
        config.load(config_file)
    task = DefineVisitsTask(config=config, butler=butler)

    # Assume the dataset type is "raw" -- this is required to allow this
    # query to filter out exposures not relevant to the specified collection.
    task.run(butler.registry.queryDataIds(["exposure"], dataId={"instrument": instr.getName()},
                                          collections=collections, datasets="raw"),
             collections=collections, processes=processes)
