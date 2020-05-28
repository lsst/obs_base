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

"""Convert a gen2 butler repo to gen 3. See
`lsst.obs.base.ConvertRepoConfig` for most of the config options.
"""

import lsst.daf.butler
import lsst.log
import lsst.utils

from ..gen2to3 import ConvertRepoTask, ConvertRepoSkyMapConfig, Rerun
from ..utils import getInstrument


def convert(repo, gen2root, instrument, skymap_name, skymap_config, calibs, reruns, config_file, transfer):
    """Implements the command line interface `butler convert` subcommand,
    should only be called by command line tools and unit test code that tests
    this function.

    Convert the gen 2 Butler repo at `gen2root` into a gen 3 repo
    living at `repo`.

    Parameters
    ----------
    repo : `str`
        URI to the gen 3 repository.
    gen2root : `str`
        URI to the gen 2 repository.
    instrument : `str`
        The fully-qualified name of the Instrument subclass being converted.
    skymap_name : `str` or None
        Name of the skymap to be converted in the repo.
    skymap_config : `str` or None
        Path to the `lsst.skymap.BaseSkyMapConfig` of the gen2 skymap to be
        converted.
    calibs : `str` or None
        Path to the gen2 calibration repository to be converted.
        If a relative path, it is assumed to be relative to `gen2root`.
    reruns : `list` [`str`] or None
        List of reruns to convert. They will be placed in the
        ``shared/INSTRUMENT/RERUN`` collection.
    config_file : `str` or None
        Path to `lsst.obs.base.ConvertRepoConfig` configuration to load
        after all default/instrument configurations.
    transfer : `str` or None
        Mode to use when transferring data into the gen3 repository.
    """
    # Allow a gen3 butler to be reused
    try:
        butlerConfig = lsst.daf.butler.Butler.makeRepo(repo)
    except FileExistsError:
        # Use the existing butler configuration
        butlerConfig = repo

    butler = lsst.daf.butler.Butler(butlerConfig)
    instr = getInstrument(instrument, butler.registry)

    convertRepoConfig = ConvertRepoTask.ConfigClass()
    instr.applyConfigOverrides(ConvertRepoTask._DefaultName, convertRepoConfig)
    convertRepoConfig.instrument = instrument
    convertRepoConfig.raws.transfer = transfer
    if skymap_name is not None:
        convertRepoConfig.skyMaps[skymap_name] = ConvertRepoSkyMapConfig()
        convertRepoConfig.skyMaps[skymap_name].load(skymap_config)
        convertRepoConfig.rootSkyMapName = skymap_name
    if config_file is not None:
        convertRepoConfig.load(config_file)

    if reruns is None:
        rerunsArg = []
    else:
        rerunsArg = [Rerun(rerun, runName=f"shared/{instr.getName()}/{rerun}",
                           chainName=f"shared/{instr.getName()}", parents=[]) for rerun in reruns]

    # create a new butler instance for running the convert repo task
    butler = lsst.daf.butler.Butler(butlerConfig, run=f"raw/{instr.getName()}")
    convertRepoTask = ConvertRepoTask(config=convertRepoConfig, butler3=butler)
    convertRepoTask.run(
        root=gen2root,
        reruns=rerunsArg,
        calibs=None if calibs is None else {calibs: f"calib/{instr.getName()}"}
    )
