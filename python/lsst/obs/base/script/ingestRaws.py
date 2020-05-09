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

import os

from lsst.daf.butler import Butler
from lsst.pipe.base.configOverrides import ConfigOverrides
from lsst.utils import doImport


def ingestRaws(repo, output_run, config=None, config_file=None, directory=None, file=None, transfer="auto",
               ingest_task="lsst.obs.base.RawIngestTask"):
    """Ingests raw frames into the butler registry

    Parameters
    ----------
    repo : `str`
        URI to the repository.
    output_run : `str`
        The path to the location, the run, where datasets should be put.
    config : `dict` [`str`, `str`] or `None`
        Key-vaule pairs to apply as overrides to the ingest config.
    config_file : `str` or `None`
        Path to a config file that contains overrides to the ingest config.
    directory : `str` or `None`
        Path to the directory containing the raws to ingest.
    file : `str` or `None`
        Path to a file containing raws to ingest.
    transfer : `str` or None
        The external data transfer type, by default "auto".
    ingest_task : `str`
        The fully qualified class name of the ingest task to use by default
        lsst.obs.base.RawIngestTask.

    Raises
    ------
    `Exception` is raised if operations on configuration object fail.
    """
    butler = Butler(repo, run=output_run)
    TaskClass = doImport(ingest_task)
    ingestConfig = TaskClass.ConfigClass()
    ingestConfig.transfer = transfer
    configOverrides = ConfigOverrides()
    if config_file is not None:
        configOverrides.addFileOverride(config_file)
    if config is not None:
        for name, value in config:
            configOverrides.addValueOverride(name, value)
    configOverrides.applyTo(ingestConfig)
    ingester = TaskClass(config=ingestConfig, butler=butler)
    files = [file] if file is not None else []
    if directory is not None:
        files.extend(
            [os.path.join(directory, f) for f in os.listdir(directory) if f.lower().endswith("fits")])
    ingester.run(files)
