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
from lsst.pipe.base import Instrument

log = logging.getLogger(__name__)


def writeCuratedCalibrations(repo, instrument, collection, labels):
    """Add an instrument's curated calibrations to the data repository.

    Parameters
    ----------
    repo : `str`
        URI to the location to create the repo.
    instrument : `str`
        The name or the fully qualified class name of an instrument.
    collection : `str` or `None`
        The path to the collection that associates datasets with validity
        ranges.
        Can be `None` in which case the collection name will be determined
        automatically.
    labels : `Sequence` [ `str` ]
        Extra strings to include in the names of collections that datasets are
        inserted directly into, and if ``collection`` is `None`, the automatic
        calibration collection name as well.

    Raises
    ------
    RuntimeError
        Raised if the instrument can not be imported, instantiated, or obtained
        from the registry.
    TypeError
        Raised if the instrument is not a subclass of
        `lsst.obs.base.Instrument`.
    """
    butler = Butler(repo, writeable=True)
    instr = Instrument.from_string(instrument, butler.registry)
    if collection is None and not labels:
        labels = instr.get_curated_calibration_labels()
        if not labels:
            raise ValueError("At least one label or --collection must be provided.")
    instr.writeCuratedCalibrations(butler, collection=collection, labels=labels)
