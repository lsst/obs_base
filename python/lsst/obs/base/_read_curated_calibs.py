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

from __future__ import annotations

__all__ = ["CuratedCalibration", "read_all"]

import glob
import os
from collections.abc import Mapping
from typing import TYPE_CHECKING, Any, Protocol, Type

import dateutil.parser

if TYPE_CHECKING:
    import datetime

    import lsst.afw.cameraGeom


class CuratedCalibration(Protocol):
    """Protocol that describes the methods needed by this class when dealing
    with curated calibration datasets."""

    @classmethod
    def readText(cls, path: str) -> CuratedCalibration:
        ...

    def getMetadata(self) -> Mapping:
        ...


def read_one_calib(
    path: tuple,
    chip_id: int,
    filter_id: str,
    calib_class: Type[CuratedCalibration],
) -> tuple[dict[datetime.datetime, CuratedCalibration], str]:
    """Read data for a particular path from the standard format at a
    particular root.

    Parameters
    ----------
    path : `tuple` [`str`]
        This tuple contains the top level of the data tree at index=0,
        and then further optional subdirectories in subsequent
        indices.  Sensor directories are assumed to be at a higher
        level than any existin filter directories.  All entries are
        expected to be lower case.
    chip_id : `int`
        The identifier for the sensor in question.  To be used in
        validation.
    filter_id : `str`
        The identifier for the filter in question.  To be used in
        validation.
    calib_class : `Any`
        The class to use to read the curated calibration text file. Must
        support the ``readText()`` method.

    Returns
    -------
    `dict`
        A dictionary of objects constructed from the appropriate factory class.
        The key is the validity start time as a `datetime` object.
    """
    files = []
    extensions = (".ecsv", ".yaml", ".json")
    for ext in extensions:
        files.extend(glob.glob(os.path.join(*path, f"*{ext}")))

    parts = os.path.split(path[0])
    instrument = os.path.split(parts[0])[1]  # convention is that these reside at <instrument>/<data_name>
    data_name = parts[1]
    data_dict: dict[datetime.datetime, Any] = {}
    for f in files:
        date_str = os.path.splitext(os.path.basename(f))[0]
        valid_start = dateutil.parser.parse(date_str)
        data_dict[valid_start] = calib_class.readText(f)
        check_metadata(data_dict[valid_start], valid_start, instrument, chip_id, filter_id, f, data_name)
    return data_dict, data_name


def check_metadata(
    obj: Any, valid_start: datetime.datetime, instrument: str, chip_id: int, filter_id: str,
    filepath: str, data_name: str
) -> None:
    """Check that the metadata is complete and self consistent

    Parameters
    ----------
    obj : object of same type as the factory
        Object to retrieve metadata from in order to compare with
        metadata inferred from the path.
    valid_start : `datetime`
        Start of the validity range for data.
    instrument : `str`
        Name of the instrument in question.
    chip_id : `int`
        Identifier of the sensor in question.
    filter_id : `str`
        Identifier of the filter in question.
    filepath : `str`
        Path of the file read to construct the data.
    data_name : `str`
        Name of the type of data being read.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If the metadata from the path and the metadata encoded
        in the path do not match for any reason.
    """
    md = obj.getMetadata()
    # It is an error if these two do not exist.
    finst = md["INSTRUME"]
    fdata_name = md["OBSTYPE"]
    # These may optionally not exist.
    fchip_id = md.get("DETECTOR", None)
    ffilter_id = md.get("FILTER", None)

    if chip_id is not None:
        fchip_id = int(fchip_id)
    if filter_id is not None:
        ffilter_id = ffilter_id.lower()
        filter_id = filter_id.lower()

    if not (
        (finst.lower(), fchip_id, ffilter_id, fdata_name.lower())
            == (instrument.lower(), chip_id, filter_id, data_name.lower())
    ):
        raise ValueError(
            f"Path and file metadata do not agree:\n"
            f"Path metadata: {instrument} {chip_id} {filter_id} {data_name}\n"
            f"File metadata: {finst} {fchip_id} {ffilter_id} {fdata_name}\n"
            f"File read from : %s\n" % (filepath)
        )


def read_all(
    root: str,
    camera: lsst.afw.cameraGeom.Camera,
    calib_class: Type[CuratedCalibration],
    required_dimensions: list,
    filters: list,
) -> tuple[dict[str, dict[datetime.datetime, CuratedCalibration]], str]:
    """Read all data from the standard format at a particular root.

    Parameters
    ----------
    root : `str`
        Path to the top level of the data tree.  This is expected to hold
        directories named after the sensor names.  They are expected to be
        lower case.
    camera : `lsst.afw.cameraGeom.Camera`
        The camera that goes with the data being read.
    calib_class : `Any`
        The class to use to read the curated calibration text file. Must
        support the ``readText()`` and ``getMetadata()`` methods.
    required_dimensions : `list` [`str`]
        Dimensions required for the calibration.
    filters : `list` [`str`]
        List of the known filters for this camera.  Used to identify
        filter-dependent calibrations.

    Returns
    -------
    dict
        A dictionary of dictionaries of objects constructed with the
        appropriate factory class. The first key is the sensor name lowered,
        and the second is the validity start time as a `datetime` object.

    Notes
    -----
    Each leaf object in the constructed dictionary has metadata associated with
    it. The detector ID may be retrieved from the DETECTOR entry of that
    metadata.
    """
    calibration_data = {}

    root = os.path.normpath(root)
    dirs = os.listdir(root)  # assumes all directories contain data
    dirs = [d for d in dirs if os.path.isdir(os.path.join(root, d))]
    if not dirs:
        raise RuntimeError(f"No data found on path {root}")

    calib_types = set()
    detector_map = {
        det.getName().lower(): det.getName() for det in camera
    }  # we assume the directories have been lowered
    filter_map = {
        filterName.lower(): filterName for filterName in filters
    }

    paths_to_search = []
    for d in dirs:
        dir_name = os.path.basename(d)
        # Catch possible mistakes:
        if 'detector' in required_dimensions:
            if dir_name not in detector_map:
                # top level directories must be detectors if they're
                # required.
                detectors = [det for det in detector_map.keys()]
                max_detectors = 10
                note_str = "knows"
                if len(detectors) > max_detectors:
                    # report example subset
                    note_str = "examples"
                    detectors = detectors[:max_detectors]
                raise RuntimeError(
                    f"Detector {dir_name} not known to supplied camera "
                    f"{camera.getName()} ({note_str}: {','.join(detectors)})"
                )
            else:
                if 'physical_filter' in required_dimensions:
                    subdirs = os.listdir(os.path.join(root, chip_name))
                    subdirs = [d for d in subdirs if os.path.isdir(os.path.join(root, chip_name, d))]
                    for sd in subdirs:
                        subdir_name = os.path.basename(sd)
                        if subdir_name not in filter_map:
                            raise RuntimeError(f"Filter {subdir_name} not known to supplied camera.")
                        else:
                            paths_to_search.append((root, dir_name, subdir_name))
                else:
                    paths_to_search.append((root, dir_name))
        else:
            if 'physical_filter' in required_dimensions:
                if dir_name not in filter_map:
                    raise RuntimeError(f"Filter {dir_name} not known to supplied camera.")
                paths_to_search.append((root, dir_name))
            else:
                paths_to_search.append((root, ))

    for path in paths_to_search:
        chip_id = None
        filter_id = None
        if 'detector' in required_dimensions:
            chip_id = camera[detector_map[path[1]]].getId()
        if 'physical_filter' in required_dimensions:
            filter_id = filter_map[path[-1]]

        calibration_data[path], calib_type = read_one_calib(path, chip_id, filter_id, calib_class)

        calib_types.add(calib_type)
        if len(calib_types) != 1:  # set.add(None) has length 1 so None is OK here.
            raise ValueError(f"Error mixing calib types: {calib_types}")

    no_data = all([v == {} for v in calibration_data.values()])
    if no_data:
        raise RuntimeError("No data to ingest")

    return calibration_data, calib_type
