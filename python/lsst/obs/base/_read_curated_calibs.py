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

import datetime
import os
import re
from collections.abc import Mapping
from typing import TYPE_CHECKING, Any, NamedTuple, Protocol

from lsst.resources import ResourcePath

if TYPE_CHECKING:
    import lsst.afw.cameraGeom
    from lsst.resources import ResourcePathExpression


class _SearchPath(NamedTuple):
    root: ResourcePath
    children: tuple[str, ...]


class CuratedCalibration(Protocol):
    """Protocol that describes the methods needed by this class when dealing
    with curated calibration datasets.
    """

    @classmethod
    def readText(cls, path: str) -> CuratedCalibration: ...

    def getMetadata(self) -> Mapping: ...


def read_one_calib(
    path: _SearchPath,
    chip_id: int | None,
    filter_name: str | None,
    calib_class: type[CuratedCalibration],
) -> tuple[dict[datetime.datetime, CuratedCalibration], str]:
    """Read data for a particular path from the standard format at a
    particular root.

    Parameters
    ----------
    path : `_SearchPath`
        This tuple contains the top level of the data tree
        and then further optional subdirectories in subsequent
        indices.  See Notes below for more details.
    chip_id : `int` or None
        The identifier for the sensor in question.  To be used in
        validation.
    filter_name : `str` or None
        The identifier for the filter in question.  To be used in
        validation.
    calib_class : `typing.Any`
        The class to use to read the curated calibration text file. Must
        support the ``readText()`` method.

    Returns
    -------
    data_dict : `dict`
        A dictionary of calibration objects constructed from the
        appropriate factory class.  The key is the validity start time
        as a `datetime.datetime` object.
    calib_type : `str`
        The type of calibrations that have been read and are included
        in the ``data_dict``.

    Notes
    -----
    Curated calibrations are read from the appropriate ``obs_ _data``
    package, and are required to have a common directory structure to
    be identified and ingested properly.  The top-level directories
    are organized by the instrument's ``policyName``.  These names are
    generally all lower-case, but this is not universally true.

    Below the top-level instrument directory, subdirectories named
    after the curated calibration type contained within, with the
    dataset_type_name forced to lowercase.  For calibrations that
    depend on the detector (i.e., the defects), the next level of
    subdirectories should contain directories named with the detector
    name, again forced to lowercase.

    For filter dependent calibrations that do not depend on the
    detector (i.e., transmission_filter), the calibrations should be
    grouped into directories named with the physical filter name
    (again, all lowercase) below the dataset_type_name directory.
    Filter dependent calibrations that do depend on the detector
    (i.e., transmission_system), have physical filter named
    directories below the detector level directories.

    The files themselves are required to be named with the validity start
    date and the appropriate file extension. The date must be in an ISO format
    that can be parsed by `datetime.datetime.fromisoformat`. This can include
    variants commonly used data packages such as ``2025-04-30T12:23:00`` or the
    more compact ``20250430T123000``.
    """
    files: list[ResourcePath] = []
    extensions = (".ecsv", ".yaml", ".json")
    search_root = path.root
    for subdir in path.children:
        search_root = search_root.join(subdir, forceDirectory=True)

    file_filter = "(" + "|".join(re.escape(ext) for ext in extensions) + ")$"
    files = list(ResourcePath.findFileResources([search_root], file_filter=file_filter))

    # Convention is that data reside in location where the final two parts of
    # the directory root are <instrument>/<calib type>.
    # os.path.split() does not like a trailing "/" directory indicator.
    parts = os.path.split(path.root.path.removesuffix("/"))
    instrument = os.path.split(parts[0])[1]  # convention is that these reside at <instrument>/<calib_type>
    calib_type = parts[1]

    data_dict: dict[datetime.datetime, Any] = {}
    for f in files:
        date_str = f.updatedExtension("").basename()
        # Assume files are using some form of ISO date string.
        valid_start = datetime.datetime.fromisoformat(date_str)
        # For now readText does not know about URIs.
        with f.as_local() as local_uri:
            data_dict[valid_start] = calib_class.readText(local_uri.ospath)
        check_metadata(data_dict[valid_start], valid_start, instrument, chip_id, filter_name, f, calib_type)
    return data_dict, calib_type


def check_metadata(
    obj: Any,
    valid_start: datetime.datetime,
    instrument: str,
    chip_id: int | None,
    filter_name: str | None,
    filepath: ResourcePath,
    calib_type: str,
) -> None:
    """Check that the metadata is complete and self consistent.

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
    filter_name : `str`
        Identifier of the filter in question.
    filepath : `lsst.resources.ResourcePath`
        Path of the file read to construct the data.
    calib_type : `str`
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
    fcalib_type = md["OBSTYPE"]
    # These may optionally not exist.
    fchip_id = md.get("DETECTOR", None)
    ffilter_name = md.get("FILTER", None)

    if chip_id is not None:
        fchip_id = int(fchip_id)
    if filter_name is not None:
        ffilter_name = ffilter_name.lower()
        filter_name = filter_name.lower()

    if not (
        (finst.lower(), fchip_id, ffilter_name, fcalib_type.lower())
        == (instrument.lower(), chip_id, filter_name, calib_type.lower())
    ):
        raise ValueError(
            "Path and file metadata do not agree:\n"
            f"Path metadata: {instrument} {chip_id} {filter_name} {calib_type}\n"
            f"File metadata: {finst} {fchip_id} {ffilter_name} {fcalib_type}\n"
            f"File read from : {filepath}\n"
        )


def read_all(
    root: ResourcePathExpression,
    camera: lsst.afw.cameraGeom.Camera,
    calib_class: type[CuratedCalibration],
    required_dimensions: list[str],
    filters: set[str],
) -> tuple[dict[_SearchPath, dict[datetime.datetime, CuratedCalibration]], str]:
    """Read all data from the standard format at a particular root.

    Parameters
    ----------
    root : `lsst.resources.ResourcePathExpression`
        URI to the top level of the data tree.  This is expected to hold
        directories named after the sensor names.  They are expected to be
        lower case.
    camera : `lsst.afw.cameraGeom.Camera`
        The camera that goes with the data being read.
    calib_class : `typing.Any`
        The class to use to read the curated calibration text file. Must
        support the ``readText()`` and ``getMetadata()`` methods.
    required_dimensions : `list` [`str`]
        Dimensions required for the calibration.
    filters : `list` [`str`]
        List of the known filters for this camera.  Used to identify
        filter-dependent calibrations.

    Returns
    -------
    calibration_data : `dict`
        A dictionary of dictionaries of calibration objects
        constructed with the appropriate factory class. The first key
        is the sensor name in lower case, and the second is the
        validity start time as a `datetime.datetime` object.
    calib_type : `str`
        The type of calibrations that have been read and are included
        in the ``data_dict``.

    Notes
    -----
    Each leaf object in the constructed dictionary has metadata associated with
    it. The detector ID may be retrieved from the DETECTOR entry of that
    metadata.
    """
    calibration_data: dict[_SearchPath, dict[datetime.datetime, CuratedCalibration]] = {}

    root_uri = ResourcePath(root, forceDirectory=True, forceAbsolute=True)

    # Read all the subdirectories. These will be things like detectors or
    # physical filters. We assume that every sub directory contains calibration
    # data. We only walk the top level.
    _, dirs, _ = next(root_uri.walk())

    # If there are no sub directories this is a global calibration.
    if not dirs:
        dirs = [""]

    calib_types = set()
    # We assume the directories have been lowered.
    detector_map = {det.getName().lower(): det.getName() for det in camera}
    filter_map = {filterName.lower().replace(" ", "_"): filterName for filterName in filters}

    paths_to_search: list[_SearchPath] = []
    for dir_name in dirs:
        # Catch possible mistakes:
        if "detector" in required_dimensions:
            if dir_name not in detector_map:
                # Top level directories must be detectors if they're
                # required.
                detectors = list(detector_map)
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
            elif "physical_filter" in required_dimensions:
                # If the calibration depends on both detector and
                # physical_filter, the subdirs here should contain the
                # filter name.
                subdir = root_uri.join(dir_name, forceDirectory=True)
                _, subdirs, _ = next(subdir.walk())
                for subdir_name in subdirs:
                    if subdir_name not in filter_map:
                        raise RuntimeError(f"Filter {subdir_name} not known to supplied camera in {subdir}.")
                    else:
                        paths_to_search.append(_SearchPath(root_uri, (dir_name, subdir_name)))
            else:
                paths_to_search.append(_SearchPath(root_uri, (dir_name,)))
        elif "physical_filter" in required_dimensions:
            # If detector is not required, but physical_filter is,
            # then the top level should contain the filter
            # directories.
            if dir_name not in filter_map:
                raise RuntimeError(f"Filter {dir_name} not known to supplied camera in {root_uri}.")
            paths_to_search.append(_SearchPath(root_uri, (dir_name,)))
        else:
            # Neither detector nor physical_filter are required, so
            # the calibration is global, and will not be found in
            # subdirectories.
            paths_to_search.append(_SearchPath(root_uri, ()))

    if not paths_to_search:
        raise RuntimeError(f"Found no data files in {root_uri}")

    calib_type = "undefined"
    for path in paths_to_search:
        chip_id = None
        filter_name = None
        if "detector" in required_dimensions:
            chip_id = camera[detector_map[path.children[0]]].getId()
        if "physical_filter" in required_dimensions:
            filter_name = filter_map[path.children[-1]]

        calibration_data[path], calib_type = read_one_calib(path, chip_id, filter_name, calib_class)

        calib_types.add(calib_type)
        if len(calib_types) != 1:  # set.add(None) has length 1 so None is OK here.
            raise ValueError(f"Error mixing calib types: {calib_types}")

    no_data = all(v == {} for v in calibration_data.values())
    if no_data:
        raise RuntimeError("No data to ingest")

    return calibration_data, calib_type
