# This file is part of obs_base.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
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
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from __future__ import annotations

__all__ = (
    "CameraBuilderUpdater",
    "CameraLoader",
    "DetectorBuilderUpdater",
    "DetectorLoader",
)

from abc import ABC, abstractmethod
from typing import Any, Iterable, Optional, Union

from lsst.afw.cameraGeom import Camera, CameraBuilder, Detector, DetectorBuilder
from lsst.daf.butler import Butler, DataId, DatasetRef, DatasetType


class CameraBuilderUpdater(ABC):
    """An interface for classes that know how to modify a `CameraBuilder`.

    Notes
    -----
    It is expected that these objects will be saved to data repositories with
    data IDs that identify only the ``instrument`` dimension.  They are
    typically versioned in `~lsst.daf.butler.CollectionType.CALIBRATION`
    collections.

    This interface should only be implemented by classes that modify state
    associated with the full camera, such as the distortion model; similar
    classes whose instances correspond to a single detector should implement
    *only* `DetectorBuilderUpdater` instead.
    """

    @abstractmethod
    def update_camera_builder(self, builder: CameraBuilder) -> None:
        """Update the given camera builder with this object's state.

        Parameters
        ----------
        builder : `CameraBuilder`
            Camera builder to update.
        """
        raise NotImplementedError()


class DetectorBuilderUpdater(ABC):
    """An interface for classes that know how to modify a `DetectorBuilder`.

    Notes
    -----
    It is expected that these objects will be saved to data repositories with
    data IDs that identify only the ``instrument`` and ``detector`` dimensions.
    They are typically versioned in
    `~lsst.daf.butler.CollectionType.CALIBRATION` collections.

    This interface should only be implemented by classes that modify state
    associated with the full camera, such as the distortion model; similar
    classes whose instances correspond to a single detector should implement
    *only* `DetectorBuilderUpdater` instead.
    """

    @property
    @abstractmethod
    def detector_id(self) -> int:
        """The integer ID of the detector this object corresponds to (`int`).

        This is assumed to be both a valid index for a `Camera` object and a
        valid ``detector`` data ID value.
        """
        raise NotImplementedError()

    @abstractmethod
    def update_detector_builder(self, builder: DetectorBuilder) -> None:
        """Update the given detector builder with this object's state.

        Parameters
        ----------
        builder : `DetectorBuilder`
            Detector builder to update.
        """
        raise NotImplementedError()


class DetectorLoader:
    """Helper class for loading detectors from data repositories, including
    overrides of versioned content.

    Parameters
    ----------
    arg : `str`, `int`, `DatasetRef`, `Detector`, \
            `DetectorBuilder`
        If a `DatasetRef`, `Detector`, or `DetectorBuilder`,
        an object that can be used or loaded to provide the detector directly
        (``camera`` will be ignored).
        If an `int` or `str`, a detector identifier to use with ``camera``.
    butler : `Butler`
        Butler client to read from.  If not initialized with the desired
        collection search path, the ``collections`` argument to `update` must
        be used.
    camera : `str`, `DatasetType`, `DatasetRef`, `Camera`, \
            `CameraBuilder`, optional
        Used to obtain a `Camera` or `CameraBuilder` as a way to get a
        `DetectorBuilder`.  If a `Camera` is given, the detector is extracted
        and then a `Detector.rebuild` is called, instead of rebuilding the
        full camera.  Defaults to "camera" (as the nominal dataset type name).
    data_id : `dict` or `lsst.daf.butler.DataCoordinate`, optional
        Data ID values that identify the ``exposure`` and/or ``visit``
        dimension(s).  May include any extended keys supported by `Butler.get`,
        such as ``day_obs``.  Should not identify the ``detector`` dimension.
    **kwargs
        Additional keyword arguments are interpreted as data ID key-value pairs
        and forwarded directly to `Butler.get`.

    Notes
    -----
    Most code should use the higher-level `Instrument.load_detector` interface
    to this functionality instead of using this class directly.
    """

    def __init__(
        self,
        arg: Union[str, int, DatasetRef, Detector, DetectorBuilder],
        /,
        butler: Butler,
        camera: Union[str, DatasetType, DatasetRef, Camera, CameraBuilder] = "camera",
        data_id: Optional[DataId] = None,
        **kwargs: Any,
    ):
        if isinstance(arg, (int, str)):
            if isinstance(camera, CameraBuilder):
                self.builder = camera[arg]
            elif isinstance(camera, Camera):
                self.builder = camera[arg].rebuild()
            elif isinstance(camera, DatasetRef):
                self.builder = butler.getDirect(camera)[arg].rebuild()
            elif isinstance(arg, (str, DatasetType)):
                self.builder = butler.get(camera, data_id)[arg].rebuild()
            else:
                raise TypeError(f"Unrecognized argument for camera: {camera!r} ({type(arg)}).")
        elif isinstance(arg, DetectorBuilder):
            self.builder = arg
        elif isinstance(arg, Detector):
            self.builder = arg.rebuild()
        elif isinstance(arg, DatasetRef):
            self.builder = butler.getDirect(arg).rebuild()
        else:
            raise TypeError(f"Unrecognized first argument: {arg!r} ({type(arg)}).")
        self.butler = butler
        # We keep data ID and kwargs for data ID separate so we can pass them
        # as-is to Butler.get, which has special handling of non-dimension keys
        # that we currently can't call directly (DM-30439).
        self._data_id = data_id
        self._data_id_kwargs = kwargs

    def update(
        self,
        *,
        collections: Optional[Iterable[str]] = None,
        **kwargs: Union[bool, str, DatasetType, DatasetRef, DetectorBuilderUpdater],
    ) -> None:
        """Apply updates to the detector.

        Parameters
        ----------
        collections : `Iterable` [ `str` ], optional
            Collections to search, in order.  If not provided, the butler used
            to initialize this loader must already have the right collections.
        **kwargs
            Named updates to apply to the detector.  Keys are typically those
            in `Instrument.detector_calibrations`, but are only used here in
            error messages.  Values are one of the following:

            - `False`: do not update this dataset (same as not passing a key).
            - `str`, `DatasetType`: load from the butler with this custom
              dataset type.
            - `DatasetRef`: load with `Butler.getDirect` (must be a resolved
              reference).
            - `DetectorBuilderUpdater`: just apply this calibration directly.

        Notes
        -----
        This method only applies the updates given as keyword arguments, while
        the `Instrument.load_detector` method attempts to apply all updates
        potentially used by that instrument.
        """
        for component, arg in kwargs.items():
            if arg is False:
                continue
            elif isinstance(arg, DetectorBuilderUpdater):
                visitor = arg
            elif isinstance(arg, DatasetRef):
                visitor = self.butler.getDirect(arg)
            elif isinstance(arg, (str, DatasetType)):
                visitor = self.butler.get(
                    arg,
                    self._data_id,
                    collections=collections,
                    detector=self.builder.getId(),
                    **self._data_id_kwargs,
                )
            else:
                raise TypeError(f"Unrecognized value for {component}: {arg!r} ({type(arg)}).")
            visitor.update_detector_builder(self.builder)

    def finish(self) -> Detector:
        """Finish updating and return the updated `Detector`."""
        return self.builder.finish()


class CameraLoader:
    """Helper class for loading detectors from data repositories, including
    overrides of versioned content.

    Parameters
    ----------
    butler : `Butler`
        Butler client to read from.  If not initialized with the desired
        collection search path, the ``collections`` argument to `update` must
        be used.
    camera : `str`, `DatasetType`, `DatasetRef`, `Camera`, `CameraBuilder`, \
            optional
        A `CameraBuilder`, a `Camera` to rebuild into one, or a butler dataset
        type or reference that can be used to load a `Camera`.
        Defaults to "camera" (as the nominal dataset type name).
    data_id : `dict` or `lsst.daf.butler.DataCoordinate`, optional
        Data ID values that identify the ``exposure`` and/or ``visit``
        dimension(s).  May include any extended keys supported by `Butler.get`,
        such as ``day_obs``.
    **kwargs
        Additional keyword arguments are interpreted as data ID key-value pairs
        and forwarded directly to `Butler.get`.

    Notes
    -----
    Most code should use the higher-level `Instrument.load_camera` interface
    to this functionality instead of using this class directly.
    """

    def __init__(
        self,
        butler: Butler,
        camera: Union[str, DatasetType, DatasetRef, Camera, CameraBuilder] = "camera",
        data_id: Optional[DataId] = None,
        **kwargs: Any,
    ):
        if isinstance(camera, CameraBuilder):
            self.builder = camera
        elif isinstance(camera, Camera):
            self.builder = camera.rebuild()
        elif isinstance(camera, DatasetRef):
            self.builder = butler.getDirect(camera).rebuild()
        elif isinstance(camera, (str, DatasetType)):
            self.builder = butler.get(camera, data_id).rebuild()
        else:
            raise TypeError(f"Unrecognized camera argument: {camera!r} ({type(camera)}).")
        self.butler = butler
        # We keep data ID and kwargs for data ID separate so we can pass them
        # as-is to Butler.get, which has special handling of non-dimension keys
        # that we currently can't call directly (DM-30439).
        self._data_id = data_id
        self._data_id_kwargs = kwargs

    def update_camera(
        self,
        *,
        collections: Optional[Iterable[str]] = None,
        **kwargs: Union[
            bool,
            str,
            DatasetType,
            DatasetRef,
            CameraBuilderUpdater,
        ],
    ) -> None:
        """Apply updates to the camera as a whole.

        Parameters
        ----------
        collections : `Iterable` [ `str` ], optional
            Collections to search, in order.  If not provided, the butler used
            to initialize this loader must already have the right collections.
        **kwargs
            Keys are typically those in `Instrument.camera_calibrations`, but
            are only used here in error messages.  Values are one of the
            following:

            - `False`: do not update this dataset (same as not passing a key).
            - `str`, `DatasetType`: load from the butler with this custom
              dataset type.
            - `DatasetRef`: load with `Butler.getDirect` (must be a resolved
              reference).
            - `CameraBuilderUpdater`: just apply this calibration directly.

        Notes
        -----
        This method only applies the updates given as keyword arguments, while
        the `Instrument.load_camera` method attempts to apply all updates
        potentially used by that instrument.
        """
        for component, arg in kwargs.items():
            if arg is False:
                continue
            elif isinstance(arg, CameraBuilderUpdater):
                visitor = arg
            elif isinstance(arg, DatasetRef):
                visitor = self.butler.getDirect(arg)
            elif isinstance(arg, (str, DatasetType)):
                visitor = self.butler.get(arg, self._data_id, collections=collections, **self._data_id_kwargs)
            else:
                raise TypeError(f"Unrecognized value for {component}: {arg!r} ({type(arg)}).")
            visitor.update_camera_builder(self.builder)

    def update_detectors(
        self,
        *,
        collections: Optional[Iterable[str]] = None,
        **kwargs: Union[
            bool,
            str,
            DatasetType,
            Iterable[DatasetRef],
            Iterable[DetectorBuilderUpdater],
        ],
    ) -> None:
        """Apply updates to the camera's detectors.

        Parameters
        ----------
        collections : `Iterable` [ `str` ], optional
            Collections to search, in order.  If not provided, the butler used
            to initialize this loader must already have the right collections.
        **kwargs
            Keys are typically those in `Instrument.detector_calibrations`, but
            are only used here in error messages.  Values are one of the
            following:

            - `False`: do not update these datasets (same as not passing a
              key).
            - `str`, `DatasetType`: load from the butler with this custom
              dataset type.
            - `Iterable` of `DatasetRef`: load with `Butler.getDirect` (must be
              resolved references).
            - `Iterable` of `DetectorBuilderUpdater`: apply these calibrations
              directly.

        Notes
        -----
        This method only applies the updates given as keyword arguments, while
        the `Instrument.load_camera` method attempts to apply all updates
        potentially used by that instrument.
        """
        for component, arg in kwargs.items():
            if arg is False:
                continue
            if isinstance(arg, (str, DatasetType)):
                # It might be better to do queryDatasets here and then call
                # getDirect, but we need DM-30439 to allow non-dimension keys
                # in the data ID first.
                for detector_builder in self.builder:
                    visitor = self.butler.get(
                        arg,
                        self._data_id,
                        collections=collections,
                        detector=detector_builder.getId(),
                        **self._data_id_kwargs,
                    )
                    visitor.update_detector_builder(detector_builder)
            elif arg is True:
                raise TypeError("'True' is not a valid value for keyword argument {component!r}.")
            else:
                for item in arg:
                    if isinstance(item, DetectorBuilderUpdater):
                        visitor = item
                    elif isinstance(item, DatasetRef):
                        visitor = self.butler.getDirect(item)
                    else:
                        raise TypeError(f"Unrecognized item in {component}: {item!r} ({type(item)}).")
                    visitor.update_detector_builder(self.builder[visitor.detector_id])

    def finish(self) -> Camera:
        """Finish updating and return the updated `Camera`."""
        return self.builder.finish()
