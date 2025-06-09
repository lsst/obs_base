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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import annotations

__all__ = ()

import dataclasses
from collections.abc import Iterable, Mapping
from typing import TYPE_CHECKING, Any

import numpy as np
from tqdm.notebook import tqdm

from lsst.afw.cameraGeom import FIELD_ANGLE, PIXELS, Camera, Detector
from lsst.afw.geom import SkyWcs, makeSkyWcs
from lsst.daf.butler import Butler, DataCoordinate, DatasetNotFoundError, DimensionRecord
from lsst.daf.butler.queries import Query
from lsst.geom import Angle, Point2D, SpherePoint, degrees, SphereTransform
from lsst.sphgeom import ConvexPolygon, UnitVector3d

from ._instrument import Instrument, loadCamera

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure


@dataclasses.dataclass
class VisitDetectorGeometry:
    detector: Detector
    raw_wcs: SkyWcs
    raw_region: ConvexPolygon
    wcs_flip_x: bool
    fit_wcs: SkyWcs | None = None
    fit_region: ConvexPolygon | None = None

    @classmethod
    def build(
        cls, detector: Detector, instrument: Instrument, boresight: SpherePoint, orientation: Angle
    ) -> VisitDetectorGeometry:
        raw_formatter = instrument.getRawFormatter({"detector": detector.getId()})
        raw_wcs = raw_formatter.makeRawSkyWcsFromBoresight(boresight, orientation, detector)
        raw_region = ConvexPolygon([sp.getVector() for sp in raw_wcs.pixelToSky(detector.getCorners(PIXELS))])
        return cls(detector, raw_wcs, raw_region, wcs_flip_x=raw_formatter.wcsFlipX)


@dataclasses.dataclass
class VisitGeometry:
    visit_record: DimensionRecord
    exposure_records: list[DimensionRecord] = dataclasses.field(default_factory=list)
    detector_region_records: dict[int, DimensionRecord] = dataclasses.field(default_factory=dict)
    detectors: dict[int, VisitDetectorGeometry] = dataclasses.field(default_factory=dict)
    instrument: Instrument | None = None
    wcs_flip_x: bool | None = None

    @classmethod
    def from_query(cls, query: Query) -> list[VisitGeometry]:
        results: dict[DataCoordinate, VisitGeometry] = {}
        for data_id in query.data_ids(["visit", "exposure"]).with_dimension_records():
            if (visit_data := results.get(data_id.subset(["visit"]))) is None:
                visit_data = VisitGeometry(data_id.records["visit"])
                results[data_id.subset(["visit"])] = visit_data
            visit_data.exposure_records.append(data_id.records["exposure"])
        for vdr in query.dimension_records("visit_detector_region"):
            results[vdr.dataId.subset(["visit"])].detector_region_records[vdr.detector] = vdr
        return sorted(results.values(), key=lambda v: v.visit_record.dataId)

    @staticmethod
    def replace_butler_records(butler: Butler, visits: Iterable[VisitGeometry]) -> None:
        all_visit_records = []
        all_exposure_records = []
        all_visit_detector_region_records = []
        for visit in visits:
            all_visit_records.append(visit.visit_record)
            all_exposure_records.extend(visit.exposure_records)
            all_visit_detector_region_records.extend(visit.detector_region_records.values())
        butler.registry.insertDimensionData("visit", *all_visit_records, replace=True)
        butler.registry.insertDimensionData("exposure", *all_exposure_records, replace=True)
        butler.registry.insertDimensionData(
            "visit_detector_region", *all_visit_detector_region_records, replace=True
        )

    @property
    def data_id(self) -> DataCoordinate:
        return self.visit_record.dataId

    @property
    def boresight(self) -> SpherePoint:
        return SpherePoint(
            Angle(self.exposure_records[0].tracking_ra, degrees),
            Angle(self.exposure_records[0].tracking_dec, degrees),
        )

    @property
    def orientation(self) -> Angle:
        return Angle(self.exposure_records[0].sky_angle, degrees)

    def get_missing_detectors(self) -> set[int]:
        return {
            detector_id
            for detector_id in self.detector_region_records.keys()
            if detector_id not in self.detectors or self.detectors[detector_id].fit_wcs is None
        }

    def load_detector_geometry(
        self,
        butler: Butler,
        dataset_types: Iterable[str] = ("visit_summary", "preliminary_visit_summary"),
        collections: str | Iterable[str] | None = None,
        camera: Camera | None = None,
    ) -> None:
        self.instrument = Instrument.fromName(self.visit_record.instrument, butler.registry)
        if camera is None:
            camera, _ = loadCamera(butler, self.exposure_records[0].dataId, collections=collections)
        for detector_id in tqdm(self.detector_region_records.keys(), "Remaking raw WCS."):
            self.detectors[detector_id] = VisitDetectorGeometry.build(
                camera[detector_id],
                self.instrument,
                self.boresight,
                self.orientation,
            )
            self.wcs_flip_x = self.detectors[detector_id].wcs_flip_x
        for dataset_type in tqdm(dataset_types, "Loading visit summaries."):
            missing_detectors = self.get_missing_detectors()
            if not missing_detectors:
                break
            try:
                visit_summary = butler.get(dataset_type, self.visit_record.dataId, collections=collections)
            except DatasetNotFoundError:
                continue
            for detector_id in missing_detectors:
                if (detector_row := visit_summary.find(detector_id)) is not None:
                    if (wcs := detector_row.getWcs()) is not None:
                        detector_geometry = self.detectors[detector_id]
                        detector_geometry.fit_wcs = wcs
                        detector_geometry.fit_region = ConvexPolygon(
                            [
                                sp.getVector()
                                for sp in wcs.pixelToSky(detector_geometry.detector.getCorners(PIXELS))
                            ]
                        )

    def fit_pointing(self) -> tuple[SpherePoint, Angle]:
        fit_vectors = []
        raw_vectors = []
        y_axis_point: SpherePoint | None = None
        for detector_geometry in self.detectors.values():
            if detector_geometry.fit_region is not None:
                fit_vectors.extend(detector_geometry.fit_region.getVertices())
                raw_vectors.extend(detector_geometry.raw_region.getVertices())
            if y_axis_point is None:
                # Define a point at (0, 1deg) in the FIELD_ANGLE system,
                # according to the raw WCS, to let us fit the rotation angle.
                # We could do this with any detector and get the same answer.
                # We're going from field angle to pixels to sky even though the
                # raw WCS goes from pixels to field angle to sky, because that
                # minimizes how much this code knows about how the rotation
                # angle is defined.
                y_axis_point = detector_geometry.raw_wcs.pixelToSky(
                    detector_geometry.detector.transform(Point2D(0.0, np.pi / 180.0), FIELD_ANGLE, PIXELS)
                )
        # Fit (in the least-squares sense) for the 3-d rotation that best maps
        # the raw_vectors onto the fit_vectors.
        # See https://igl.ethz.ch/projects/ARAP/svd_rot.pdf for derivation.
        fit_matrix = np.array(fit_vectors)
        raw_matrix = np.array(raw_vectors)
        rotation = SphereTransform.fit_unit_vectors(raw_matrix, fit_matrix)
        new_boresight = rotation(self.boresight)
        new_y_axis_point = rotation(y_axis_point)
        new_orientation = Angle(90, degrees) - new_boresight.bearingTo(new_y_axis_point)
        if self.wcs_flip_x:
            # TODO: test this logic branch on an appropriate camera!  ComCam
            # doesn't exercise it.
            new_orientation = -new_orientation
        return new_boresight, new_orientation

    def fix_records(self) -> None:
        boresight, orientation = self.fit_pointing()
        self.exposure_records = [
            self._make_updated_dimension_record(
                record,
                tracking_ra=boresight.getRa().asDegrees(),
                tracking_dec=boresight.getDec().asDegrees(),
                sky_angle=orientation.asDegrees(),
            )
            for record in self.exposure_records
        ]
        all_vertices: list[UnitVector3d] = []
        for detector_id, detector_geometry in self.detectors.items():
            if detector_geometry.fit_region is not None:
                new_region = detector_geometry.fit_region
            else:
                fallback_detector_geometry = VisitDetectorGeometry.build(
                    detector_geometry.detector,
                    self.instrument,
                    boresight,
                    orientation,
                )
                new_region = fallback_detector_geometry.raw_region
            self.detector_region_records[detector_id] = self._make_updated_dimension_record(
                self.detector_region_records[detector_id], region=new_region
            )
            all_vertices.extend(new_region.getVertices())
        self.visit_record = self._make_updated_dimension_record(
            self.visit_record, region=ConvexPolygon.convexHull(all_vertices)
        )

    @staticmethod
    def _make_updated_dimension_record(original: DimensionRecord, **kwargs: Any) -> DimensionRecord:
        as_dict = original.toDict()
        as_dict.update(**kwargs)
        return type(original)(**as_dict)

    def plot(
        self,
        *,
        axes: Axes | None = None,
        figure: Figure | None = None,
        proj: SkyWcs | None = None,
        figure_kwargs: Mapping[str, Any] | None = None,
        fit_region_style: Mapping[str, Any] | None = None,
        raw_region_style: Mapping[str, Any] | None = None,
        fallback_region_style: Mapping[str, Any] | None = None,
        raw_axes_style: Mapping[str, Any] | None = None,
        fallback_axes_style: Mapping[str, Any] | None = None,
        plot_all_fallback: bool = False,
        detector_record_style: Mapping[str, Any] | None = None,
        visit_record_style: Mapping[str, Any] | None = None,
    ) -> tuple[Axes, SkyWcs]:
        if axes is None:
            if proj is None:
                proj = makeSkyWcs(Point2D(0.0, 0.0), self.boresight, np.identity(2))
            if figure is None:
                from lsst.utils.plotting import make_figure

                figure = make_figure(**(figure_kwargs or {}))
            axes = figure.add_subplot()
            axes.set_title(str(self.visit_record.dataId))
            axes.set_xlabel("local right ascension (deg)")
            axes.set_ylabel("local declination (deg)")
        elif proj is None:
            raise TypeError("proj must be provided if axes is provided.")
        if raw_axes_style is not None:
            self._plot_field_angle_axes(axes, proj, Angle(0.5, degrees), **raw_axes_style)
        if fallback_axes_style is not None or fallback_region_style is not None:
            new_boresight, new_orientation = self.fit_pointing()
            fallback_detector_geometry: VisitDetectorGeometry | None = None
            if fallback_region_style is not None:
                for original_detector_geometry in self.detectors.values():
                    if plot_all_fallback or original_detector_geometry.fit_wcs is None:
                        fallback_detector_geometry = VisitDetectorGeometry.build(
                            original_detector_geometry.detector,
                            self.instrument,
                            new_boresight,
                            new_orientation,
                        )
                        self._plot_polygon(
                            axes, proj, fallback_detector_geometry.raw_region, **fallback_region_style
                        )
            if fallback_region_style is not None:
                if fallback_detector_geometry is None:
                    fallback_detector_geometry = VisitDetectorGeometry.build(
                        next(iter(self.detectors.values())).detector,
                        self.instrument,
                        new_boresight,
                        new_orientation,
                    )
                self._plot_field_angle_axes(
                    axes,
                    proj,
                    Angle(0.5, degrees),
                    detector_geometry=fallback_detector_geometry,
                    **fallback_axes_style,
                )
        for detector_geometry in self.detectors.values():
            if detector_geometry.fit_region is not None and fit_region_style is not None:
                self._plot_polygon(axes, proj, detector_geometry.fit_region, **fit_region_style)
            if raw_region_style is not None:
                self._plot_polygon(axes, proj, detector_geometry.raw_region, **raw_region_style)
        if detector_record_style is not None:
            for record in self.detector_region_records.values():
                self._plot_polygon(axes, proj, record.region, **detector_record_style)
        if visit_record_style is not None:
            self._plot_polygon(axes, proj, self.visit_record.region, **visit_record_style)
        return axes, proj

    @staticmethod
    def _plot_polygon(
        axes: Axes,
        proj: SkyWcs,
        polygon: ConvexPolygon,
        **style: Any,
    ) -> None:
        vertices = np.array(proj.skyToPixel([SpherePoint(v) for v in polygon.getVertices()]))
        axes.fill(vertices[:, 0], vertices[:, 1], **style)

    def _plot_field_angle_axes(
        self,
        axes: Axes,
        proj: SkyWcs,
        size: Angle,
        detector_geometry: VisitDetectorGeometry | None = None,
        **style: Any,
    ) -> None:
        if detector_geometry is None:
            detector_geometry = next(iter(self.detectors.values()))
        origin_sky = detector_geometry.raw_wcs.pixelToSky(
            detector_geometry.detector.transform(Point2D(0.0, 0.0), FIELD_ANGLE, PIXELS)
        )
        y_axis_sky = detector_geometry.raw_wcs.pixelToSky(
            detector_geometry.detector.transform(Point2D(0.0, size.asRadians()), FIELD_ANGLE, PIXELS)
        )
        x_axis_sky = detector_geometry.raw_wcs.pixelToSky(
            detector_geometry.detector.transform(Point2D(size.asRadians(), 0.0), FIELD_ANGLE, PIXELS)
        )
        points = np.array(proj.skyToPixel([x_axis_sky, origin_sky, y_axis_sky]))
        axes.plot(points[:, 0], points[:, 1], **style)
        axes.text(points[0, 0], points[0, 1], "+X")
        axes.text(points[2, 0], points[2, 1], "+Y")
