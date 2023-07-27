from __future__ import annotations

import dataclasses
import math
from collections.abc import Callable, Iterable

import tqdm
import astropy.time
from lsst.geom import Box2D
from lsst.sphgeom import ConvexPolygon
from lsst.afw.image import ExposureFitsReader
from lsst.daf.butler import Butler, DatasetRef, DatasetType, FileDataset, Timespan, DataCoordinate
from lsst.obs.base.formatters.fitsExposure import FitsExposureFormatter


class VisitSimulationImageReader:
    """A helper class for simulation ingest that reads metadata from on-disk
    `lsst.afw.image.ExposureF` FITS files that are simulations of {visit,
    detector} processed images.

    Parameters
    ----------
    filename : `str`
        Path to a file to be ingested.

    Notes
    -----
    This class can be subclassed in order to support files that lack some
    standard exposure components, by reimplementing methods to obtain the
    necessary information by parsing the filename or reading an external file.

    Once instance of this class is constructed for each file, but most methods
    are only called on the first file encountered for each visit.
    """

    def __init__(self, filename: str):
        self.filename = filename
        self.reader = ExposureFitsReader(filename)
        self.visit_info = self.reader.readVisitInfo()

    def get_visit_id(self) -> int:
        """Return the unique integer ID for this visit.

        This is shared by all per-detector files from a single visit.
        """
        return self.reader.readExposureId()

    def get_visit_name(self) -> str:
        """Return the alternate string identifer for this visit.

        This is shared by all per-detector files from a single visit.
        """
        return str(self.get_visit_id())

    def get_detector_id(self) -> int:
        """Return the unique integer ID for this detector.

        This must correspond to a detector ID already defined by the
        appropriate `lsst.obs.base.Instrument` class.
        """
        return self.reader.readDetector().getId()

    def get_physical_filter(self) -> str:
        """Return the name of the physical filter for this visit.

        This must correspond to a filter name already defined by the
        appropriate `lsst.obs.base.Instrument` class.  It is assumed to be the
        same for all detectors with the same ``visit_id``.  It may not be
        `None`.
        """
        return self.reader.readFilter().physicalLabel

    def get_spatial_region(self) -> ConvexPolygon:
        """Return the spatial region on the sky for this ``{visit, detector}``
        combination.

        This may not be `None`.
        """
        wcs = self.reader.readWcs()
        pixel_bbox = Box2D(self.reader.readBBox())
        sky_corners = wcs.pixelToSky(pixel_bbox.getCorners())
        return ConvexPolygon([corner.getVector() for corner in sky_corners])

    def get_timespan(self) -> Timespan | None:
        """Return the temporal bounds of this visit.

        This is assumed to be the same for all detectors with the same
        ``visit_id``.  It may be `None`.
        """
        if self.visit_info is None or not self.visit_info.getDate().isValid():
            return None
        mid = self.visit_info.getDate().toAstropy()
        if (total_s := self.get_exposure_time()) is None:
            return None
        half = astropy.time.TimeDelta(0.5 * total_s, format="sec", scale="tai")
        return Timespan(mid - half, mid + half)

    def get_exposure_time(self) -> float | None:
        """Return the duration of this visit in seconds.

        This is assumed to be the same for all detectors with the same
        ``visit_id``.  It may be `None`.
        """
        if self.visit_info is None:
            return None
        if math.isnan(result := self.visit_info.getExposureTime()):
            return None
        return result

    def get_science_program(self) -> str | None:
        """Return the science program of this visit (e.g. survey/simulation
        name).

        This is assumed to be the same for all detectors with the same
        ``visit_id``.  It may be `None`.
        """
        if self.visit_info is None:
            return None
        if not (result := self.visit_info.getScienceProgram()):
            return None
        return result

    def get_observation_reason(self) -> str | None:
        """Return the purpose of this visit (e.g. "science", "flat", etc.).

        This is assumed to be the same for all detectors with the same
        ``visit_id``.  It may be `None`, but defaults to "science" (the
        conventional value for standard, on-sky data) when not provided by the
        file's embedded visit information.
        """
        if self.visit_info is None:
            return "science"
        if not (result := self.visit_info.getObservationReason()):
            return "science"
        return result

    # There are a few metadata fields missing here, but it's not as obvious how
    # to get them from VisitInfo, and hopefully we just won't need them.


@dataclasses.dataclass
class VisitIngestionBundle:
    """A helper class for simulation ingest that holds already-extracted visit
    metadata and spatial regions.

    Notes
    -----
    One instance of this class is constructed for each visit, with one entry
    in its `detectors` attribute for each file associated with that visit.
    """

    visit_id: int
    """The unique integer ID for this visit."""

    name: str
    """The alternate string identifier for this visit."""

    physical_filter: str
    """The name of the physical filter for this visit."""

    timespan: Timespan | None = None
    """The temporal bounds of this visit."""

    exposure_time: float | None = None
    """The duration of this visit in seconds."""

    science_program: str | None = None
    """The science program of this visit (e.g. survey/simulation name)."""

    observation_reason: str | None = None
    """The purpose of this visit (e.g. "science", "flat", etc.)."""

    detectors: dict[int, tuple[ConvexPolygon, str]] = dataclasses.field(default_factory=dict)
    """Mapping of detectors in this visit.

    Keys are detector IDs, values are tuples of the per-{visit, detector}
    region and the filename.
    """


def ingest_visit_simulation_images(
    butler: Butler,
    instrument: str,
    filenames: Iterable[str],
    collection: str,
    transfer: str | None = "auto",
    dataset_type_name: str = "calexp",
    storage_class: str = "ExposureF",
    make_reader: Callable[[str], VisitSimulationImageReader] = VisitSimulationImageReader,
):
    """Ingest simulated ``{visit, detector}`` images into a butler, with
    simplified metadata.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
        Client for the butler reository. Must be read-write.
    instrument : `str`
        Short name for the instrument (the value returned by
        `lsst.pipe.base.Instrument.getName`).  This must correspond to an
        instrument that has already been registered with the data repository,
        defining all detectors and physical filters referenced by the simulated
        images.
    filenames : `~collections.abc.Iterable` [ `str` ]
        Iterable of filenames to process.  All files will be processed before
        and data repository writes occur to allow files from the same visit to
        be grouped together.  It may be necessary to chunk up very large lists
        of files into multiple calls (keeping visits together) to keep memory
        usage under control.
    collection : `str`
        Butler `~lsst.daf.butler.CollectionType.RUN` collection to ingest into.
        Will be created if it does not exist.
    transfer : `str`, optional
        How to transfer files to the data repository.  Common options are
        'copy', 'move', 'hardlink', and 'direct'; the latter leaves the file in
        its current location and just stores the absolute path, and guarantees
        that butler operations will never delete the file.
    dataset_type_name : `str`, optional
        Name of the dataset type to ingest files as.
    storage_class : `str`, optional
        Storage class (an object related to the image's in-memory Python type,
        defined in butler configuration) for this dataset type.
    make_reader : callable
        A callable that is given a filename and returns an instance of
        `VisitSimulationImageReader`.  This can be a subclass type if that
        subclass can be constructed from just a filename, or a more
        sophisticated factory (which could for example construct a subclass
        with additional metadata it loaded from a another location).

    Notes
    -----
    This function does not define 'exposure' metadata at all for simplicity; it
    assumes that the simulations represent images after the snaps have been
    combined into a single visit.  This means it cannot be used to ingest
    immediate post-ISR CCD images that mimick those produced by the real
    pipeline, as these are per-snap, not per-visit.

    This function also assumes that all files for a single visit will be
    ingested together, as it needs to group all files for each visit together
    to compute the spatial extent of that visit.
    """
    # Iterate over all files, read metadata, and group by visit.
    bundles: dict[int, VisitIngestionBundle] = {}
    for filename in tqdm.tqdm(filenames, desc="Reading file metadata."):
        reader = make_reader(filename)
        visit_id = reader.get_visit_id()
        if (bundle := bundles.get(visit_id)) is None:
            bundle = VisitIngestionBundle(
                visit_id=visit_id,
                name=reader.get_visit_name(),
                physical_filter=reader.get_physical_filter(),
                timespan=reader.get_timespan(),
                exposure_time=reader.get_exposure_time(),
                science_program=reader.get_science_program(),
                observation_reason=reader.get_observation_reason(),
            )
            bundles[visit_id] = bundle
        bundle.detectors[reader.get_detector_id()] = (reader.get_spatial_region(), filename)
    dataset_type = DatasetType(
        dataset_type_name, {"visit", "detector"}, storageClass=storage_class, universe=butler.dimensions
    )
    butler.registry.registerDatasetType(dataset_type)
    butler.registry.registerRun(collection)
    # Iterate over per-visit bundles to generate full-focal-plane visit
    # regions, insert dimension records, and actually ingest files.
    for visit_id, bundle in tqdm.tqdm(bundles.items(), desc="Ingesting visits."):
        vertices = []
        for detector_polygon, _ in bundle.detectors.values():
            vertices.extend(detector_polygon.getVertices())
        visit_polygon = ConvexPolygon.convexHull(vertices)
        visit_record = butler.dimensions["visit"].RecordClass(
            instrument=instrument,
            id=visit_id,
            name=bundle.name,
            physical_filter=bundle.physical_filter,
            timespan=bundle.timespan,
            exposure_time=bundle.exposure_time,
            observation_reason=bundle.observation_reason,
            science_program=bundle.science_program,
            region=visit_polygon,
        )
        # Actually insert the visit dimension record, but check to see if it
        # already exists.
        if not butler.registry.syncDimensionData("visit", visit_record):
            # If necessary we can make this case work by computing the union
            # of current region and the new one, but that's a lot more work;
            # let's hope we don't need that.
            raise NotImplementedError(
                f"Visit {visit_id} was already present; partial-visit ingests not supported."
            )
        # Insert dimension records that hold the spatial region for each
        # detector (everything else is visit-wide).
        region_records = [
            butler.dimensions["visit_detector_region"].RecordClass(
                instrument=instrument,
                detector=detector_id,
                visit=visit_id,
                region=detector_polygon,
            )
            for detector_id, (detector_polygon, filename) in bundle.detectors.items()
        ]
        butler.registry.insertDimensionData("visit_detector_region", *region_records)
        # Actually ingest the image files.
        file_datasets = [
            FileDataset(
                refs=[
                    DatasetRef(
                        dataset_type,
                        DataCoordinate.standardize(
                            visit=visit_id,
                            detector=detector_id,
                            instrument=instrument,
                            universe=butler.dimensions,
                        ),
                        run=collection,
                    )
                ],
                path=filename,
                formatter=FitsExposureFormatter,
            )
            for detector_id, (_, filename) in bundle.detectors.items()
        ]
        butler.ingest(*file_datasets, transfer=transfer)


def main() -> None:
    """Run simplified ingest from the command-line."""
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("root", help="Butler repository root or config file.")
    parser.add_argument("instrument", help="Short name of the instrument, as it appears in data IDs.")
    parser.add_argument("collection", help="RUN collection to ingest into.")
    parser.add_argument("filenames", nargs="+", help="Files to ingest.")
    namespace = parser.parse_args()
    butler = Butler(namespace.root, writeable=True)
    ingest_visit_simulation_images(butler, namespace.instrument, namespace.filenames, namespace.collection)


if __name__ == "__main__":
    main()
