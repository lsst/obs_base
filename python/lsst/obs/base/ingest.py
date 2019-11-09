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


__all__ = ("RawIngestTask", "RawIngestConfig", "makeTransferChoiceField")

import os.path
import itertools
from dataclasses import dataclass
from typing import List, Dict, Iterator, Iterable, Type, Optional, Any, Mapping
from collections import defaultdict
from multiprocessing import Pool

from astro_metadata_translator import ObservationInfo, fix_header, merge_headers
from lsst.utils import doImport
from lsst.afw.fits import readMetadata
from lsst.daf.butler import (
    Butler,
    DataCoordinate,
    DatasetRef,
    DatasetType,
    DimensionRecord,
    FileDataset,
)
from lsst.obs.base.instrument import makeExposureRecordFromObsInfo, makeVisitRecordFromObsInfo
from lsst.geom import Box2D
from lsst.pex.config import Config, Field, ChoiceField
from lsst.pipe.base import Task
from lsst.sphgeom import ConvexPolygon

from .fitsRawFormatterBase import FitsRawFormatterBase


@dataclass
class RawFileData:
    """Structure that holds information about a single raw file, used during
    ingest.
    """

    dataId: DataCoordinate
    """Data ID for this file (`lsst.daf.butler.DataCoordinate`).

    This may be a minimal `~lsst.daf.butler.DataCoordinate` base instance, or
    a complete `~lsst.daf.butler.ExpandedDataCoordinate`.
    """

    obsInfo: ObservationInfo
    """Standardized observation metadata extracted directly from the file
    headers (`astro_metadata_translator.ObservationInfo`).
    """

    region: ConvexPolygon
    """Region on the sky covered by this file, possibly with padding
    (`lsst.sphgeom.ConvexPolygon`).
    """

    filename: str
    """Name of the file this information was extracted from (`str`).

    This is the path prior to ingest, not the path after ingest.
    """

    FormatterClass: Type[FitsRawFormatterBase]
    """Formatter class that should be used to ingest this file and compute
    a spatial region for it (`type`; as subclass of `FitsRawFormatterBase`).
    """


@dataclass
class RawExposureData:
    """Structure that holds information about a complete raw exposure, used
    during ingest.
    """

    dataId: DataCoordinate
    """Data ID for this exposure (`lsst.daf.butler.DataCoordinate`).

    This may be a minimal `~lsst.daf.butler.DataCoordinate` base instance, or
    a complete `~lsst.daf.butler.ExpandedDataCoordinate`.
    """

    files: List[RawFileData]
    """List of structures containing file-level information.
    """

    records: Optional[Dict[str, List[DimensionRecord]]] = None
    """Dictionary containing `DimensionRecord` instances that must be inserted
    into the `~lsst.daf.butler.Registry` prior to file-level ingest (`dict`).

    Keys are the names of dimension elements ("exposure" and optionally "visit"
    and "visit_detector_region"), while values are lists of `DimensionRecord`.

    May be `None` during some ingest steps.
    """


def makeTransferChoiceField(doc="How to transfer files (None for no transfer).", default=None):
    """Create a Config field with options for how to transfer files between
    data repositories.

    The allowed options for the field are exactly those supported by
    `lsst.daf.butler.Datastore.ingest`.

    Parameters
    ----------
    doc : `str`
        Documentation for the configuration field.

    Returns
    -------
    field : `lsst.pex.config.ChoiceField`
        Configuration field.
    """
    return ChoiceField(
        doc=doc,
        dtype=str,
        allowed={"move": "move",
                 "copy": "copy",
                 "hardlink": "hard link",
                 "symlink": "symbolic (soft) link"},
        optional=True,
        default=default
    )


class RawIngestConfig(Config):
    transfer = makeTransferChoiceField()
    padRegionAmount = Field(
        dtype=int,
        default=0,
        doc="Pad an image with specified number of pixels before calculating region"
    )
    instrument = Field(
        doc=("Fully-qualified Python name of the `Instrument` subclass to "
             "associate with all raws."),
        dtype=str,
        optional=False,
        default=None,
    )


class RawIngestTask(Task):
    """Driver Task for ingesting raw data into Gen3 Butler repositories.

    This Task is intended to be runnable from the command-line, but it doesn't
    meet the other requirements of CmdLineTask or PipelineTask, and wouldn't
    gain much from being one.  It also wouldn't really be appropriate as a
    subtask of a CmdLineTask or PipelineTask; it's a Task essentially just to
    leverage the logging and configurability functionality that provides.

    Each instance of `RawIngestTask` writes to the same Butler.  Each
    invocation of `RawIngestTask.run` ingests a list of files.

    Parameters
    ----------
    config : `RawIngestConfig`
        Configuration for the task.
    butler : `~lsst.daf.butler.Butler`
        Butler instance.  Ingested Datasets will be created as part of
        ``butler.run`` and associated with its Collection.
    kwds
        Additional keyword arguments are forwarded to the `lsst.pipe.base.Task`
        constructor.

    Other keyword arguments are forwarded to the Task base class constructor.
    """

    ConfigClass = RawIngestConfig

    _DefaultName = "ingest"

    def getDatasetType(self):
        """Return the DatasetType of the Datasets ingested by this Task.
        """
        return DatasetType("raw", ("instrument", "detector", "exposure"), "Exposure",
                           universe=self.butler.registry.dimensions)

    def __init__(self, config: Optional[RawIngestConfig] = None, *, butler: Butler, **kwds: Any):
        super().__init__(config, **kwds)
        self.butler = butler
        self.universe = self.butler.registry.dimensions
        self.instrument = doImport(self.config.instrument)()
        # For now, we get a nominal Camera from the Instrument.
        # In the future, we may want to load one from a Butler calibration
        # collection that's appropriate for the observation timestamp of
        # the exposure.
        self.camera = self.instrument.getCamera()
        self.datasetType = self.getDatasetType()

    def extractMetadata(self, filename: str) -> RawFileData:
        """Extract and process metadata from a single raw file.

        Parameters
        ----------
        filename : `str`
            Path to the file.

        Returns
        -------
        data : `RawFileData`
            A structure containing the metadata extracted from the file,
            as well as the original filename.  All fields will be populated,
            but the `RawFileData.dataId` attribute will be a minimal
            (unexpanded) `DataCoordinate` instance.
        """
        phdu = readMetadata(filename, 0)
        header = merge_headers([phdu, readMetadata(filename)], mode="overwrite")
        fix_header(header)
        obsInfo = ObservationInfo(header)
        dataId = DataCoordinate.standardize(instrument=obsInfo.instrument,
                                            exposure=obsInfo.exposure_id,
                                            detector=obsInfo.detector_num,
                                            universe=self.universe)
        if obsInfo.instrument != self.instrument.getName():
            raise ValueError(f"Incorrect instrument (expected {self.instrument.getName()}, "
                             f"got {obsInfo.instrument}) for file {filename}.")
        FormatterClass = self.instrument.getRawFormatter(dataId)
        if obsInfo.visit_id is not None and obsInfo.tracking_radec is not None:
            formatter = FormatterClass.fromMetadata(metadata=header, obsInfo=obsInfo)
            visitInfo = formatter.makeVisitInfo()
            detector = self.camera[obsInfo.detector_num]
            wcs = formatter.makeWcs(visitInfo, detector)
            pixBox = Box2D(detector.getBBox())
            if self.config.padRegionAmount > 0:
                pixBox.grow(self.config.padRegionAmount)
            pixCorners = pixBox.getCorners()
            sphCorners = [wcs.pixelToSky(point).getVector() for point in pixCorners]
            region = ConvexPolygon(sphCorners)
        else:
            region = None
        return RawFileData(obsInfo=obsInfo, region=region, filename=filename,
                           FormatterClass=FormatterClass, dataId=dataId)

    def groupByExposure(self, files: Iterable[RawFileData]) -> List[RawExposureData]:
        """Group an iterable of `RawFileData` by exposure.

        Parameters
        ----------
        files : iterable of `RawFileData`
            File-level information to group.

        Returns
        -------
        exposures : `list` of `RawExposureData`
            A list of structures that group the file-level information by
            exposure.  The `RawExposureData.records` attributes of elements
            will be `None`, but all other fields will be populated.  The
            `RawExposureData.dataId` attributes will be minimal (unexpanded)
            `DataCoordinate` instances.
        """
        exposureDimensions = self.universe["exposure"].graph
        byExposure = defaultdict(list)
        for f in files:
            byExposure[f.dataId.subset(exposureDimensions)].append(f)

        return [RawExposureData(dataId=dataId, files=exposureFiles)
                for dataId, exposureFiles in byExposure.items()]

    def collectDimensionRecords(self, exposure: RawExposureData) -> RawExposureData:
        """Collect the `DimensionRecord` instances that must be inserted into
        the `~lsst.daf.butler.Registry` before an exposure's raw files may be.

        Parameters
        ----------
        exposure : `RawExposureData`
            A structure containing information about the exposure to be
            ingested.  Should be considered consumed upon return.

        Returns
        -------
        exposure : `RawExposureData`
            An updated version of the input structure, with
            `RawExposureData.records` populated.
        """
        firstFile = exposure.files[0]
        VisitDetectorRegionRecordClass = self.universe["visit_detector_region"].RecordClass
        exposure.records = {
            "exposure": [makeExposureRecordFromObsInfo(firstFile.obsInfo, self.universe)],
        }
        if firstFile.obsInfo.visit_id is not None:
            exposure.records["visit_detector_region"] = []
            visitVertices = []
            for file in exposure.files:
                if file.obsInfo.visit_id != firstFile.obsInfo.visit_id:
                    raise ValueError(f"Inconsistent visit/exposure relationship for "
                                     f"exposure {firstFile.obsInfo.exposure_id} between "
                                     f"{file.filename} and {firstFile.filename}: "
                                     f"{file.obsInfo.visit_id} != {firstFile.obsInfo.visit_id}.")
                if file.region is None:
                    self.log.warn("No region found for visit=%s, detector=%s.", file.obsInfo.visit_id,
                                  file.obsInfo.detector_num)
                    continue
                visitVertices.extend(file.region.getVertices())
                exposure.records["visit_detector_region"].append(
                    VisitDetectorRegionRecordClass.fromDict({
                        "instrument": file.obsInfo.instrument,
                        "visit": file.obsInfo.visit_id,
                        "detector": file.obsInfo.detector_num,
                        "region": file.region,
                    })
                )
            if visitVertices:
                visitRegion = ConvexPolygon(visitVertices)
            else:
                self.log.warn("No region found for visit=%s.", file.obsInfo.visit_id,
                              file.obsInfo.detector_num)
                visitRegion = None
            exposure.records["visit"] = [
                makeVisitRecordFromObsInfo(firstFile.obsInfo, self.universe, region=visitRegion)
            ]
        return exposure

    def expandDataIds(self, data: RawExposureData) -> RawExposureData:
        """Expand the data IDs associated with a raw exposure to include
        additional metadata records.

        Parameters
        ----------
        exposure : `RawExposureData`
            A structure containing information about the exposure to be
            ingested.  Must have `RawExposureData.records` populated. Should
            be considered consumed upon return.

        Returns
        -------
        exposure : `RawExposureData`
            An updated version of the input structure, with
            `RawExposureData.dataId` and nested `RawFileData.dataId` attributes
            containing `~lsst.daf.butler.ExpandedDataCoordinate` instances.
        """
        hasVisit = "visit" in data.records
        # We start by expanded the exposure-level data ID; we won't use that
        # directly in file ingest, but this lets us do some database lookups
        # once per exposure instead of once per file later.
        data.dataId = self.butler.registry.expandDataId(
            data.dataId,
            # We pass in the records we'll be inserting shortly so they aren't
            # looked up from the database.  We do expect instrument and filter
            # records to be retrieved from the database here (though the
            # Registry may cache them so there isn't a lookup every time).
            records={
                "exposure": data.records["exposure"][0],
                "visit": data.records["visit"][0] if hasVisit else None,
            }
        )
        # Now we expand the per-file (exposure+detector) data IDs.  This time
        # we pass in the records we just retrieved from the exposure data ID
        # expansion as well as the visit_detector_region record, if there is
        # one.
        vdrRecords = data.records["visit_detector_region"] if hasVisit else itertools.repeat(None)
        for file, vdrRecord in zip(data.files, vdrRecords):
            file.dataId = self.butler.registry.expandDataId(
                file.dataId,
                records=dict(data.dataId.records, visit_detector_region=vdrRecord)
            )
        return data

    def prep(self, files, pool: Optional[Pool] = None, processes: int = 1) -> Iterator[RawExposureData]:
        """Perform all ingest preprocessing steps that do not involve actually
        modifying the database.

        Parameters
        ----------
        files : iterable over `str` or path-like objects
            Paths to the files to be ingested.  Will be made absolute
            if they are not already.
        pool : `multiprocessing.Pool`, optional
            If not `None`, a process pool with which to parallelize some
            operations.
        processes : `int`, optional
            The number of processes to use.  Ignored if ``pool`` is not `None`.

        Yields
        ------
        exposure : `RawExposureData`
            Data structures containing dimension records, filenames, and data
            IDs to be ingested (one structure for each exposure).
        """
        if pool is None and processes > 1:
            pool = Pool(processes)
        mapFunc = map if pool is None else pool.imap_unordered

        # Extract metadata and build per-detector regions.
        fileData: Iterator[RawFileData] = mapFunc(self.extractMetadata, files)

        # Use that metadata to group files (and extracted metadata) by
        # exposure.  Never parallelized because it's intrinsically a gather
        # step.
        exposureData: List[RawExposureData] = self.groupByExposure(fileData)

        # The next few operations operate on RawExposureData instances (one at
        # a time) in-place and then return the modified instance.  We call them
        # as pass-throughs instead of relying on the arguments we pass in to
        # have been modified because in the parallel case those arguments are
        # going to be pickled and unpickled, and I'm not certain
        # multiprocessing is careful enough with that for output arguments to
        # work.  We use the same variable names to reflect the fact that we
        # consider the arguments to have been consumed/invalidated.

        # Extract DimensionRecords from the metadata that will need to be
        # inserted into the Registry before the raw datasets themselves are
        # ingested.
        exposureData: Iterator[RawExposureData] = mapFunc(self.collectDimensionRecords, exposureData)

        # Expand the data IDs to include all dimension metadata; we need this
        # because we may need to generate path templates that rely on that
        # metadata.
        # This is the first step that involves actual database calls (but just
        # SELECTs), so if there's going to be a problem with connections vs.
        # multiple processes, or lock contention (in SQLite) slowing things
        # down, it'll happen here.
        return mapFunc(self.expandDataIds, exposureData)

    def insertDimensionData(self, records: Mapping[str, List[DimensionRecord]]):
        """Insert dimension records for one or more exposures.

        Parameters
        ----------
        records : `dict` mapping `str` to `list`
            Dimension records to be inserted, organized as a mapping from
            dimension name to a list of records for that dimension.  This
            may be a single `RawExposureData.records` dict, or an aggregate
            for multiple exposures created by concatenating the value lists
            of those dictionaries.

        Returns
        -------
        refs : `list` of `lsst.daf.butler.DatasetRef`
            Dataset references for ingested raws.
        """
        # TODO: This currently assumes that either duplicate inserts of
        # visit records are ignored, or there is exactly one visit per
        # exposure.  I expect us to switch up the visit-exposure
        # relationship and hence rewrite some of this code before that
        # becomes a practical problem.
        # Iterate over dimensions explicitly to order for foreign key
        # relationships.
        for dimension in ("visit", "exposure", "visit_detector_region"):
            recordsForDimension = records.get(dimension)
            if recordsForDimension:
                # TODO: once Registry has options to ignore or replace
                # existing dimension records with the same primary keys
                # instead of aborting on conflicts, add configuration
                # options and logic to use them.
                self.butler.registry.insertDimensionData(dimension, *recordsForDimension)

    def ingestExposureDatasets(self, exposure: RawExposureData, butler: Optional[Butler] = None
                               ) -> List[DatasetRef]:
        """Ingest all raw files in one exposure.

        Parameters
        ----------
        exposure : `RawExposureData`
            A structure containing information about the exposure to be
            ingested.  Must have `RawExposureData.records` populated and all
            data ID attributes expanded.
        butler : `lsst.daf.butler.Butler`, optional
            Butler to use for ingest.  If not provided, ``self.butler`` will
            be used.

        Returns
        -------
        refs : `list` of `lsst.daf.butler.DatasetRef`
            Dataset references for ingested raws.
        """
        if butler is None:
            butler = self.butler
        datasets = [FileDataset(path=os.path.abspath(file.filename),
                                ref=DatasetRef(self.datasetType, file.dataId),
                                formatter=file.FormatterClass)
                    for file in exposure.files]
        butler.ingest(*datasets, transfer=self.config.transfer)
        return [dataset.ref for dataset in datasets]

    def run(self, files, pool: Optional[Pool] = None, processes: int = 1):
        """Ingest files into a Butler data repository.

        This creates any new exposure or visit Dimension entries needed to
        identify the ingested files, creates new Dataset entries in the
        Registry and finally ingests the files themselves into the Datastore.
        Any needed instrument, detector, and physical_filter Dimension entries
        must exist in the Registry before `run` is called.

        Parameters
        ----------
        files : iterable over `str` or path-like objects
            Paths to the files to be ingested.  Will be made absolute
            if they are not already.
        pool : `multiprocessing.Pool`, optional
            If not `None`, a process pool with which to parallelize some
            operations.
        processes : `int`, optional
            The number of processes to use.  Ignored if ``pool`` is not `None`.

        Returns
        -------
        refs : `list` of `lsst.daf.butler.DatasetRef`
            Dataset references for ingested raws.

        Notes
        -----
        This method inserts all records (dimensions and datasets) for an
        exposure within a transaction, guaranteeing that partial exposures
        are never ingested.
        """
        exposureData = self.prep(files, pool=pool, processes=processes)
        # Up to this point, we haven't modified the data repository at all.
        # Now we finally do that, with one transaction per exposure.  This is
        # not parallelized at present because the performance of this step is
        # limited by the database server.  That may or may not change in the
        # future once we increase our usage of bulk inserts and reduce our
        # usage of savepoints; we've tried to get everything but the database
        # operations done in advance to reduce the time spent inside
        # transactions.
        self.butler.registry.registerDatasetType(self.datasetType)
        refs = []
        for exposure in exposureData:
            with self.butler.transaction():
                self.insertDimensionData(exposure.records)
                refs.extend(self.ingestExposureDatasets(exposure))
        return refs
