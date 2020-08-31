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


__all__ = ("RawIngestTask", "RawIngestConfig", "makeTransferChoiceField")

import os.path
from dataclasses import dataclass, InitVar
from typing import List, Iterator, Iterable, Type, Optional, Any
from collections import defaultdict
from multiprocessing import Pool

from astro_metadata_translator import ObservationInfo, fix_header, merge_headers
from lsst.afw.fits import readMetadata
from lsst.daf.butler import (
    Butler,
    CollectionType,
    DataCoordinate,
    DatasetRef,
    DatasetType,
    DimensionRecord,
    DimensionUniverse,
    FileDataset,
)
from lsst.pex.config import Config, ChoiceField
from lsst.pipe.base import Task

from ._instrument import Instrument, makeExposureRecordFromObsInfo
from ._fitsRawFormatterBase import FitsRawFormatterBase


@dataclass
class RawFileDatasetInfo:
    """Structure that holds information about a single dataset within a
    raw file.
    """

    dataId: DataCoordinate
    """Data ID for this file (`lsst.daf.butler.DataCoordinate`).
    """

    obsInfo: ObservationInfo
    """Standardized observation metadata extracted directly from the file
    headers (`astro_metadata_translator.ObservationInfo`).
    """


@dataclass
class RawFileData:
    """Structure that holds information about a single raw file, used during
    ingest.
    """

    datasets: List[RawFileDatasetInfo]
    """The information describing each dataset within this raw file.
    (`list` of `RawFileDatasetInfo`)
    """

    filename: str
    """Name of the file this information was extracted from (`str`).

    This is the path prior to ingest, not the path after ingest.
    """

    FormatterClass: Type[FitsRawFormatterBase]
    """Formatter class that should be used to ingest this file (`type`; as
    subclass of `FitsRawFormatterBase`).
    """

    instrumentClass: Type[Instrument]
    """The `Instrument` class associated with this file."""


@dataclass
class RawExposureData:
    """Structure that holds information about a complete raw exposure, used
    during ingest.
    """

    dataId: DataCoordinate
    """Data ID for this exposure (`lsst.daf.butler.DataCoordinate`).
    """

    files: List[RawFileData]
    """List of structures containing file-level information.
    """

    universe: InitVar[DimensionUniverse]
    """Set of all known dimensions.
    """

    record: Optional[DimensionRecord] = None
    """The exposure `DimensionRecord` that must be inserted into the
    `~lsst.daf.butler.Registry` prior to file-level ingest (`DimensionRecord`).
    """

    def __post_init__(self, universe: DimensionUniverse):
        # We don't care which file or dataset we read metadata from, because
        # we're assuming they'll all be the same; just use the first ones.
        self.record = makeExposureRecordFromObsInfo(self.files[0].datasets[0].obsInfo, universe)


def makeTransferChoiceField(doc="How to transfer files (None for no transfer).", default="auto"):
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
                 "auto": "choice will depend on datastore",
                 "link": "hard link falling back to symbolic link",
                 "hardlink": "hard link",
                 "symlink": "symbolic (soft) link",
                 "relsymlink": "relative symbolic link",
                 },
        optional=True,
        default=default
    )


class RawIngestConfig(Config):
    transfer = makeTransferChoiceField()


class RawIngestTask(Task):
    """Driver Task for ingesting raw data into Gen3 Butler repositories.

    Parameters
    ----------
    config : `RawIngestConfig`
        Configuration for the task.
    butler : `~lsst.daf.butler.Butler`
        Writeable butler instance, with ``butler.run`` set to the appropriate
        `~lsst.daf.butler.CollectionType.RUN` collection for these raw
        datasets.
    **kwargs
        Additional keyword arguments are forwarded to the `lsst.pipe.base.Task`
        constructor.

    Notes
    -----
    Each instance of `RawIngestTask` writes to the same Butler.  Each
    invocation of `RawIngestTask.run` ingests a list of files.
    """

    ConfigClass = RawIngestConfig

    _DefaultName = "ingest"

    def getDatasetType(self):
        """Return the DatasetType of the datasets ingested by this Task.
        """
        return DatasetType("raw", ("instrument", "detector", "exposure"), "Exposure",
                           universe=self.butler.registry.dimensions)

    def __init__(self, config: Optional[RawIngestConfig] = None, *, butler: Butler, **kwargs: Any):
        config.validate()  # Not a CmdlineTask nor PipelineTask, so have to validate the config here.
        super().__init__(config, **kwargs)
        self.butler = butler
        self.universe = self.butler.registry.dimensions
        self.datasetType = self.getDatasetType()

        # Import all the instrument classes so that we ensure that we
        # have all the relevant metadata translators loaded.
        Instrument.importAll(self.butler.registry)

    @classmethod
    # WARNING: this method hardcodes the parameters to pipe.base.Task.__init__.
    # Nobody seems to know a way to delegate them to Task code.
    def _makeTask(cls, config: RawIngestConfig, butler: Butler, name: str, parentTask: Task):
        """Construct a RawIngestTask using only positional arguments.

        Parameters
        ----------
        All parameters are as for `RawIngestTask`.
        """
        return cls(config=config, butler=butler, name=name, parentTask=parentTask)

    # Overrides Task.__reduce__
    def __reduce__(self):
        return (self._makeTask, (self.config, self.butler, self._name, self._parentTask))

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

        Notes
        -----
        Assumes that there is a single dataset associated with the given
        file.  Instruments using a single file to store multiple datasets
        must implement their own version of this method.
        """
        # Manually merge the primary and "first data" headers here because we
        # do not know in general if an input file has set INHERIT=T.
        phdu = readMetadata(filename, 0)
        header = merge_headers([phdu, readMetadata(filename)], mode="overwrite")
        fix_header(header)
        datasets = [self._calculate_dataset_info(header, filename)]

        # The data model currently assumes that whilst multiple datasets
        # can be associated with a single file, they must all share the
        # same formatter.
        instrument = Instrument.fromName(datasets[0].dataId["instrument"], self.butler.registry)
        FormatterClass = instrument.getRawFormatter(datasets[0].dataId)

        return RawFileData(datasets=datasets, filename=filename,
                           FormatterClass=FormatterClass,
                           instrumentClass=instrument)

    def _calculate_dataset_info(self, header, filename):
        """Calculate a RawFileDatasetInfo from the supplied information.

        Parameters
        ----------
        header : `Mapping`
            Header from the dataset.
        filename : `str`
            Filename to use for error messages.

        Returns
        -------
        dataset : `RawFileDatasetInfo`
            The dataId, and observation information associated with this
            dataset.
        """
        obsInfo = ObservationInfo(header)
        dataId = DataCoordinate.standardize(instrument=obsInfo.instrument,
                                            exposure=obsInfo.exposure_id,
                                            detector=obsInfo.detector_num,
                                            universe=self.universe)
        return RawFileDatasetInfo(obsInfo=obsInfo, dataId=dataId)

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
            exposure. All fields will be populated.  The
            `RawExposureData.dataId` attributes will be minimal (unexpanded)
            `DataCoordinate` instances.
        """
        exposureDimensions = self.universe["exposure"].graph
        byExposure = defaultdict(list)
        for f in files:
            # Assume that the first dataset is representative for the file
            byExposure[f.datasets[0].dataId.subset(exposureDimensions)].append(f)

        return [RawExposureData(dataId=dataId, files=exposureFiles, universe=self.universe)
                for dataId, exposureFiles in byExposure.items()]

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
            updated to data IDs for which `DataCoordinate.hasRecords` returns
            `True`.
        """
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
                self.butler.registry.dimensions["exposure"]: data.record,
            }
        )
        # Now we expand the per-file (exposure+detector) data IDs.  This time
        # we pass in the records we just retrieved from the exposure data ID
        # expansion.
        for file in data.files:
            for dataset in file.datasets:
                dataset.dataId = self.butler.registry.expandDataId(
                    dataset.dataId,
                    records=dict(data.dataId.records)
                )
        return data

    def prep(self, files, *, pool: Optional[Pool] = None, processes: int = 1) -> Iterator[RawExposureData]:
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

        # The next operation operates on RawExposureData instances (one at
        # a time) in-place and then returns the modified instance.  We call it
        # as a pass-through instead of relying on the arguments we pass in to
        # have been modified because in the parallel case those arguments are
        # going to be pickled and unpickled, and I'm not certain
        # multiprocessing is careful enough with that for output arguments to
        # work.

        # Expand the data IDs to include all dimension metadata; we need this
        # because we may need to generate path templates that rely on that
        # metadata.
        # This is the first step that involves actual database calls (but just
        # SELECTs), so if there's going to be a problem with connections vs.
        # multiple processes, or lock contention (in SQLite) slowing things
        # down, it'll happen here.
        return mapFunc(self.expandDataIds, exposureData)

    def ingestExposureDatasets(self, exposure: RawExposureData, *, run: Optional[str] = None
                               ) -> List[DatasetRef]:
        """Ingest all raw files in one exposure.

        Parameters
        ----------
        exposure : `RawExposureData`
            A structure containing information about the exposure to be
            ingested.  Must have `RawExposureData.records` populated and all
            data ID attributes expanded.
        run : `str`, optional
            Name of a RUN-type collection to write to, overriding
            ``self.butler.run``.

        Returns
        -------
        refs : `list` of `lsst.daf.butler.DatasetRef`
            Dataset references for ingested raws.
        """
        datasets = [FileDataset(path=os.path.abspath(file.filename),
                                refs=[DatasetRef(self.datasetType, d.dataId) for d in file.datasets],
                                formatter=file.FormatterClass)
                    for file in exposure.files]
        self.butler.ingest(*datasets, transfer=self.config.transfer, run=run)
        return [ref for dataset in datasets for ref in dataset.refs]

    def run(self, files, *, pool: Optional[Pool] = None, processes: int = 1, run: Optional[str] = None):
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
        run : `str`, optional
            Name of a RUN-type collection to write to, overriding
            the default derived from the instrument name.

        Returns
        -------
        refs : `list` of `lsst.daf.butler.DatasetRef`
            Dataset references for ingested raws.

        Notes
        -----
        This method inserts all datasets for an exposure within a transaction,
        guaranteeing that partial exposures are never ingested.  The exposure
        dimension record is inserted with `Registry.syncDimensionData` first
        (in its own transaction), which inserts only if a record with the same
        primary key does not already exist.  This allows different files within
        the same exposure to be incremented in different runs.
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
        runs = set()
        for exposure in exposureData:
            self.butler.registry.syncDimensionData("exposure", exposure.record)
            # Override default run if nothing specified explicitly
            if run is None:
                instrumentClass = exposure.files[0].instrumentClass
                this_run = instrumentClass.makeDefaultRawIngestRunName()
            else:
                this_run = run
            if this_run not in runs:
                self.butler.registry.registerCollection(this_run, type=CollectionType.RUN)
                runs.add(this_run)
            with self.butler.transaction():
                refs.extend(self.ingestExposureDatasets(exposure, run=this_run))
        return refs
