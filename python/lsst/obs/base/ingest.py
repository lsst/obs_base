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
from typing import Callable, List, Iterator, Iterable, Tuple, Type, Optional, Any, Union
from collections import defaultdict
from multiprocessing import Pool

from astro_metadata_translator import ObservationInfo, merge_headers
from astro_metadata_translator.indexing import read_sidecar, read_index
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
    Formatter,
)
from lsst.pex.config import Config, ChoiceField, Field
from lsst.pipe.base import Task, timeMethod

from ._instrument import Instrument, makeExposureRecordFromObsInfo
from ._fitsRawFormatterBase import FitsRawFormatterBase


def _do_nothing(*args, **kwargs) -> None:
    """A function that accepts anything and does nothing, for use as a default
    in callback arguments.
    """
    pass


def _log_msg_counter(noun: Union[int, Iterable]) -> Tuple[int, str]:
    """Function to count the iterable and return the count and plural modifier.

    Parameters
    ----------
    noun : Iterable or `int`
        Thing to count. If given an integer it is assumed to be the count
        to use to calculate modifier.

    Returns
    -------
    num : `int`
        Number of items found in ``noun``.
    modifier : `str`
        Character to add to the end of a string referring to these items
        to indicate whether it was a single item or not. Returns empty
        string if there is one item or "s" otherwise.
    """
    if isinstance(noun, int):
        num = noun
    else:
        num = len(noun)
    return num, "" if num == 1 else "s"


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

    instrumentClass: Optional[Type[Instrument]]
    """The `Instrument` class associated with this file. Can be `None`
    if ``datasets`` is an empty list."""


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
                 "direct": "use URI to ingested file directly in datastore",
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
    failFast = Field(
        dtype=bool,
        default=False,
        doc="If True, stop ingest as soon as any problem is encountered with any file. "
            "Otherwise problems files will be skipped and logged and a report issued at completion.",
    )


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
    on_success : `Callable`, optional
        A callback invoked when all of the raws associated with an exposure
        are ingested.  Will be passed a list of `FileDataset` objects, each
        containing one or more resolved `DatasetRef` objects.  If this callback
        raises it will interrupt the entire ingest process, even if
        `RawIngestConfig.failFast` is `False`.
    on_metadata_failure : `Callable`, optional
        A callback invoked when a failure occurs trying to translate the
        metadata for a file.  Will be passed the filename and the exception, in
        that order, as positional arguments.  Guaranteed to be called in an
        ``except`` block, allowing the callback to re-raise or replace (with
        ``raise ... from``) to override the task's usual error handling (before
        `RawIngestConfig.failFast` logic occurs).
    on_ingest_failure : `Callable`, optional
        A callback invoked when dimension record or dataset insertion into the
        database fails for an exposure.  Will be passed a `RawExposureData`
        instance and the exception, in that order, as positional arguments.
        Guaranteed to be called in an ``except`` block, allowing the callback
        to re-raise or replace (with ``raise ... from``) to override the task's
        usual error handling (before `RawIngestConfig.failFast` logic occurs).
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

    def __init__(self, config: Optional[RawIngestConfig] = None, *, butler: Butler,
                 on_success: Callable[[List[FileDataset]], Any] = _do_nothing,
                 on_metadata_failure: Callable[[str, Exception], Any] = _do_nothing,
                 on_ingest_failure: Callable[[RawExposureData, Exception], Any] = _do_nothing,
                 **kwargs: Any):
        config.validate()  # Not a CmdlineTask nor PipelineTask, so have to validate the config here.
        super().__init__(config, **kwargs)
        self.butler = butler
        self.universe = self.butler.registry.dimensions
        self.datasetType = self.getDatasetType()
        self._on_success = on_success
        self._on_metadata_failure = on_metadata_failure
        self._on_ingest_failure = on_ingest_failure

        # Import all the instrument classes so that we ensure that we
        # have all the relevant metadata translators loaded.
        Instrument.importAll(self.butler.registry)

    def _reduce_kwargs(self):
        # Add extra parameters to pickle
        return dict(**super()._reduce_kwargs(), butler=self.butler, on_success=self._on_success,
                    on_metadata_failure=self._on_metadata_failure, on_ingest_failure=self._on_ingest_failure)

    def _determine_instrument_formatter(self, dataId, filename):
        """Determine the instrument and formatter class.

        Parameters
        ----------
        dataId : `lsst.daf.butler.DataCoordinate`
            The dataId associated with this dataset.
        filename : `str`
            Filename used for error reporting.

        Returns
        -------
        instrument : `Instrument` or `None`
            Instance of the `Instrument` associated with this dataset. `None`
            indicates that the instrument could not be determined.
        formatterClass : `type`
            Class to be used as the formatter for this dataset.
        """
        # The data model currently assumes that whilst multiple datasets
        # can be associated with a single file, they must all share the
        # same formatter.
        try:
            instrument = Instrument.fromName(dataId["instrument"], self.butler.registry)
        except LookupError as e:
            self._on_metadata_failure(filename, e)
            self.log.warning("Instrument %s for file %s not known to registry",
                             dataId["instrument"], filename)
            if self.config.failFast:
                raise RuntimeError(f"Instrument {dataId['instrument']} for"
                                   f" file {filename} not known to registry") from e
            FormatterClass = Formatter
            # Indicate that we could not work out the instrument.
            instrument = None
        else:
            FormatterClass = instrument.getRawFormatter(dataId)
        return instrument, FormatterClass

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
            (unexpanded) `~lsst.daf.butler.DataCoordinate` instance. The
            ``instrumentClsas`` field will be `None` if there is a problem
            with metadata extraction.

        Notes
        -----
        Assumes that there is a single dataset associated with the given
        file.  Instruments using a single file to store multiple datasets
        must implement their own version of this method.

        By default the method will catch all exceptions unless the ``failFast``
        configuration item is `True`.  If an error is encountered the
        `_on_metadata_failure()` method will be called. If no exceptions
        result and an error was encountered the returned object will have
        a null-instrument class and no datasets.

        This method supports sidecar JSON files which can be used to
        extract metadata without having to read the data file itself.
        The sidecar file is always used if found.
        """

        sidecar_fail_msg = ""
        try:
            root, ext = os.path.splitext(filename)
            sidecar_file = root + ".json"
            if os.path.exists(sidecar_file):
                header = read_sidecar(sidecar_file)
                sidecar_fail_msg = " (via sidecar)"
            else:
                # Read the metadata from the data file itself
                # Manually merge the primary and "first data" headers here
                # because we do not know in general if an input file has
                # set INHERIT=T.
                phdu = readMetadata(filename, 0)
                header = merge_headers([phdu, readMetadata(filename)], mode="overwrite")
            datasets = [self._calculate_dataset_info(header, filename)]
        except Exception as e:
            self.log.debug("Problem extracting metadata from %s%s: %s", filename, sidecar_fail_msg, e)
            # Indicate to the caller that we failed to read
            datasets = []
            formatterClass = Formatter
            instrument = None
            self._on_metadata_failure(filename, e)
            if self.config.failFast:
                raise RuntimeError("Problem extracting metadata for file "
                                   f"{filename}{sidecar_fail_msg}") from e
        else:
            self.log.debug("Extracted metadata for file %s%s", filename, sidecar_fail_msg)
            # The data model currently assumes that whilst multiple datasets
            # can be associated with a single file, they must all share the
            # same formatter.
            instrument, formatterClass = self._determine_instrument_formatter(datasets[0].dataId, filename)
            if instrument is None:
                datasets = []

        return RawFileData(datasets=datasets, filename=filename,
                           FormatterClass=formatterClass,
                           instrumentClass=instrument)

    def _calculate_dataset_info(self, header, filename):
        """Calculate a RawFileDatasetInfo from the supplied information.

        Parameters
        ----------
        header : `Mapping` or `astro_metadata_translator.ObservationInfo`
            Header from the dataset or previously-translated content.
        filename : `str`
            Filename to use for error messages.

        Returns
        -------
        dataset : `RawFileDatasetInfo`
            The dataId, and observation information associated with this
            dataset.
        """
        # To ensure we aren't slowed down for no reason, explicitly
        # list here the properties we need for the schema
        # Use a dict with values a boolean where True indicates
        # that it is required that we calculate this property.
        ingest_subset = {
            "altaz_begin": False,
            "boresight_rotation_coord": False,
            "boresight_rotation_angle": False,
            "dark_time": False,
            "datetime_begin": True,
            "datetime_end": True,
            "detector_num": True,
            "exposure_group": False,
            "exposure_id": True,
            "exposure_time": True,
            "instrument": True,
            "tracking_radec": False,
            "object": False,
            "observation_counter": False,
            "observation_id": True,
            "observation_reason": False,
            "observation_type": True,
            "observing_day": False,
            "physical_filter": True,
            "science_program": False,
            "visit_id": False,
        }

        if isinstance(header, ObservationInfo):
            obsInfo = header
            missing = []
            # Need to check the required properties are present
            for property, required in ingest_subset.items():
                if not required:
                    continue
                value = getattr(obsInfo, property)
                if value is None:
                    missing.append(property)
            if missing:
                raise ValueError(f"Requested required properties are missing from file {filename}:"
                                 f" {missing} (via JSON)")

        else:
            obsInfo = ObservationInfo(header, pedantic=False, filename=filename,
                                      required={k for k in ingest_subset if ingest_subset[k]},
                                      subset=set(ingest_subset))

        dataId = DataCoordinate.standardize(instrument=obsInfo.instrument,
                                            exposure=obsInfo.exposure_id,
                                            detector=obsInfo.detector_num,
                                            universe=self.universe)
        return RawFileDatasetInfo(obsInfo=obsInfo, dataId=dataId)

    def locateAndReadIndexFiles(self, files):
        """Given a list of files, look for index files and read them.

        Index files can either be explicitly in the list of files to
        ingest, or else located in the same directory as a file to ingest.
        Index entries are always used if present.

        Parameters
        ----------
        files : iterable over `str` or path-like objects
            Paths to the files to be ingested.  Will be made absolute
            if they are not already.

        Returns
        -------
        index : `dict` [`str`, Any]
            Merged contents of all relevant index files found. These can
            be explicitly specified index files or ones found in the
            directory alongside a data file to be ingested.
        updated_files : `set[str]`
            Updated list of the input files with entries removed that were
            found listed in an index file.
        bad_index_files: `set[str]`
            Files that looked like index files but failed to read properly.
        """
        # Create set of input files to allow easy removal of entries and
        # convert to absolute path for easy comparison with index content.
        # Do not convert to real paths since we have to assume that index
        # files are in this location and not the location which it links to.
        files_set = {os.path.abspath(f) for f in files}

        # Index files must be named this
        index_root_file = "_index.json"

        # Group the files by directory
        files_by_directory = defaultdict(set)

        for path in files_set:
            directory, file_in_dir = os.path.split(path)
            files_by_directory[directory].add(file_in_dir)

        # All the metadata read from index files with keys of full path
        index_entries = {}

        # Index files we failed to read
        bad_index_files = set()

        # Any good index files that were found and used
        good_index_files = set()

        # Look for index files in those directories
        for directory, files_in_directory in files_by_directory.items():
            possible_index_file = os.path.join(directory, index_root_file)
            if os.path.exists(possible_index_file):
                # If we are explicitly requesting an index file the
                # messages should be different
                index_msg = "inferred"
                is_implied = True
                if index_root_file in files_in_directory:
                    index_msg = "explicit"
                    is_implied = False

                # Try to read the index file and catch and report any
                # problems.
                try:
                    index = read_index(possible_index_file, force_dict=True)
                except Exception as e:
                    # For now only trigger the callback if the index file
                    # was asked for explicitly.
                    if not is_implied:
                        self._on_metadata_failure(possible_index_file, e)
                    if self.config.failFast:
                        raise RuntimeError(f"Problem reading index file from {index_msg} "
                                           f"location {possible_index_file}") from e
                    bad_index_files.add(possible_index_file)
                    continue

                self.log.debug("Extracted index metadata from %s file %s", index_msg, possible_index_file)
                good_index_files.add(possible_index_file)

                # Go through the index adding entries for files.
                # If we have non-index files in this directory marked for
                # ingest we should only get index information for those.
                # If the index file was explicit we use all entries.
                if is_implied:
                    files_to_ingest = files_in_directory
                else:
                    files_to_ingest = set(index)

                # Copy relevant metadata into a single dict for all index
                # entries.
                for file_in_dir in files_to_ingest:
                    # Skip an explicitly specified index file.
                    if file_in_dir == index_root_file:
                        continue
                    if file_in_dir in index:
                        file = os.path.abspath(os.path.join(directory, file_in_dir))
                        if file in index_entries:
                            # ObservationInfo overrides raw metadata
                            if isinstance(index[file_in_dir], ObservationInfo) \
                                    and not isinstance(index_entries[file], ObservationInfo):
                                self.log.warning("File %s already specified in an index file but overriding"
                                                 " with ObservationInfo content from %s",
                                                 file, possible_index_file)
                            else:
                                self.log.warning("File %s already specified in an index file, "
                                                 "ignoring content from %s", file, possible_index_file)
                                # Do nothing in this case
                                continue

                        index_entries[file] = index[file_in_dir]

        # Remove files from list that have index entries and also
        # any files that we determined to be explicit index files
        # or any index files that we failed to read.
        files = files_set - set(index_entries) - good_index_files - bad_index_files

        return index_entries, files, good_index_files, bad_index_files

    def processIndexEntries(self, index_entries):
        """Convert index entries to RawFileData.

        Parameters
        ----------
        index_entries : `dict` [`str`, Any]
            Dict indexed by name of file to ingest and with keys either
            raw metadata or translated `ObservationInfo`.

        Returns
        -------
        data :  `RawFileData`
            A structure containing the metadata extracted from the file,
            as well as the original filename.  All fields will be populated,
            but the `RawFileData.dataId` attribute will be a minimal
            (unexpanded) `~lsst.daf.butler.DataCoordinate` instance.
        """
        fileData = []
        for filename, metadata in index_entries.items():
            datasets = [self._calculate_dataset_info(metadata, filename)]
            instrument, formatterClass = self._determine_instrument_formatter(datasets[0].dataId, filename)
            if instrument is None:
                datasets = []
            fileData.append(RawFileData(datasets=datasets, filename=filename,
                                        FormatterClass=formatterClass, instrumentClass=instrument))
        return fileData

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
            `~lsst.daf.butler.DataCoordinate` instances.
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
            updated to data IDs for which
            `~lsst.daf.butler.DataCoordinate.hasRecords` returns `True`.
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

    def prep(self, files, *, pool: Optional[Pool] = None, processes: int = 1
             ) -> Tuple[Iterator[RawExposureData], List[str]]:
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

        Returns
        -------
        exposures : `Iterator` [ `RawExposureData` ]
            Data structures containing dimension records, filenames, and data
            IDs to be ingested (one structure for each exposure).
        bad_files : `list` of `str`
            List of all the files that could not have metadata extracted.
        """
        if pool is None and processes > 1:
            pool = Pool(processes)
        mapFunc = map if pool is None else pool.imap_unordered

        def _partition_good_bad(file_data: Iterable[RawFileData]) -> Tuple[List[RawFileData], List[str]]:
            """Filter out bad files and return good with list of bad."""
            good_files = []
            bad_files = []
            for fileDatum in file_data:
                if not fileDatum.datasets:
                    bad_files.append(fileDatum.filename)
                else:
                    good_files.append(fileDatum)
            return good_files, bad_files

        # Look for index files and read them
        # There should be far fewer index files than data files
        index_entries, files, good_index_files, bad_index_files = self.locateAndReadIndexFiles(files)
        if bad_index_files:
            self.log.info("Failed to read the following index files:"),
            for bad in sorted(bad_index_files):
                self.log.info("- %s", bad)

        # Now convert all the index file entries to standard form for ingest
        indexFileData = self.processIndexEntries(index_entries)
        if indexFileData:
            indexFileData, bad_index_file_data = _partition_good_bad(indexFileData)
            self.log.info("Successfully extracted metadata for %d file%s found in %d index file%s"
                          " with %d failure%s",
                          *_log_msg_counter(indexFileData),
                          *_log_msg_counter(good_index_files),
                          *_log_msg_counter(bad_index_file_data))

        # Extract metadata and build per-detector regions.
        # This could run in a subprocess so collect all output
        # before looking at failures.
        fileData: Iterator[RawFileData] = mapFunc(self.extractMetadata, files)

        # Filter out all the failed reads and store them for later
        # reporting
        fileData, bad_files = _partition_good_bad(fileData)
        self.log.info("Successfully extracted metadata from %d file%s with %d failure%s",
                      *_log_msg_counter(fileData),
                      *_log_msg_counter(bad_files))

        # Combine with data from index files
        if indexFileData:
            fileData.extend(indexFileData)
            bad_files.extend(bad_index_file_data)

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
        return mapFunc(self.expandDataIds, exposureData), bad_files

    def ingestExposureDatasets(self, exposure: RawExposureData, *, run: Optional[str] = None
                               ) -> List[FileDataset]:
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
        datasets : `list` of `lsst.daf.butler.FileDataset`
            Per-file structures identifying the files ingested and their
            dataset representation in the data repository.
        """
        datasets = [FileDataset(path=os.path.abspath(file.filename),
                                refs=[DatasetRef(self.datasetType, d.dataId) for d in file.datasets],
                                formatter=file.FormatterClass)
                    for file in exposure.files]
        self.butler.ingest(*datasets, transfer=self.config.transfer, run=run)
        return datasets

    @timeMethod
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
        exposureData, bad_files = self.prep(files, pool=pool, processes=processes)
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
        n_exposures = 0
        n_exposures_failed = 0
        n_ingests_failed = 0
        for exposure in exposureData:

            self.log.debug("Attempting to ingest %d file%s from exposure %s:%s",
                           *_log_msg_counter(exposure.files),
                           exposure.record.instrument, exposure.record.obs_id)

            try:
                self.butler.registry.syncDimensionData("exposure", exposure.record)
            except Exception as e:
                self._on_ingest_failure(exposure, e)
                n_exposures_failed += 1
                self.log.warning("Exposure %s:%s could not be registered: %s",
                                 exposure.record.instrument, exposure.record.obs_id, e)
                if self.config.failFast:
                    raise e
                continue

            # Override default run if nothing specified explicitly
            if run is None:
                instrumentClass = exposure.files[0].instrumentClass
                this_run = instrumentClass.makeDefaultRawIngestRunName()
            else:
                this_run = run
            if this_run not in runs:
                self.butler.registry.registerCollection(this_run, type=CollectionType.RUN)
                runs.add(this_run)
            try:
                with self.butler.transaction():
                    datasets_for_exposure = self.ingestExposureDatasets(exposure, run=this_run)
            except Exception as e:
                self._on_ingest_failure(exposure, e)
                n_ingests_failed += 1
                self.log.warning("Failed to ingest the following for reason: %s", e)
                for f in exposure.files:
                    self.log.warning("- %s", f.filename)
                if self.config.failFast:
                    raise e
                continue
            else:
                self._on_success(datasets_for_exposure)
                for dataset in datasets_for_exposure:
                    refs.extend(dataset.refs)

            # Success for this exposure
            n_exposures += 1
            self.log.info("Exposure %s:%s ingested successfully",
                          exposure.record.instrument, exposure.record.obs_id)

        had_failure = False

        if bad_files:
            had_failure = True
            self.log.warning("Could not extract observation metadata from the following:")
            for f in bad_files:
                self.log.warning("- %s", f)

        self.log.info("Successfully processed data from %d exposure%s with %d failure%s from exposure"
                      " registration and %d failure%s from file ingest.",
                      *_log_msg_counter(n_exposures),
                      *_log_msg_counter(n_exposures_failed),
                      *_log_msg_counter(n_ingests_failed))
        if n_exposures_failed > 0 or n_ingests_failed > 0:
            had_failure = True
        self.log.info("Ingested %d distinct Butler dataset%s", *_log_msg_counter(refs))

        if had_failure:
            raise RuntimeError("Some failures encountered during ingestion")

        return refs
