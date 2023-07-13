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

import json
import re
from collections import defaultdict
from collections.abc import Callable, Iterable, Iterator, MutableMapping, Sized
from dataclasses import InitVar, dataclass
from multiprocessing import Pool
from typing import Any, ClassVar

from astro_metadata_translator import MetadataTranslator, ObservationInfo, merge_headers
from astro_metadata_translator.indexing import process_index_data, process_sidecar_data
from lsst.afw.fits import readMetadata
from lsst.daf.butler import (
    Butler,
    CollectionType,
    DataCoordinate,
    DatasetIdGenEnum,
    DatasetRef,
    DatasetType,
    DimensionRecord,
    DimensionUniverse,
    FileDataset,
    Formatter,
    Progress,
)
from lsst.pex.config import ChoiceField, Config, Field
from lsst.pipe.base import Instrument, Task
from lsst.resources import ResourcePath, ResourcePathExpression
from lsst.utils.timer import timeMethod

from ._instrument import makeExposureRecordFromObsInfo

# multiprocessing.Pool is actually a function, not a type, and the real type
# isn't exposed, so we can't used it annotations, so we'll just punt on it via
# this alias instead.
PoolType = Any


def _do_nothing(*args: Any, **kwargs: Any) -> None:
    """Do nothing.

    This is a function that accepts anything and does nothing.
    For use as a default in callback arguments.
    """
    pass


def _log_msg_counter(noun: int | Sized) -> tuple[int, str]:
    """Count the iterable and return the count and plural modifier.

    Parameters
    ----------
    noun : `Sized` or `int`
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

    Examples
    --------
    .. code-block:: python

       log.warning("Found %d file%s", *_log_msg_counter(nfiles))
    """
    if isinstance(noun, int):
        num = noun
    else:
        num = len(noun)
    return num, "" if num == 1 else "s"


@dataclass
class RawFileDatasetInfo:
    """Information about a single dataset within a raw file."""

    dataId: DataCoordinate
    """Data ID for this file (`lsst.daf.butler.DataCoordinate`)."""

    obsInfo: ObservationInfo
    """Standardized observation metadata extracted directly from the file
    headers (`astro_metadata_translator.ObservationInfo`).
    """


@dataclass
class RawFileData:
    """Information about a single raw file, used during ingest."""

    datasets: list[RawFileDatasetInfo]
    """The information describing each dataset within this raw file.
    (`list` of `RawFileDatasetInfo`)
    """

    filename: ResourcePath
    """URI of the file this information was extracted from (`str`).

    This is the path prior to ingest, not the path after ingest.
    """

    FormatterClass: type[Formatter]
    """Formatter class that should be used to ingest this file (`type`; as
    subclass of `~lsst.daf.butler.Formatter`).
    """

    instrument: Instrument | None
    """The `Instrument` instance associated with this file. Can be `None`
    if ``datasets`` is an empty list."""


@dataclass
class RawExposureData:
    """Information about a complete raw exposure, used during ingest."""

    dataId: DataCoordinate
    """Data ID for this exposure (`lsst.daf.butler.DataCoordinate`).
    """

    files: list[RawFileData]
    """List of structures containing file-level information.
    """

    universe: InitVar[DimensionUniverse]
    """Set of all known dimensions.
    """

    record: DimensionRecord
    """The exposure `DimensionRecord` that must be inserted into the
    `~lsst.daf.butler.Registry` prior to file-level ingest
    (`~lsst.daf.butler.DimensionRecord`).
    """

    dependencyRecords: dict[str, DimensionRecord]
    """Additional records that must be inserted into the
    `~lsst.daf.butler.Registry` prior to ingesting the exposure ``record``
    (e.g., to satisfy foreign key constraints), indexed by the dimension name.
    """


def makeTransferChoiceField(
    doc: str = "How to transfer files (None for no transfer).", default: str = "auto"
) -> ChoiceField:
    """Create a Config field with options for transferring data between repos.

    The allowed options for the field are exactly those supported by
    `lsst.daf.butler.Datastore.ingest`.

    Parameters
    ----------
    doc : `str`
        Documentation for the configuration field.
    default : `str`, optional
        Default transfer mode for the field.

    Returns
    -------
    field : `lsst.pex.config.ChoiceField`
        Configuration field.
    """
    return ChoiceField(
        doc=doc,
        dtype=str,
        allowed={
            "move": "move",
            "copy": "copy",
            "auto": "choice will depend on datastore",
            "direct": "use URI to ingested file directly in datastore",
            "link": "hard link falling back to symbolic link",
            "hardlink": "hard link",
            "symlink": "symbolic (soft) link",
            "relsymlink": "relative symbolic link",
        },
        optional=True,
        default=default,
    )


class RawIngestConfig(Config):
    """Configuration class for RawIngestTask."""

    transfer = makeTransferChoiceField()
    failFast: Field[bool] = Field(
        dtype=bool,
        default=False,
        doc="If True, stop ingest as soon as any problem is encountered with any file. "
        "Otherwise problem files will be skipped and logged and a report issued at completion.",
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
        metadata for a file.  Will be passed the URI and the exception, in
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

    ConfigClass: ClassVar[type[Config]] = RawIngestConfig

    _DefaultName: ClassVar[str] = "ingest"

    def getDatasetType(self) -> DatasetType:
        """Return the default DatasetType of the datasets ingested by this
        Task.

        Returns
        -------
        datasetType : `DatasetType`
            The default dataset type to use for the data being ingested. This
            is only used if the relevant `~lsst.pipe.base.Instrument` does not
            define an override.
        """
        return DatasetType(
            "raw",
            ("instrument", "detector", "exposure"),
            "Exposure",
            universe=self.butler.dimensions,
        )

    # Mypy can not determine that the config passed to super() is this type.
    config: RawIngestConfig

    def __init__(
        self,
        config: RawIngestConfig,
        *,
        butler: Butler,
        on_success: Callable[[list[FileDataset]], Any] = _do_nothing,
        on_metadata_failure: Callable[[ResourcePath, Exception], Any] = _do_nothing,
        on_ingest_failure: Callable[[RawExposureData, Exception], Any] = _do_nothing,
        **kwargs: Any,
    ):
        config.validate()  # Not a CmdlineTask nor PipelineTask, so have to validate the config here.
        super().__init__(config, **kwargs)
        self.butler = butler
        self.universe = self.butler.dimensions
        self.datasetType = self.getDatasetType()
        self._on_success = on_success
        self._on_metadata_failure = on_metadata_failure
        self._on_ingest_failure = on_ingest_failure
        self.progress = Progress("obs.base.RawIngestTask")

        # Import all the instrument classes so that we ensure that we
        # have all the relevant metadata translators loaded.
        Instrument.importAll(self.butler.registry)

    def _reduce_kwargs(self) -> dict[str, Any]:
        # Add extra parameters to pickle.
        return dict(
            **super()._reduce_kwargs(),
            butler=self.butler,
            on_success=self._on_success,
            on_metadata_failure=self._on_metadata_failure,
            on_ingest_failure=self._on_ingest_failure,
        )

    def _determine_instrument_formatter(
        self, dataId: DataCoordinate, filename: ResourcePath
    ) -> tuple[Instrument | None, type[Formatter]]:
        """Determine the instrument and formatter class.

        Parameters
        ----------
        dataId : `lsst.daf.butler.DataCoordinate`
            The dataId associated with this dataset.
        filename : `lsst.resources.ResourcePath`
            URI of file used for error reporting.

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
            instrument = Instrument.fromName(dataId["instrument"], self.butler.registry)  # type: ignore
        except LookupError as e:
            self._on_metadata_failure(filename, e)
            self.log.warning(
                "Instrument %s for file %s not known to registry", dataId["instrument"], filename
            )
            if self.config.failFast:
                raise RuntimeError(
                    f"Instrument {dataId['instrument']} for file {filename} not known to registry"
                ) from e
            FormatterClass = Formatter
            # Indicate that we could not work out the instrument.
            instrument = None
        else:
            assert instrument is not None, "Should be guaranted by fromName succeeding."
            FormatterClass = instrument.getRawFormatter(dataId)
        return instrument, FormatterClass

    def extractMetadata(self, filename: ResourcePath) -> RawFileData:
        """Extract and process metadata from a single raw file.

        Parameters
        ----------
        filename : `lsst.resources.ResourcePath`
            URI to the file.

        Returns
        -------
        data : `RawFileData`
            A structure containing the metadata extracted from the file,
            as well as the original filename.  All fields will be populated,
            but the `RawFileData.dataId` attribute will be a minimal
            (unexpanded) `~lsst.daf.butler.DataCoordinate` instance. The
            ``instrument`` field will be `None` if there is a problem
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
        sidecar_fail_msg = ""  # Requires prepended space when set.
        try:
            sidecar_file = filename.updatedExtension(".json")
            if sidecar_file.exists():
                content = json.loads(sidecar_file.read())
                headers = [process_sidecar_data(content)]
                sidecar_fail_msg = " (via sidecar)"
            else:
                # Read the metadata from the data file itself.

                # For remote files download the entire file to get the
                # header. This is very inefficient and it would be better
                # to have some way of knowing where in the file the headers
                # are and to only download those parts of the file.
                with filename.as_local() as local_file:
                    # Read the primary. This might be sufficient.
                    header = readMetadata(local_file.ospath, 0)

                    try:
                        # Try to work out a translator class early.
                        translator_class = MetadataTranslator.determine_translator(
                            header, filename=str(filename)
                        )
                    except ValueError:
                        # Primary header was not sufficient (maybe this file
                        # has been compressed or is a MEF with minimal
                        # primary). Read second header and merge with primary.
                        header = merge_headers([header, readMetadata(local_file.ospath, 1)], mode="overwrite")

                    # Try again to work out a translator class, letting this
                    # fail.
                    translator_class = MetadataTranslator.determine_translator(header, filename=str(filename))

                    # Request the headers to use for ingest
                    headers = list(translator_class.determine_translatable_headers(local_file.ospath, header))

            # Add each header to the dataset list
            datasets = [self._calculate_dataset_info(h, filename) for h in headers]

        except Exception as e:
            self.log.debug("Problem extracting metadata from %s%s: %s", filename, sidecar_fail_msg, e)
            # Indicate to the caller that we failed to read.
            datasets = []
            formatterClass = Formatter
            instrument = None
            self._on_metadata_failure(filename, e)
            if self.config.failFast:
                raise RuntimeError(
                    f"Problem extracting metadata for file {filename}{sidecar_fail_msg}"
                ) from e
        else:
            self.log.debug("Extracted metadata for file %s%s", filename, sidecar_fail_msg)
            # The data model currently assumes that whilst multiple datasets
            # can be associated with a single file, they must all share the
            # same formatter.
            instrument, formatterClass = self._determine_instrument_formatter(datasets[0].dataId, filename)
            if instrument is None:
                datasets = []

        return RawFileData(
            datasets=datasets,
            filename=filename,
            # MyPy wants this to be a non-abstract class, which is not true
            # for the error case where instrument is None and datasets=[].
            FormatterClass=formatterClass,  # type: ignore
            instrument=instrument,
        )

    @classmethod
    def getObservationInfoSubsets(cls) -> tuple[set, set]:
        """Return subsets of fields in the `ObservationInfo` that we care
        about.

        These fields will be used in constructing an exposure record.

        Returns
        -------
        required : `set`
            Set of `ObservationInfo` field names that are required.
        optional : `set`
            Set of `ObservationInfo` field names we will use if they are
            available.
        """
        # Marking the new properties "group_counter_*" and
        # "has_simulated_content" as required, assumes that we either
        # recreate any existing index/sidecar files that include translated
        # values, or else allow astro_metadata_translator to fill in
        # defaults.
        required = {
            "datetime_begin",
            "datetime_end",
            "detector_num",
            "exposure_id",
            "exposure_time",
            "group_counter_end",
            "group_counter_start",
            "has_simulated_content",
            "instrument",
            "observation_id",
            "observation_type",
            "physical_filter",
        }
        optional = {
            "altaz_begin",
            "boresight_rotation_coord",
            "boresight_rotation_angle",
            "dark_time",
            "exposure_group",
            "tracking_radec",
            "object",
            "observation_counter",
            "observation_reason",
            "observing_day",
            "science_program",
            "visit_id",
        }
        return required, optional

    def _calculate_dataset_info(
        self, header: MutableMapping[str, Any] | ObservationInfo, filename: ResourcePath
    ) -> RawFileDatasetInfo:
        """Calculate a RawFileDatasetInfo from the supplied information.

        Parameters
        ----------
        header : Mapping or `astro_metadata_translator.ObservationInfo`
            Header from the dataset or previously-translated content.
        filename : `lsst.resources.ResourcePath`
            Filename to use for error messages.

        Returns
        -------
        dataset : `RawFileDatasetInfo`
            The dataId, and observation information associated with this
            dataset.
        """
        required, optional = self.getObservationInfoSubsets()
        if isinstance(header, ObservationInfo):
            obsInfo = header
            missing = []
            # Need to check the required properties are present.
            for property in required:
                # getattr does not need to be protected because it is using
                # the defined list above containing properties that must exist.
                value = getattr(obsInfo, property)
                if value is None:
                    missing.append(property)
            if missing:
                raise ValueError(
                    f"Requested required properties are missing from file {filename}: {missing} (via JSON)"
                )

        else:
            obsInfo = ObservationInfo(
                header,
                pedantic=False,
                filename=str(filename),
                required=required,
                subset=required | optional,
            )

        dataId = DataCoordinate.standardize(
            instrument=obsInfo.instrument,
            exposure=obsInfo.exposure_id,
            detector=obsInfo.detector_num,
            universe=self.universe,
        )
        return RawFileDatasetInfo(obsInfo=obsInfo, dataId=dataId)

    def locateAndReadIndexFiles(
        self, files: Iterable[ResourcePath]
    ) -> tuple[dict[ResourcePath, Any], list[ResourcePath], set[ResourcePath], set[ResourcePath]]:
        """Given a list of files, look for index files and read them.

        Index files can either be explicitly in the list of files to
        ingest, or else located in the same directory as a file to ingest.
        Index entries are always used if present.

        Parameters
        ----------
        files : iterable over `lsst.resources.ResourcePath`
            URIs to the files to be ingested.

        Returns
        -------
        index : `dict` [`ResourcePath`, Any]
            Merged contents of all relevant index files found. These can
            be explicitly specified index files or ones found in the
            directory alongside a data file to be ingested.
        updated_files : `list` of `ResourcePath`
            Updated list of the input files with entries removed that were
            found listed in an index file. Order is not guaranteed to
            match the order of the files given to this routine.
        good_index_files: `set` [ `ResourcePath` ]
            Index files that were successfully read.
        bad_index_files: `set` [ `ResourcePath` ]
            Files that looked like index files but failed to read properly.
        """
        # Convert the paths to absolute for easy comparison with index content.
        # Do not convert to real paths since we have to assume that index
        # files are in this location and not the location which it links to.
        files = tuple(f.abspath() for f in files)

        # Index files must be named this.
        index_root_file = "_index.json"

        # Group the files by directory.
        files_by_directory = defaultdict(set)

        for path in files:
            directory, file_in_dir = path.split()
            files_by_directory[directory].add(file_in_dir)

        # All the metadata read from index files with keys of full path.
        index_entries: dict[ResourcePath, Any] = {}

        # Index files we failed to read.
        bad_index_files = set()

        # Any good index files that were found and used.
        good_index_files = set()

        # Look for index files in those directories.
        for directory, files_in_directory in files_by_directory.items():
            possible_index_file = directory.join(index_root_file)
            if possible_index_file.exists():
                # If we are explicitly requesting an index file the
                # messages should be different.
                index_msg = "inferred"
                is_implied = True
                if index_root_file in files_in_directory:
                    index_msg = "explicit"
                    is_implied = False

                # Try to read the index file and catch and report any
                # problems.
                try:
                    content = json.loads(possible_index_file.read())
                    index = process_index_data(content, force_dict=True)
                    # mypy should in theory know that this is a mapping
                    # from the overload type annotation of process_index_data.
                    assert isinstance(index, MutableMapping)
                except Exception as e:
                    # Only trigger the callback if the index file
                    # was asked for explicitly. Triggering on implied file
                    # might be surprising.
                    if not is_implied:
                        self._on_metadata_failure(possible_index_file, e)
                    if self.config.failFast:
                        raise RuntimeError(
                            f"Problem reading index file from {index_msg} location {possible_index_file}"
                        ) from e
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
                    # This should never happen because an explicit index
                    # file will force ingest of all files in the index
                    # and not use the explicit file list. If somehow
                    # this is not true we continue. Raising an exception
                    # seems like the wrong thing to do since this is harmless.
                    if file_in_dir == index_root_file:
                        self.log.info(
                            "Logic error found scanning directory %s. Please file ticket.", directory
                        )
                        continue
                    if file_in_dir in index:
                        file = directory.join(file_in_dir)
                        if file in index_entries:
                            # ObservationInfo overrides raw metadata
                            if isinstance(index[file_in_dir], ObservationInfo) and not isinstance(
                                index_entries[file], ObservationInfo
                            ):
                                self.log.warning(
                                    "File %s already specified in an index file but overriding"
                                    " with ObservationInfo content from %s",
                                    file,
                                    possible_index_file,
                                )
                            else:
                                self.log.warning(
                                    "File %s already specified in an index file, ignoring content from %s",
                                    file,
                                    possible_index_file,
                                )
                                # Do nothing in this case
                                continue

                        index_entries[file] = index[file_in_dir]

        # Remove files from list that have index entries and also
        # any files that we determined to be explicit index files
        # or any index files that we failed to read.
        filtered = set(files) - set(index_entries) - good_index_files - bad_index_files

        # The filtered list loses the initial order. Retaining the order
        # is good for testing but does have a cost if there are many
        # files when copying the good values out. A dict would have faster
        # lookups (using the files as keys) but use more memory.
        ordered = [f for f in filtered if f in files]

        return index_entries, ordered, good_index_files, bad_index_files

    def processIndexEntries(self, index_entries: dict[ResourcePath, Any]) -> list[RawFileData]:
        """Convert index entries to RawFileData.

        Parameters
        ----------
        index_entries : `dict` [`ResourcePath`, Any]
            Dict indexed by name of file to ingest and with keys either
            raw metadata or translated
            `~astro_metadata_translator.ObservationInfo`.

        Returns
        -------
        data : `list` [ `RawFileData` ]
            Structures containing the metadata extracted from the file,
            as well as the original filename.  All fields will be populated,
            but the `RawFileData.dataId` attributes will be minimal
            (unexpanded) `~lsst.daf.butler.DataCoordinate` instances.
        """
        fileData = []
        for filename, metadata in index_entries.items():
            try:
                datasets = [self._calculate_dataset_info(metadata, filename)]
            except Exception as e:
                self.log.debug("Problem extracting metadata for file %s found in index file: %s", filename, e)
                datasets = []
                formatterClass = Formatter
                instrument = None
                self._on_metadata_failure(filename, e)
                if self.config.failFast:
                    raise RuntimeError(
                        f"Problem extracting metadata for file {filename} found in index file"
                    ) from e
            else:
                instrument, formatterClass = self._determine_instrument_formatter(
                    datasets[0].dataId, filename
                )
                if instrument is None:
                    datasets = []
            fileData.append(
                RawFileData(
                    datasets=datasets,
                    filename=filename,
                    # MyPy wants this to be a non-abstract class, which is not
                    # true for the error case where instrument is None and
                    # datasets=[].
                    FormatterClass=formatterClass,  # type: ignore
                    instrument=instrument,
                )
            )
        return fileData

    def groupByExposure(self, files: Iterable[RawFileData]) -> list[RawExposureData]:
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
            # Assume that the first dataset is representative for the file.
            byExposure[f.datasets[0].dataId.subset(exposureDimensions)].append(f)

        return [
            RawExposureData(
                dataId=dataId,
                files=exposureFiles,
                universe=self.universe,
                record=self.makeExposureRecord(exposureFiles[0].datasets[0].obsInfo, self.universe),
                dependencyRecords=self.makeDependencyRecords(
                    exposureFiles[0].datasets[0].obsInfo, self.universe
                ),
            )
            for dataId, exposureFiles in byExposure.items()
        ]

    def makeExposureRecord(
        self, obsInfo: ObservationInfo, universe: DimensionUniverse, **kwargs: Any
    ) -> DimensionRecord:
        """Construct a registry record for an exposure.

        This is a method that subclasses will often want to customize. This can
        often be done by calling this base class implementation with additional
        ``kwargs``.

        Parameters
        ----------
        obsInfo : `ObservationInfo`
            Observation details for (one of the components of) the exposure.
        universe : `DimensionUniverse`
            Set of all known dimensions.
        **kwargs
            Additional field values for this record.

        Returns
        -------
        record : `DimensionRecord`
            The exposure record that must be inserted into the
            `~lsst.daf.butler.Registry` prior to file-level ingest.
        """
        return makeExposureRecordFromObsInfo(obsInfo, universe, **kwargs)

    def makeDependencyRecords(
        self, obsInfo: ObservationInfo, universe: DimensionUniverse
    ) -> dict[str, DimensionRecord]:
        """Construct dependency records.

        These dependency records will be inserted into the
        `~lsst.daf.butler.Registry` before the exposure records, because they
        are dependencies of the exposure. This allows an opportunity to satisfy
        foreign key constraints that exist because of dimensions related to the
        exposure.

        This is a method that subclasses may want to customize, if they've
        added dimensions that relate to an exposure.

        Parameters
        ----------
        obsInfo : `ObservationInfo`
            Observation details for (one of the components of) the exposure.
        universe : `DimensionUniverse`
            Set of all known dimensions.

        Returns
        -------
        records : `dict` [`str`, `DimensionRecord`]
            The records to insert, indexed by dimension name.
        """
        return {}

    def expandDataIds(self, data: RawExposureData) -> RawExposureData:
        """Expand the data IDs associated with a raw exposure.

        This adds the metadata records.

        Parameters
        ----------
        exposure : `RawExposureData`
            A structure containing information about the exposure to be
            ingested.  Must have `RawExposureData.record` populated. Should
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
            records={"exposure": data.record},
        )
        # Now we expand the per-file (exposure+detector) data IDs.  This time
        # we pass in the records we just retrieved from the exposure data ID
        # expansion.
        for file in data.files:
            for dataset in file.datasets:
                dataset.dataId = self.butler.registry.expandDataId(
                    dataset.dataId, records=data.dataId.records
                )
        return data

    def prep(
        self, files: Iterable[ResourcePath], *, pool: PoolType | None = None
    ) -> tuple[Iterator[RawExposureData], list[ResourcePath]]:
        """Perform all non-database-updating ingest preprocessing steps.

        Parameters
        ----------
        files : iterable over `str` or path-like objects
            Paths to the files to be ingested.  Will be made absolute
            if they are not already.
        pool : `multiprocessing.Pool`, optional
            If not `None`, a process pool with which to parallelize some
            operations.

        Returns
        -------
        exposures : `Iterator` [ `RawExposureData` ]
            Data structures containing dimension records, filenames, and data
            IDs to be ingested (one structure for each exposure).
        bad_files : `list` of `str`
            List of all the files that could not have metadata extracted.
        """
        mapFunc = map if pool is None else pool.imap_unordered

        def _partition_good_bad(
            file_data: Iterable[RawFileData],
        ) -> tuple[list[RawFileData], list[ResourcePath]]:
            """Filter out bad files and return good with list of bad."""
            good_files = []
            bad_files = []
            for fileDatum in self.progress.wrap(file_data, desc="Reading image metadata"):
                if not fileDatum.datasets:
                    bad_files.append(fileDatum.filename)
                else:
                    good_files.append(fileDatum)
            return good_files, bad_files

        # Look for index files and read them.
        # There should be far fewer index files than data files.
        index_entries, files, good_index_files, bad_index_files = self.locateAndReadIndexFiles(files)
        if bad_index_files:
            self.log.info("Failed to read the following explicitly requested index files:")
            for bad in sorted(bad_index_files):
                self.log.info("- %s", bad)

        # Now convert all the index file entries to standard form for ingest.
        processed_bad_index_files: list[ResourcePath] = []
        indexFileData = self.processIndexEntries(index_entries)
        if indexFileData:
            indexFileData, processed_bad_index_files = _partition_good_bad(indexFileData)
            self.log.info(
                "Successfully extracted metadata for %d file%s found in %d index file%s with %d failure%s",
                *_log_msg_counter(indexFileData),
                *_log_msg_counter(good_index_files),
                *_log_msg_counter(processed_bad_index_files),
            )

        # Extract metadata and build per-detector regions.
        # This could run in a subprocess so collect all output
        # before looking at failures.
        fileData: Iterator[RawFileData] = mapFunc(self.extractMetadata, files)

        # Filter out all the failed reads and store them for later
        # reporting.
        good_file_data, bad_files = _partition_good_bad(fileData)
        self.log.info(
            "Successfully extracted metadata from %d file%s with %d failure%s",
            *_log_msg_counter(good_file_data),
            *_log_msg_counter(bad_files),
        )

        # Combine with data from index files.
        good_file_data.extend(indexFileData)
        bad_files.extend(processed_bad_index_files)
        bad_files.extend(bad_index_files)

        # Use that metadata to group files (and extracted metadata) by
        # exposure.  Never parallelized because it's intrinsically a gather
        # step.
        exposureData: list[RawExposureData] = self.groupByExposure(good_file_data)

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

    def ingestExposureDatasets(
        self,
        exposure: RawExposureData,
        datasetType: DatasetType,
        *,
        run: str,
        skip_existing_exposures: bool = False,
        track_file_attrs: bool = True,
    ) -> list[FileDataset]:
        """Ingest all raw files in one exposure.

        Parameters
        ----------
        exposure : `RawExposureData`
            A structure containing information about the exposure to be
            ingested.  Must have `RawExposureData.records` populated and all
            data ID attributes expanded.
        datasetType : `DatasetType`
            The dataset type associated with this exposure.
        run : `str`
            Name of a RUN-type collection to write to.
        skip_existing_exposures : `bool`, optional
            If `True` (`False` is default), skip raws that have already been
            ingested (i.e. raws for which we already have a dataset with the
            same data ID in the target collection, even if from another file).
            Note that this is much slower than just not passing
            already-ingested files as inputs, because we still need to read and
            process metadata to identify which exposures to search for.  It
            also will not work reliably if multiple processes are attempting to
            ingest raws from the same exposure concurrently, in that different
            processes may still attempt to ingest the same raw and conflict,
            causing a failure that prevents other raws from the same exposure
            from being ingested.
        track_file_attrs : `bool`, optional
            Control whether file attributes such as the size or checksum should
            be tracked by the datastore. Whether this parameter is honored
            depends on the specific datastore implementation.

        Returns
        -------
        datasets : `list` of `lsst.daf.butler.FileDataset`
            Per-file structures identifying the files ingested and their
            dataset representation in the data repository.
        """
        if skip_existing_exposures:
            existing = {
                ref.dataId
                for ref in self.butler.registry.queryDatasets(
                    datasetType,
                    collections=[run],
                    dataId=exposure.dataId,
                )
            }
        else:
            existing = set()

        # Raw files are preferentially ingested using a UUID derived from
        # the collection name and dataId.
        if self.butler.registry.supportsIdGenerationMode(DatasetIdGenEnum.DATAID_TYPE_RUN):
            mode = DatasetIdGenEnum.DATAID_TYPE_RUN
        else:
            mode = DatasetIdGenEnum.UNIQUE

        datasets = []
        for file in exposure.files:
            refs = [
                DatasetRef(datasetType, d.dataId, run=run, id_generation_mode=mode)
                for d in file.datasets
                if d.dataId not in existing
            ]
            if refs:
                datasets.append(
                    FileDataset(path=file.filename.abspath(), refs=refs, formatter=file.FormatterClass)
                )

        self.butler.ingest(
            *datasets,
            transfer=self.config.transfer,
            record_validation_info=track_file_attrs,
        )
        return datasets

    def ingestFiles(
        self,
        files: Iterable[ResourcePath],
        *,
        pool: PoolType | None = None,
        processes: int = 1,
        run: str | None = None,
        skip_existing_exposures: bool = False,
        update_exposure_records: bool = False,
        track_file_attrs: bool = True,
    ) -> tuple[list[DatasetRef], list[ResourcePath], int, int, int]:
        """Ingest files into a Butler data repository.

        This creates any new exposure or visit Dimension entries needed to
        identify the ingested files, creates new Dataset entries in the
        Registry and finally ingests the files themselves into the Datastore.
        Any needed instrument, detector, and physical_filter Dimension entries
        must exist in the Registry before `run` is called.

        Parameters
        ----------
        files : iterable over `lsst.resources.ResourcePath`
            URIs to the files to be ingested.
        pool : `multiprocessing.Pool`, optional
            If not `None`, a process pool with which to parallelize some
            operations.
        processes : `int`, optional
            The number of processes to use.  Ignored if ``pool`` is not `None`.
        run : `str`, optional
            Name of a RUN-type collection to write to, overriding
            the default derived from the instrument name.
        skip_existing_exposures : `bool`, optional
            If `True` (`False` is default), skip raws that have already been
            ingested (i.e. raws for which we already have a dataset with the
            same data ID in the target collection, even if from another file).
            Note that this is much slower than just not passing
            already-ingested files as inputs, because we still need to read and
            process metadata to identify which exposures to search for.  It
            also will not work reliably if multiple processes are attempting to
            ingest raws from the same exposure concurrently, in that different
            processes may still attempt to ingest the same raw and conflict,
            causing a failure that prevents other raws from the same exposure
            from being ingested.
        update_exposure_records : `bool`, optional
            If `True` (`False` is default), update existing exposure records
            that conflict with the new ones instead of rejecting them.  THIS IS
            AN ADVANCED OPTION THAT SHOULD ONLY BE USED TO FIX METADATA THAT IS
            KNOWN TO BE BAD.  This should usually be combined with
            ``skip_existing_exposures=True``.
        track_file_attrs : `bool`, optional
            Control whether file attributes such as the size or checksum should
            be tracked by the datastore. Whether this parameter is honored
            depends on the specific datastore implentation.

        Returns
        -------
        refs : `list` of `lsst.daf.butler.DatasetRef`
            Dataset references for ingested raws.
        bad_files : `list` of `ResourcePath`
            Given paths that could not be ingested.
        n_exposures : `int`
            Number of exposures successfully ingested.
        n_exposures_failed : `int`
            Number of exposures that failed when inserting dimension data.
        n_ingests_failed : `int`
            Number of exposures that failed when ingesting raw datasets.
        """
        created_pool = False
        if pool is None and processes > 1:
            pool = Pool(processes)
            created_pool = True

        try:
            exposureData, bad_files = self.prep(files, pool=pool)
        finally:
            if created_pool and pool:
                # The pool is not needed any more so close it if we created
                # it to ensure we clean up resources.
                pool.close()
                pool.join()

        # Up to this point, we haven't modified the data repository at all.
        # Now we finally do that, with one transaction per exposure.  This is
        # not parallelized at present because the performance of this step is
        # limited by the database server.  That may or may not change in the
        # future once we increase our usage of bulk inserts and reduce our
        # usage of savepoints; we've tried to get everything but the database
        # operations done in advance to reduce the time spent inside
        # transactions.
        refs = []
        runs = set()
        datasetTypes: dict[str, DatasetType] = {}
        n_exposures = 0
        n_exposures_failed = 0
        n_ingests_failed = 0
        for exposure in self.progress.wrap(exposureData, desc="Ingesting raw exposures"):
            assert exposure.record is not None, "Should be guaranteed by prep()"
            self.log.debug(
                "Attempting to ingest %d file%s from exposure %s:%s",
                *_log_msg_counter(exposure.files),
                exposure.record.instrument,
                exposure.record.obs_id,
            )

            try:
                for name, record in exposure.dependencyRecords.items():
                    self.butler.registry.syncDimensionData(name, record, update=update_exposure_records)
                inserted_or_updated = self.butler.registry.syncDimensionData(
                    "exposure",
                    exposure.record,
                    update=update_exposure_records,
                )
            except Exception as e:
                self._on_ingest_failure(exposure, e)
                n_exposures_failed += 1
                self.log.warning(
                    "Exposure %s:%s could not be registered: %s",
                    exposure.record.instrument,
                    exposure.record.obs_id,
                    e,
                )
                if self.config.failFast:
                    raise e
                continue

            if isinstance(inserted_or_updated, dict):
                # Exposure is in the registry and we updated it, so
                # syncDimensionData returned a dict.
                self.log.info(
                    "Exposure %s:%s was already present, but columns %s were updated.",
                    exposure.record.instrument,
                    exposure.record.obs_id,
                    str(list(inserted_or_updated.keys())),
                )

            # Determine the instrument so we can work out the dataset type.
            instrument = exposure.files[0].instrument
            assert (
                instrument is not None
            ), "file should have been removed from this list by prep if instrument could not be found"

            if raw_definition := getattr(instrument, "raw_definition", None):
                datasetTypeName, dimensions, storageClass = raw_definition
                if not (datasetType := datasetTypes.get(datasetTypeName)):
                    datasetType = DatasetType(
                        datasetTypeName, dimensions, storageClass, universe=self.butler.dimensions
                    )
            else:
                datasetType = self.datasetType
            if datasetType.name not in datasetTypes:
                self.butler.registry.registerDatasetType(datasetType)
                datasetTypes[datasetType.name] = datasetType

            # Override default run if nothing specified explicitly.
            if run is None:
                this_run = instrument.makeDefaultRawIngestRunName()
            else:
                this_run = run
            if this_run not in runs:
                self.butler.registry.registerCollection(this_run, type=CollectionType.RUN)
                runs.add(this_run)
            try:
                datasets_for_exposure = self.ingestExposureDatasets(
                    exposure,
                    datasetType=datasetType,
                    run=this_run,
                    skip_existing_exposures=skip_existing_exposures,
                    track_file_attrs=track_file_attrs,
                )
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

            # Success for this exposure.
            n_exposures += 1
            self.log.info(
                "Exposure %s:%s ingested successfully", exposure.record.instrument, exposure.record.obs_id
            )

        return refs, bad_files, n_exposures, n_exposures_failed, n_ingests_failed

    @timeMethod
    def run(
        self,
        files: Iterable[ResourcePathExpression],
        *,
        pool: PoolType | None = None,
        processes: int = 1,
        run: str | None = None,
        file_filter: str | re.Pattern = r"\.fit[s]?\b",
        group_files: bool = True,
        skip_existing_exposures: bool = False,
        update_exposure_records: bool = False,
        track_file_attrs: bool = True,
    ) -> list[DatasetRef]:
        """Ingest files into a Butler data repository.

        This creates any new exposure or visit Dimension entries needed to
        identify the ingested files, creates new Dataset entries in the
        Registry and finally ingests the files themselves into the Datastore.
        Any needed instrument, detector, and physical_filter Dimension entries
        must exist in the Registry before `run` is called.

        Parameters
        ----------
        files : iterable `lsst.resources.ResourcePath`, `str` or path-like
            Paths to the files to be ingested.  Can refer to directories.
            Will be made absolute if they are not already.
        pool : `multiprocessing.Pool`, optional
            If not `None`, a process pool with which to parallelize some
            operations.
        processes : `int`, optional
            The number of processes to use.  Ignored if ``pool`` is not `None`.
        run : `str`, optional
            Name of a RUN-type collection to write to, overriding
            the default derived from the instrument name.
        file_filter : `str` or `re.Pattern`, optional
            Pattern to use to discover files to ingest within directories.
            The default is to search for FITS files. The regex applies to
            files within the directory.
        group_files : `bool`, optional
            Group files by directory if they have been discovered in
            directories. Will not affect files explicitly provided.
        skip_existing_exposures : `bool`, optional
            If `True` (`False` is default), skip raws that have already been
            ingested (i.e. raws for which we already have a dataset with the
            same data ID in the target collection, even if from another file).
            Note that this is much slower than just not passing
            already-ingested files as inputs, because we still need to read and
            process metadata to identify which exposures to search for.  It
            also will not work reliably if multiple processes are attempting to
            ingest raws from the same exposure concurrently, in that different
            processes may still attempt to ingest the same raw and conflict,
            causing a failure that prevents other raws from the same exposure
            from being ingested.
        update_exposure_records : `bool`, optional
            If `True` (`False` is default), update existing exposure records
            that conflict with the new ones instead of rejecting them.  THIS IS
            AN ADVANCED OPTION THAT SHOULD ONLY BE USED TO FIX METADATA THAT IS
            KNOWN TO BE BAD.  This should usually be combined with
            ``skip_existing_exposures=True``.
        track_file_attrs : `bool`, optional
            Control whether file attributes such as the size or checksum should
            be tracked by the datastore. Whether this parameter is honored
            depends on the specific datastore implentation.

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
        the same exposure to be ingested in different runs.
        """
        refs = []
        bad_files = []
        n_exposures = 0
        n_exposures_failed = 0
        n_ingests_failed = 0
        if group_files:
            for group in ResourcePath.findFileResources(files, file_filter, group_files):
                new_refs, bad, n_exp, n_exp_fail, n_ingest_fail = self.ingestFiles(
                    group,
                    pool=pool,
                    processes=processes,
                    run=run,
                    skip_existing_exposures=skip_existing_exposures,
                    update_exposure_records=update_exposure_records,
                    track_file_attrs=track_file_attrs,
                )
                refs.extend(new_refs)
                bad_files.extend(bad)
                n_exposures += n_exp
                n_exposures_failed += n_exp_fail
                n_ingests_failed += n_ingest_fail
        else:
            refs, bad_files, n_exposures, n_exposures_failed, n_ingests_failed = self.ingestFiles(
                ResourcePath.findFileResources(files, file_filter, group_files),
                pool=pool,
                processes=processes,
                run=run,
                skip_existing_exposures=skip_existing_exposures,
                update_exposure_records=update_exposure_records,
            )

        had_failure = False

        if bad_files:
            had_failure = True
            self.log.warning("Could not extract observation metadata from the following:")
            for f in bad_files:
                self.log.warning("- %s", f)

        self.log.info(
            "Successfully processed data from %d exposure%s with %d failure%s from exposure"
            " registration and %d failure%s from file ingest.",
            *_log_msg_counter(n_exposures),
            *_log_msg_counter(n_exposures_failed),
            *_log_msg_counter(n_ingests_failed),
        )
        if n_exposures_failed > 0 or n_ingests_failed > 0:
            had_failure = True
        self.log.info("Ingested %d distinct Butler dataset%s", *_log_msg_counter(refs))

        if had_failure:
            raise RuntimeError("Some failures encountered during ingestion")

        return refs
