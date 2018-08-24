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


__all__ = ("RawIngestTask", "RawIngestConfig", "VisitInfoRawIngestTask")

import os.path
from abc import ABCMeta, abstractmethod

from lsst.afw.image import readMetadata
from lsst.daf.butler import DatasetType, StorageClassFactory, Run
from lsst.daf.butler.instrument import makeExposureEntryFromVisitInfo, makeVisitEntryFromVisitInfo
from lsst.pex.config import Config, Field, ChoiceField
from lsst.pipe.base import Task


class IngestConflictError(RuntimeError):
    pass


class RawIngestConfig(Config):
    transfer = ChoiceField(
        ("How to transfer files (None for no transfer)."),
        dtype=str,
        allowed={"move": "move",
                 "copy": "copy",
                 "hardlink": "hard link",
                 "symlink": "symbolic (soft) link"},
        optional=True,
    )
    conflict = ChoiceField(
        ("What to do if a raw Dataset with the same data ID as an "
         "ingested file already exists in the Butler's Collection."),
        dtype=str,
        allowed={"ignore": ("Do not add the new file to the Collection.  If "
                            "'stash' is not None, the new file will be "
                            "ingested into the stash Collection instead."),
                 "fail": ("Raise RuntimeError if a conflict is encountered "
                          "(which may then be caught if onError == 'continue')."),
                 },
        optional=False,
        default="ignore",
    )
    stash = Field(
        "Name of an alternate Collection to hold Datasets that lose conflicts.",
        dtype=str,
        default=None,
    )
    onError = ChoiceField(
        "What to do if an error (including fatal conflicts) occurs.",
        dtype=str,
        allowed={"continue": "Warn and continue with the next file.",
                 "break": ("Stop processing immediately, but leave "
                           "already-ingested datasets in the repository."),
                 "rollback": ("Stop processing and attempt to remove aleady-"
                              "ingested datasets from the repository."),
                 },
        optional=False,
        default="continue",
    )


class RawIngestTask(Task, metaclass=ABCMeta):
    """Driver Task for ingesting raw data into Gen3 Butler repositories.

    This Task is intended to be runnable from the command-line, but it doesn't
    meet the other requirements of CmdLineTask or PipelineTask, and wouldn't
    gain much from being one.  It also wouldn't really be appropriate as a
    subtask of a CmdLineTask or PipelineTask; it's a Task essentially just to
    leverage the logging and configurability functionality that provides.

    Each instance of `RawIngestTask` writes to the same Butler and maintains a
    cache of DataUnit entries that have already been added to or extracted
    from its Registry.  Each invocation of `RawIngestTask.run` ingests a list
    of files (possibly semi-atomically; see `RawIngestConfig.onError`).

    RawIngestTask should be subclassed to specialize ingest for the actual
    structure of raw data files produced by a particular camera. Subclasses
    must either provide populated `MetadataReader` instances in the
    `dataIdReader`, `visitReader`, and `exposureReader` class attributes, or
    alternate implementations of the `extractDataId`, `extractVisit`, and
    `extractExposure` methods that do not use those attributes (each
    attribute-method pair may be handled differently).  Subclasses may also
    wish to override `getFormatter` and/or (rarely) `getDatasetType`.  We do
    not anticipate overriding `run`, `ensureDataUnits`, `ingestFile`, or
    `processFile` to ever be necessary.

    Parameters
    ----------
    config : `RawIngestConfig`
        Configuration for whether/how to transfer files and how to handle
        conflicts and errors.
    butler : `~lsst.daf.butler.Butler`
        Butler instance.  Ingested Datasets will be created as part of
        ``butler.run`` and associated with its Collection.

    Other keyword arguments are forwarded to the Task base class constructor.
    """

    ConfigClass = RawIngestConfig

    _DefaultName = "ingest"

    @classmethod
    def getDatasetType(cls):
        """Return the DatasetType of the Datasets ingested by this Task.
        """
        return DatasetType("raw", ("Camera", "Sensor", "Exposure"),
                           StorageClassFactory().getStorageClass("Exposure"))

    def __init__(self, config=None, *, butler, **kwds):
        super().__init__(config, **kwds)
        self.butler = butler
        self.datasetType = self.getDatasetType()
        self.units = tuple(butler.registry.getDataUnitDefinition(k)
                           for k in ("Camera", "Sensor", "PhysicalFilter", "Visit", "Exposure", ))
        # Nested dictionary of form {<unit-name>: {<primary-key-tuple>: {<field>: <value>}}}, where:
        #  - <unit-name> is a DataUnit name (e.g. Camera, Exposure)
        #  - <primary-key-tuple> is a tuple of values that correspond to the [compound] primary
        #    key for that DataUnit.  (TODO: make these DataId objects on DM-15034).
        #  - <field> is the name of a column in the table for this DataUnit.
        #  - <value> is the value of that field.
        # The {<field>: <value>} dict is called an "entry" in this class and in Registry methods.
        self.unitEntryCache = {k.name: {} for k in self.units}
        # (Possibly) create a Run object for the "stash": where we put datasets
        # that lose conflicts.  Note that this doesn't actually add this Run
        # to the Registry; we only do that on first use.
        self.stashRun = Run(self.config.stash) if self.config.stash is not None else None

    def run(self, files):
        """Ingest files into a Butler data repository.

        This creates any new Exposure or Visit DataUnit entries needed to
        identify the ingested files, creates new Dataset entries in the
        Registry and finally ingests the files themselves into the Datastore.
        Any needed Camera, Sensor, and PhysicalFilter DataUnit entries must
        exist in the Registry before `run` is called.

        Parameters
        ----------
        files : iterable over `str` or path-like objects
            Paths to the files to be ingested.  Will be made absolute
            if they are not already.
        """
        self.butler.registry.registerDatasetType(self.getDatasetType())
        if self.config.onError == "rollback":
            with self.butler.transaction():
                for file in files:
                    self.processFile(os.path.abspath(file))
        elif self.config.onError == "break":
            for file in files:
                self.processFile(os.path.abspath(file))
        elif self.config.onError == "continue":
            for file in files:
                try:
                    self.processFile(os.path.abspath(file))
                except Exception as err:
                    self.log.warnf("Error processing '{}': {}", file, err)

    def readHeaders(self, file):
        """Read and return any relevant headers from the given file.

        The default implementation simply reads the header of the first
        non-empty HDU, so it always returns a single-element list.

        Parameters
        ----------
        file : `str` or path-like object
            Absolute path to the file to be ingested.

        Returns
        -------
        headers : `list` of `~lsst.daf.base.PropertyList`
            Single-element list containing the header of the first
            non-empty HDU.
        """
        return [readMetadata(file)]

    def ensureDataUnits(self, file):
        """Extract metadata from a raw file and add Exposure and Visit
        DataUnit entries.

        Any needed Camera, Sensor, and PhysicalFilter DataUnit entries must
        exist in the Registry before `run` is called.

        Parameters
        ----------
        file : `str` or path-like object
            Absolute path to the file to be ingested.

        Returns
        -------
        headers : `list` of `~lsst.daf.base.PropertyList`
            Result of calling `readHeaders`.
        dataId : `dict`
            Data ID dictionary, as returned by `extractDataId`.
        """
        headers = self.readHeaders(file)

        # Extract a dictionary with structure {<link-name>: <value>} where:
        #  - <link-name> is the name of a DataUnit link to the Dataset table,
        #    usually a DataUnit primary key field (e.g. 'camera' or 'visit').
        #  - <value> is the value of that field
        dataId = self.extractDataId(file, headers)
        dataId.setdefault("physical_filter", None)
        dataId.setdefault("visit", None)

        # Locate or extract additional DataUnit metadata, producing a nested
        # dict with structure {<unit-name>: {<field>: <value>}}.  This is the
        # same content as self.unitEntryCache, but without the middle layer,
        # because this contains only the entries associated with this
        # particular file.
        associatedUnitEntries = {}
        for unit in self.units:
            # Start by looking in the Task's cache of unit entries, which is keyed by a tuple.
            unitPrimaryKeyTuple = tuple(dataId[f] for f in unit.primaryKey)
            if any(v is None for v in unitPrimaryKeyTuple):
                # This DataUnit isn't actually applicable for this file; move
                # on. Could be a calibration Exposure that doesn't have a
                # Visit, for example.
                associatedUnitEntries[unit.name] = None
                continue
            unitEntryDict = self.unitEntryCache[unit.name].get(unitPrimaryKeyTuple, None)
            if unitEntryDict is None:
                # Next look in the Registry, which is keyed by a dataId-like dict
                unitPrimaryKeyDict = {f: dataId[f] for f in unit.primaryKey}
                unitEntryDict = self.butler.registry.findDataUnitEntry(unit.name, unitPrimaryKeyDict)
                if unitEntryDict is None:
                    # If we haven't found it, either raise an exception or extract that information
                    # from the headers (and possibly the filename).
                    if unit.name == "Visit":
                        extractMethod = self.extractVisitEntry
                    elif unit.name == "Exposure":
                        extractMethod = self.extractExposureEntry
                    else:
                        raise LookupError("{} with keys {} not found; must be present in Registry prior "
                                          "to ingest.".format(unit.name, unitPrimaryKeyDict))
                    unitEntryDict = extractMethod(file, headers, dataId=dataId.copy(),
                                                  associated=associatedUnitEntries)
                    # Add the entry into the Registry.
                    self.butler.registry.addDataUnitEntry(unit.name, unitEntryDict)
                # Add the entry into the cache.
                self.unitEntryCache[unit.name][unitPrimaryKeyTuple] = unitEntryDict
            associatedUnitEntries[unit.name] = unitEntryDict

        return headers, dataId

    def ingestFile(self, file, headers, dataId, run=None):
        """Ingest a single raw file into the repository.

        All necessary DataUnit entres must already be present.

        This method is not transactional; it must be wrapped in a
        ``with self.butler.transaction` block to make per-file ingest
        atomic.

        Parameters
        ----------
        file : `str` or path-like object
            Absolute path to the file to be ingested.
        headers : `list` of `~lsst.daf.base.PropertyList`
            Result of calling `readHeaders`.
        dataId : `dict`
            Data ID dictionary, as returned by `extractDataId`.
        run : `~lsst.daf.butler.Run`, optional
            Run to add the Dataset to; defaults to ``self.butler.run``.
        """
        if run is None:
            run = self.butler.run

        # Add a Dataset entry to the Registry.
        try:
            # We use transactional=False here (a kwarg added by the
            # @transactional decorator) to keep the conflict exception from
            # starting a higher-level rollback - if we catch this exception,
            # we don't want to have already started rolling back the ingest of
            # *previous* files when config.onError=='rollback' but
            # config.confict=='ignore'.
            ref = self.butler.registry.addDataset(self.datasetType, dataId, run=run,
                                                  transactional=False, recursive=True)
        except ValueError:
            raise IngestConflictError("Ingest conflict on {} {}".format(file, dataId))

        # Ingest it into the Datastore.
        self.butler.datastore.ingest(file, ref, formatter=self.getFormatter(file, headers, dataId),
                                     transfer=self.config.transfer)
        return None

    def processFile(self, file):
        """Ingest a single raw data file after extacting metadata.

        This creates any new Exposure or Visit DataUnit entries needed to
        identify the ingest file, creates a new Dataset entry in the
        Registry and finally ingests the file itself into the Datastore.
        Any needed Camera, Sensor, and PhysicalFilter DataUnit entries must
        exist in the Registry before `run` is called.

        Parameters
        ----------
        file : `str` or path-like object
            Absolute path to the file to be ingested.
        """
        headers, dataId = self.ensureDataUnits(file)
        # We want ingesting a single file to be atomic even if we are
        # not trying to ingest the list of files atomically.
        with self.butler.transaction():
            try:
                self.ingestFile(file, headers, dataId)
                return
            except IngestConflictError:
                if self.config.conflict == "fail":
                    raise
        if self.config.conflict == "ignore":
            if self.stashRun is not None:
                if self.stashRun.id is None:
                    self.butler.registry.ensureRun(self.stashRun)
                self.log.infof("Conflict on {} ({}); ingesting to stash '{}' instead.",
                               dataId, file, self.config.stash)
                with self.butler.transaction():
                    self.ingestFile(file, headers, dataId, run=self.stashRun)
            else:
                self.log.infof("Conflict on {} ({}); ignoring.", dataId, file)

    @abstractmethod
    def extractDataId(self, file, headers):
        """Return the Data ID dictionary that should be used to label a file.

        Parameters
        ----------
        file : `str` or path-like object
            Absolute path to the file being ingested (prior to any transfers).
        headers : `list` of `~lsst.daf.base.PropertyList`
            All headers returned by `readHeaders()`.

        Returns
        -------
        dataId : `dict`
            Must include "camera", "sensor", and "exposure" keys. If the
            Exposure is associated with a PhysicalFilter and/or Visit,
            "physical_filter" and "visit" keys should be provided as well
            (respectively).
        """
        raise NotImplementedError("Must be implemented by subclasses.")

    @abstractmethod
    def extractVisitEntry(self, file, headers, dataId, associated):
        """Create a Visit DataUnit entry from raw file metadata.

        Parameters
        ----------
        file : `str` or path-like object
            Absolute path to the file being ingested (prior to any transfers).
        headers : `list` of `~lsst.daf.base.PropertyList`
            All headers returned by `readHeaders()`.
        dataId : `dict`
            The data ID for this file.  Implementations are permitted to
            modify this dictionary (generally by stripping off "sensor" and
            "exposure" and adding new metadata key-value pairs) and return it.
        associated : `dict`
            A dictionary containing other associated DataUnit entries.
            Guaranteed to have "Camera", "Sensor",  and "PhysicalFilter" keys,
            but the last may map to ``None`` if `extractDataId` either did not
            contain a "physical_filter" key or mapped it to ``None``.
            Subclasses may add new keys to this dict to pass arbitrary data to
            `extractExposureEntry` (`extractVisitEntry` is always called
            first), but note that when a Visit is comprised of multiple
            Exposures, `extractVisitEntry` may not be called at all.

        Returns
        -------
        entry : `dict`
            Dictionary corresponding to an Visit database table row.
            Must have all non-null columns in the Visit table as keys.
        """
        raise NotImplementedError("Must be implemented by subclasses.")

    @abstractmethod
    def extractExposureEntry(self, file, headers, dataId, associated):
        """Create an Exposure DataUnit entry from raw file metadata.

        Parameters
        ----------
        file : `str` or path-like object
            Absolute path to the file being ingested (prior to any transfers).
        headers : `list` of `~lsst.daf.base.PropertyList`
            All headers returned by `readHeaders()`.
        dataId : `dict`
            The data ID for this file.  Implementations are permitted to
            modify this dictionary (generally by stripping off "sensor" and
            adding new metadata key-value pairs) and return it.
        associated : `dict`
            A dictionary containing other associated DataUnit entries.
            Guaranteed to have "Camera", "Sensor", "PhysicalFilter", and
            "Visit" keys, but the latter two may map to ``None`` if
            `extractDataId` did not contain keys for these or mapped them to
            ``None``.  May also contain additional keys added by
            `extractVisitEntry`.

        Returns
        -------
        entry : `dict`
            Dictionary corresponding to an Exposure database table row.
            Must have all non-null columns in the Exposure table as keys.
        """
        raise NotImplementedError("Must be implemented by subclasses.")

    def getFormatter(self, file, headers, dataId):
        """Return the Formatter that should be used to read this file after
        ingestion.

        The default implementation returns None, which uses the formatter
        configured for this DatasetType/StorageClass in the Butler.
        """
        return None


class VisitInfoRawIngestTask(RawIngestTask):
    """An intermediate base class of RawIngestTask for cameras that already
    implement constructing a `afw.image.VisitInfo` object from raw data.

    Subclasses must provide (at least) implementations of `extractDataId` and
    the new `makeVisitInfo` method; the latter is used to provide concrete
    implementations of `extractVisitEntry` and `extractExposureEntry`.
    """

    @abstractmethod
    def makeVisitInfo(self, headers, exposureId):
        """Return an `afw.image.VisitInfo` object from the given header and ID.

        Parameters
        ----------
        headers : `list` of `~lsst.daf.base.PropertyList`
            All headers returned by `readHeaders()`.
        exposureId : `int`
            Integer ID to pass to the `VisitInfo` constructor.
        """
        raise NotImplementedError("Must be implemented by subclasses.")

    def extractVisitEntry(self, file, headers, dataId, associated):
        """Create a Visit DataUnit entry from raw file metadata.

        Parameters
        ----------
        file : `str` or path-like object
            Absolute path to the file being ingested (prior to any transfers).
        headers : `list` of `~lsst.daf.base.PropertyList`
            All headers returned by `readHeaders()`.
        dataId : `dict`
            The data ID for this file.  Implementations are permitted to
            modify this dictionary (generally by stripping off "sensor" and
            "exposure" and adding new metadata key-value pairs) and return it.
        associated : `dict`
            A dictionary containing other associated DataUnit entries.
            Guaranteed to have "Camera", "Sensor",  and "PhysicalFilter" keys,
            but the last may map to ``None`` if `extractDataId` either did not
            contain a "physical_filter" key or mapped it to ``None``.
            Also adds a "VisitInfo" key containing an `afw.image.VisitInfo`
            object for use by `extractExposureEntry`.

        Returns
        -------
        entry : `dict`
            Dictionary corresponding to an Visit database table row.
            Must have all non-null columns in the Visit table as keys.
        """
        visitInfo = self.makeVisitInfo(headers, exposureId=dataId["exposure"])
        associated["VisitInfo"] = visitInfo
        del dataId["sensor"]
        del dataId["exposure"]
        return makeVisitEntryFromVisitInfo(dataId, visitInfo)

    def extractExposureEntry(self, file, headers, dataId, associated):
        """Create an Exposure DataUnit entry from raw file metadata.

        Parameters
        ----------
        file : `str` or path-like object
            Absolute path to the file being ingested (prior to any transfers).
        headers : `list` of `~lsst.daf.base.PropertyList`
            All headers returned by `readHeaders()`.
        dataId : `dict`
            The data ID for this file.  Implementations are permitted to
            modify this dictionary (generally by stripping off "sensor" and
            adding new metadata key-value pairs) and return it.
        associated : `dict`
            A dictionary containing other associated DataUnit entries.
            Guaranteed to have "Camera", "Sensor", "PhysicalFilter", and
            "Visit" keys, but the latter two may map to ``None`` if
            `extractDataId` did not contain keys for these or mapped them to
            ``None``.  May also contain additional keys added by
            `extractVisitEntry`.

        Returns
        -------
        entry : `dict`
            Dictionary corresponding to an Exposure database table row.
            Must have all non-null columns in the Exposure table as keys.
        """
        try:
            visitInfo = associated["VisitInfo"]
        except KeyError:
            visitInfo = self.makeVisitInfo(headers, exposureId=dataId["exposure"])
        del dataId["sensor"]
        return makeExposureEntryFromVisitInfo(dataId, visitInfo)
