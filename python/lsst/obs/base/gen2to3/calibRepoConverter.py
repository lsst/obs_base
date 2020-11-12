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

__all__ = ["CalibRepoConverter"]

from collections import defaultdict
import os
import sqlite3
from typing import TYPE_CHECKING, Dict, Iterator, List, Mapping, Sequence, Tuple, Optional

import astropy.time
import astropy.units as u

from lsst.daf.butler import DataCoordinate, FileDataset, Timespan
from .repoConverter import RepoConverter
from .repoWalker import RepoWalker

if TYPE_CHECKING:
    from lsst.daf.butler import DatasetType, StorageClass, FormatterParameter
    from .repoWalker.scanner import PathElementHandler
    from ..cameraMapper import CameraMapper
    from ..mapping import Mapping as CameraMapperMapping  # disambiguate from collections.abc.Mapping


class CalibRepoConverter(RepoConverter):
    """A specialization of `RepoConverter` for calibration repositories.

    Parameters
    ----------
    mapper : `CameraMapper`
        Gen2 mapper for the data repository.  The root associated with the
        mapper is ignored and need not match the root of the repository.
    labels : `Sequence` [ `str` ]
        Strings injected into the names of the collections that calibration
        datasets are written and certified into (forwarded as the ``extra``
        argument to `Instrument` methods that generate collection names and
        write curated calibrations).
    **kwargs
        Additional keyword arguments are forwarded to (and required by)
        `RepoConverter`.
    """

    def __init__(self, *, mapper: CameraMapper, labels: Sequence[str] = (), **kwargs):
        super().__init__(run=None, **kwargs)
        self.mapper = mapper
        self.collection = self.task.instrument.makeCalibrationCollectionName(*labels)
        self._labels = tuple(labels)
        self._datasetTypes = set()

    def isDatasetTypeSpecial(self, datasetTypeName: str) -> bool:
        # Docstring inherited from RepoConverter.
        return datasetTypeName in self.instrument.getCuratedCalibrationNames()

    def iterMappings(self) -> Iterator[Tuple[str, CameraMapperMapping]]:
        # Docstring inherited from RepoConverter.
        yield from self.mapper.calibrations.items()

    def makeRepoWalkerTarget(self, datasetTypeName: str, template: str, keys: Dict[str, type],
                             storageClass: StorageClass, formatter: FormatterParameter = None,
                             targetHandler: Optional[PathElementHandler] = None,
                             ) -> RepoWalker.Target:
        # Docstring inherited from RepoConverter.
        target = RepoWalker.Target(
            datasetTypeName=datasetTypeName,
            storageClass=storageClass,
            template=template,
            keys=keys,
            instrument=self.task.instrument.getName(),
            universe=self.task.registry.dimensions,
            formatter=formatter,
            targetHandler=targetHandler,
            translatorFactory=self.task.translatorFactory,
        )
        self._datasetTypes.add(target.datasetType)
        return target

    def _queryGen2CalibRegistry(self, db: sqlite3.Connection, datasetType: DatasetType, calibDate: str
                                ) -> Iterator[sqlite3.Row]:
        """Query the Gen2 calibration registry for the validity ranges and
        optionally detectors and filters associated with the given dataset type
        and ``calibDate``.

        Parameters
        ----------
        db : `sqlite3.Connection`
            DBAPI connection to the Gen2 ``calibRegistry.sqlite3`` file.
        datasetType : `DatasetType`
            Gen3 dataset type being queried.
        calibDate : `str`
            String extracted from the ``calibDate`` template entry in Gen2
            filenames.

        Yields
        ------
        row : `sqlite3.Row`
            SQLite result object; will have ``validStart`` and ``validEnd``
            columns, may have a detector column (named
            ``self.task.config.ccdKey``) and/or a ``filter`` column, depending
            on whether ``datasetType.dimensions`` includes ``detector`` and
            ``physical_filter``, respectively.
        """
        fields = ["validStart", "validEnd"]
        if "detector" in datasetType.dimensions.names:
            fields.append(self.task.config.ccdKey)
        else:
            fields.append(f"NULL AS {self.task.config.ccdKey}")
        if "physical_filter" in datasetType.dimensions.names:
            fields.append("filter")
        else:
            assert "band" not in datasetType.dimensions.names
            fields.append("NULL AS filter")
        tables = self.mapper.mappings[datasetType.name].tables
        if tables is None or len(tables) == 0:
            self.task.log.warn("Could not extract calibration ranges for %s in %s; "
                               "no tables in Gen2 mapper.",
                               datasetType.name, self.root, tables[0])
            return
        query = f"SELECT DISTINCT {', '.join(fields)} FROM {tables[0]} WHERE calibDate = ?;"
        try:
            results = db.execute(query, (calibDate,))
        except sqlite3.OperationalError as e:
            self.task.log.warn("Could not extract calibration ranges for %s in %s from table %s: %r",
                               datasetType.name, self.root, tables[0], e)
            return
        yield from results

    def _finish(self, datasets: Mapping[DatasetType, Mapping[Optional[str], List[FileDataset]]]):
        # Docstring inherited from RepoConverter.
        # Read Gen2 calibration repository and extract validity ranges for
        # all datasetType + calibDate combinations we ingested.
        calibFile = os.path.join(self.root, "calibRegistry.sqlite3")
        # If the registry file does not exist this indicates a problem.
        # We check explicitly because sqlite will try to create the
        # missing file if it can.
        if not os.path.exists(calibFile):
            raise RuntimeError("Attempting to convert calibrations but no registry database"
                               f" found in {self.root}")

        # Initially we collate timespans for each dataId + dataset type
        # combination. This allows us to check for small gaps or overlaps
        # inherent in the ambiguous usage of validity ranges in gen2
        timespansByDataId = defaultdict(list)

        db = sqlite3.connect(calibFile)
        db.row_factory = sqlite3.Row

        for datasetType, datasetsByCalibDate in datasets.items():
            if not datasetType.isCalibration():
                continue
            gen2keys = {}
            if "detector" in datasetType.dimensions.names:
                gen2keys[self.task.config.ccdKey] = int
            if "physical_filter" in datasetType.dimensions.names:
                gen2keys["filter"] = str
            translator = self.instrument.makeDataIdTranslatorFactory().makeMatching(
                datasetType.name,
                gen2keys,
                instrument=self.instrument.getName()
            )
            for calibDate, datasetsForCalibDate in datasetsByCalibDate.items():
                assert calibDate is not None, ("datasetType.isCalibration() is set by "
                                               "the presence of calibDate in the Gen2 template")
                # Build a mapping that lets us find DatasetRefs by data ID,
                # for this DatasetType and calibDate.  We know there is only
                # one ref for each data ID (given DatasetType and calibDate as
                # well).
                refsByDataId = {}
                for dataset in datasetsForCalibDate:
                    refsByDataId.update((ref.dataId, ref) for ref in dataset.refs)
                # Query the Gen2 calibration repo for the validity ranges for
                # this DatasetType and calibDate, and look up the appropriate
                # refs by data ID.
                for row in self._queryGen2CalibRegistry(db, datasetType, calibDate):
                    # For validity times we use TAI as some gen2 repos have
                    # validity dates very far in the past or future.
                    timespan = Timespan(
                        astropy.time.Time(row["validStart"], format="iso", scale="tai"),
                        astropy.time.Time(row["validEnd"], format="iso", scale="tai"),
                    )
                    # Make a Gen2 data ID from query results.
                    gen2id = {}
                    if "detector" in datasetType.dimensions.names:
                        gen2id[self.task.config.ccdKey] = row[self.task.config.ccdKey]
                    if "physical_filter" in datasetType.dimensions.names:
                        gen2id["filter"] = row["filter"]
                    # Translate that to Gen3.
                    gen3id, _ = translator(gen2id)
                    dataId = DataCoordinate.standardize(gen3id, graph=datasetType.dimensions)
                    ref = refsByDataId.get(dataId)
                    if ref is not None:
                        # Validity ranges must not overlap for the same dataID
                        # datasetType combination. Use that as a primary
                        # key and store the timespan and ref in a tuple
                        # as the value for later timespan validation.
                        timespansByDataId[(ref.dataId, ref.datasetType.name)].append((timespan, ref))
                    else:
                        # The Gen2 calib registry mentions this dataset, but it
                        # isn't included in what we've ingested.  This might
                        # sometimes be a problem, but it should usually
                        # represent someone just trying to convert a subset of
                        # the Gen2 repo, so I don't think it's appropriate to
                        # warn or even log at info, since in that case there
                        # may be a _lot_ of these messages.
                        self.task.log.debug(
                            "Gen2 calibration registry entry has no dataset: %s for calibDate=%s, %s.",
                            datasetType.name, calibDate, dataId
                        )

        # Analyze the timespans to check for overlap problems
        # Gaps of a day should be closed since we assume differing
        # conventions in gen2 repos.

        # We need to correct any validity range issues and store the
        # results in a dict-of-lists keyed by Timespan, since
        # Registry.certify operates on one Timespan and multiple refs at a
        # time.
        refsByTimespan = defaultdict(list)

        # A day with a bit of fuzz to indicate the largest gap we will close
        max_gap = astropy.time.TimeDelta(1.001, format="jd", scale="tai")

        # Since in many cases the validity ranges are relevant for multiple
        # dataset types and dataIds we don't want to over-report and so
        # cache the messages for later.
        info_messages = set()
        warn_messages = set()
        for timespans in timespansByDataId.values():
            # Sort all the timespans and check overlaps
            sorted_timespans = sorted(timespans, key=lambda x: x[0])
            timespan_prev, ref_prev = sorted_timespans.pop(0)
            for timespan, ref in sorted_timespans:
                # See if we have a suspicious gap
                delta = timespan.begin - timespan_prev.end
                abs_delta = abs(delta)
                if abs_delta > 0 and abs_delta < max_gap:
                    if delta > 0:
                        # Gap between timespans
                        msg = f"Calibration validity gap closed from {timespan_prev.end} to {timespan.begin}"
                        info_messages.add(msg)
                    else:
                        # Overlap of timespans
                        msg = f"Calibration validity overlap of {abs(delta).to(u.s)} removed for period " \
                            f"{timespan.begin} to {timespan_prev.end}"
                        warn_messages.add(msg)

                    self.task.log.debug("Correcting validity range for %s with end %s",
                                        ref_prev, timespan_prev.end)

                    # Assume this gap is down to convention in gen2.
                    # We have to adjust the previous timespan to fit
                    # since we always trust validStart.
                    timespan_prev = Timespan(begin=timespan_prev.begin,
                                             end=timespan.begin)
                # Store the previous timespan and ref since it has now
                # been verified
                refsByTimespan[timespan_prev].append(ref_prev)

                # And update the previous values for the next iteration
                timespan_prev = timespan
                ref_prev = ref

            # Store the final timespan/ref pair
            refsByTimespan[timespan_prev].append(ref_prev)

        # Issue any pending log messages we have recorded
        for msg in sorted(info_messages):
            self.task.log.info(msg)
        for msg in sorted(warn_messages):
            self.task.log.warn(msg)

        # Done reading from Gen2, time to certify into Gen3.
        for timespan, refs in refsByTimespan.items():
            self.task.registry.certify(self.collection, refs, timespan)

    def getRun(self, datasetTypeName: str, calibDate: Optional[str] = None) -> str:
        # Docstring inherited from RepoConverter.
        if calibDate is None:
            return super().getRun(datasetTypeName)
        else:
            return self.instrument.makeCalibrationCollectionName(*self._labels, "gen2",
                                                                 calibDate)

    # Class attributes that will be shadowed by public instance attributes;
    # defined here only for documentation purposes.

    mapper: CameraMapper
    """Gen2 mapper associated with this repository.
    """
