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

__all__ = ("Translator", "TranslatorFactory", "KeyHandler", "CopyKeyHandler", "ConstantKeyHandler",
           "CalibKeyHandler", "BandToPhysicalFilterKeyHandler", "PhysicalFilterToBandKeyHandler",
           "makeCalibrationLabel")

import itertools
from typing import Optional, Any, Dict, Tuple, FrozenSet, Iterable, List
from abc import ABCMeta, abstractmethod

from lsst.log import Log
from lsst.skymap import BaseSkyMap


def makeCalibrationLabel(datasetTypeName: str, calibDate: str, ccd: Optional[int] = None,
                         filter: Optional[str] = None) -> str:
    """Make a Gen3 calibration_label string corresponding to a Gen2 data ID.

    Parameters
    ----------
    datasetTypeName : `str`
        Name of the dataset type this calibration label identifies.
    calibDate : `str`
        Date string used in the Gen2 template.
    ccd : `int`, optional
        Detector ID used in the Gen2 template.
    filter : `str`, optional
        Filter used in the Gen2 template.

    Returns
    -------
    label : `str`
        Calibration label string.
    """
    # TODO: this function is probably HSC-specific, but I don't know how other
    # obs calib registries behave so I don't know (yet) how to generalize it.
    elements = [datasetTypeName, calibDate]
    if ccd is not None:
        elements.append(f"{ccd:03d}")
    if filter is not None:
        elements.append(filter)
    return "gen2/{}".format("_".join(elements))


class KeyHandler(metaclass=ABCMeta):
    """Base class for Translator helpers that each handle just one Gen3 Data
    ID key.

    Parameters
    ----------
    dimension : `str`
        Name of the Gen3 dimension (data ID key) populated by
        this handler (e.g. "visit" or "band").
    """
    def __init__(self, dimension: str):
        self.dimension = dimension

    __slots__ = ("dimension",)

    def __repr__(self):
        return f"{type(self).__name__}({self.dimension}, ...)"

    def translate(self, gen2id: dict, gen3id: dict,
                  skyMap: Optional[BaseSkyMap], skyMapName: Optional[str],
                  datasetTypeName: str):
        """Update a Gen3 data ID dict with a single key-value pair from a Gen2
        data ID.

        This method is implemented by the base class and is not expected to
        be re-implemented by subclasses.

        Parameters
        ----------
        gen2id: `dict`
            Gen2 data ID from which to draw key-value pairs from.
        gen3id: `dict`
            Gen3 data ID to update in-place.
        skyMap: `BaseSkyMap`, optional
            SkyMap that defines the tracts and patches used in the Gen2 data
            ID, if any.
        skyMapName: `str`
            Name of the Gen3 skymap dimension that defines the tracts and
            patches used in the Gen3 data ID.
        datasetTypeName: `str`
            Name of the dataset type.
        """
        gen3id[self.dimension] = self.extract(gen2id, skyMap=skyMap, skyMapName=skyMapName,
                                              datasetTypeName=datasetTypeName)

    @abstractmethod
    def extract(self, gen2id: dict, skyMap: Optional[BaseSkyMap], skyMapName: Optional[str],
                datasetTypeName: str) -> Any:
        """Extract a Gen3 data ID value from a Gen2 data ID.

        Parameters
        ----------
        gen2id: `dict`
            Gen2 data ID from which to draw key-value pairs from.
        skyMap: `BaseSkyMap`, optional
            SkyMap that defines the tracts and patches used in the Gen2 data
            ID, if any.
        skyMapName: `str`
            Name of the Gen3 skymap dimension that defines the tracts and
            patches used in the Gen3 data ID.
        datasetTypeName: `str`
            Name of the dataset type.
        """
        raise NotImplementedError()


class ConstantKeyHandler(KeyHandler):
    """A KeyHandler that adds a constant key-value pair to the Gen3 data ID.

    Parameters
    ----------
    dimension : `str`
        Name of the Gen3 dimension (data ID key) populated by
        this handler (e.g. "visit" or "band").
    value : `object`
        Data ID value.
    """
    def __init__(self, dimension: str, value: Any):
        super().__init__(dimension)
        self.value = value

    __slots__ = ("value",)

    def extract(self, gen2id: dict, skyMap: Optional[BaseSkyMap], skyMapName: Optional[str],
                datasetTypeName: str) -> Any:
        # Docstring inherited from KeyHandler.extract.
        return self.value


class CopyKeyHandler(KeyHandler):
    """A KeyHandler that simply copies a value from a Gen3 data ID.

    Parameters
    ----------
    dimension : `str`
        Name of the Gen3 dimension produced by this handler.
    dtype : `type`, optional
        If not `None`, the type that values for this key must be an
        instance of.
    """
    def __init__(self, dimension: str, gen2key: Optional[str] = None,
                 dtype: Optional[type] = None):
        super().__init__(dimension)
        self.gen2key = gen2key if gen2key is not None else dimension
        self.dtype = dtype

    __slots__ = ("gen2key", "dtype")

    def __str__(self):
        return f"{type(self).__name__}({self.gen2key}, {self.dtype})"

    def extract(self, gen2id: dict, skyMap: Optional[BaseSkyMap], skyMapName: Optional[str],
                datasetTypeName: str) -> Any:
        # Docstring inherited from KeyHandler.extract.
        r = gen2id[self.gen2key]
        if self.dtype is not None:
            try:
                r = self.dtype(r)
            except ValueError as err:
                raise TypeError(
                    f"'{r}' is not a valid value for {self.dimension}; "
                    f"expected {self.dtype.__name__}, got {type(r).__name__}."
                ) from err
        return r


class PatchKeyHandler(KeyHandler):
    """A KeyHandler for skymap patches.
    """
    def __init__(self):
        super().__init__("patch")

    __slots__ = ()

    def extract(self, gen2id: dict, skyMap: Optional[BaseSkyMap], skyMapName: Optional[str],
                datasetTypeName: str) -> Any:
        # Docstring inherited from KeyHandler.extract.
        tract = gen2id["tract"]
        tractInfo = skyMap[tract]
        x, y = gen2id["patch"].split(",")
        patchInfo = tractInfo[int(x), int(y)]
        return tractInfo.getSequentialPatchIndex(patchInfo)


class SkyMapKeyHandler(KeyHandler):
    """A KeyHandler for skymaps."""
    def __init__(self):
        super().__init__("skymap")

    __slots__ = ()

    def extract(self, gen2id: dict, skyMap: Optional[BaseSkyMap], skyMapName: Optional[str],
                datasetTypeName: str) -> Any:
        # Docstring inherited from KeyHandler.extract.
        return skyMapName


class CalibKeyHandler(KeyHandler):
    """A KeyHandler for master calibration datasets.
    """
    __slots__ = ("ccdKey",)

    def __init__(self, ccdKey="ccd"):
        self.ccdKey = ccdKey
        super().__init__("calibration_label")

    def extract(self, gen2id: dict, skyMap: Optional[BaseSkyMap], skyMapName: Optional[str],
                datasetTypeName: str) -> Any:
        # Docstring inherited from KeyHandler.extract.
        return makeCalibrationLabel(datasetTypeName, gen2id["calibDate"],
                                    ccd=gen2id.get(self.ccdKey), filter=gen2id.get("filter"))


class PhysicalFilterToBandKeyHandler(KeyHandler):
    """KeyHandler for gen2 ``filter`` keys that match ``physical_filter``
    keys in gen3 but should be mapped to ``band``.

    Note that multiple physical filter can potentially map to one abstract
    filter, so be careful to only use this translator on obs packages where
    there is a one-to-one mapping.
    """

    __slots__ = ("_map",)

    def __init__(self, filterDefinitions):
        super().__init__("band")
        self._map = {d.physical_filter: d.band for d in filterDefinitions
                     if d.physical_filter is not None}

    def extract(self, gen2id, *args, **kwargs):
        physical = gen2id["filter"]
        return self._map.get(physical, physical)


class BandToPhysicalFilterKeyHandler(KeyHandler):
    """KeyHandler for gen2 ``filter`` keys that match ``band``
    keys in gen3 but should be mapped to ``physical_filter``.

    Note that one abstract filter can potentially map to multiple physical
    filters, so be careful to only use this translator on obs packages where
    there is a one-to-one mapping.
    """

    __slots__ = ("_map",)

    def __init__(self, filterDefinitions):
        super().__init__("physical_filter")
        self._map = {d.band: d.physical_filter for d in filterDefinitions
                     if d.band is not None}

    def extract(self, gen2id, *args, **kwargs):
        abstract = gen2id["filter"]
        return self._map.get(abstract, abstract)


class TranslatorFactory:
    """A class that manages a system of rules for translating Gen2 data IDs
    to Gen3 data IDs, and uses these to construct translators for particular
    dataset types.

    Parameters
    ----------
    log : `lsst.log.Log`, optional
        A logger for diagnostic messages.
    """
    def __init__(self, log: Optional[Log] = None):
        # The rules used to match KeyHandlers when constructing a Translator.
        self._rules: Dict[
            Optional[str],  # instrument name (or None to match any)
            Dict[
                Optional[str],  # dataset type name (or None to match any)
                # gen2keys, handler, consume
                List[Tuple[FrozenSet[str], KeyHandler, bool]]
            ]
        ] = {
            None: {
                None: []
            }
        }
        self._addDefaultRules()
        if log is None:
            log = Log.getLogger("obs.base.gen2to3.TranslatorFactory")
        self.log = log

    def __str__(self):
        lines = []
        for instrumentName, nested in self._rules.items():
            if instrumentName is None:
                instrumentName = "[any instrument]"
            for datasetTypeName, rules in nested.items():
                if datasetTypeName is None:
                    datasetTypeName = "[any dataset type]"
                lines.append(f"{instrumentName} + {datasetTypeName}:")
                for gen2keys, handler, consume in rules:
                    consumed = " (consumed)" if consume else ""
                    lines.append(f"   {gen2keys}{consumed}: {handler}")
        return "\n".join(lines)

    def addRule(self, handler: KeyHandler, instrument: Optional[str] = None,
                datasetTypeName: Optional[str] = None, gen2keys: Iterable[str] = (),
                consume: bool = True):
        """Add a KeyHandler and an associated matching rule.

        Parameters
        ----------
        handler : `KeyHandler`
            A KeyHandler instance to add to a Translator when this rule
            matches.
        instrument : `str`
            Gen3 instrument name the Gen2 repository must be associated with
            for this rule to match, or None to match any instrument.
        datasetTypeName : `str`
            Name of the DatasetType this rule matches, or None to match any
            DatasetType.
        gen2Keys : sequence
            Sequence of Gen2 data ID keys that must all be present for this
            rule to match.
        consume : `bool` or `tuple`
            If True (default), remove all entries in gen2keys from the set of
            keys being matched to in order to prevent less-specific handlers
            from matching them.
            May also be a `tuple` listing only the keys to consume.
        """
        # Ensure consume is always a frozenset, so we can process it uniformly
        # from here on.
        if consume is True:
            consume = frozenset(gen2keys)
        elif consume:
            consume = frozenset(consume)
        else:
            consume = frozenset()
        # find the rules for this instrument, or if we haven't seen it before,
        # add a nested dictionary that matches any DatasetType name and then
        # append this rule.
        rulesForInstrument = self._rules.setdefault(instrument, {None: []})
        rulesForInstrumentAndDatasetType = rulesForInstrument.setdefault(datasetTypeName, [])
        rulesForInstrumentAndDatasetType.append((frozenset(gen2keys), handler, consume))

    def _addDefaultRules(self):
        """Add translator rules that should always be present, and don't depend
        at all on the instrument whose datasets are being converted.

        This is called by `TranslatorFactory` construction.
        """
        # Add "skymap" to Gen3 ID if Gen2 ID has a "tract" key.
        self.addRule(SkyMapKeyHandler(), gen2keys=("tract",), consume=False)

        # Add "skymap" to Gen3 ID if DatasetType is one of a few specific ones
        for coaddName in ("deep", "goodSeeing", "psfMatched", "dcr"):
            self.addRule(SkyMapKeyHandler(), datasetTypeName=f"{coaddName}Coadd_skyMap")

        # Translate Gen2 str patch IDs to Gen3 sequential integers.
        self.addRule(PatchKeyHandler(), gen2keys=("patch",))

        # Translate any "filter" values that appear alongside "tract" to
        # "band".  This is _not_ the right choice for instruments
        # that use "physical_filter" values for coadds in Gen2 (like HSC);
        # those will need to add a rule that invokes
        # PhysicalFilterToBandKey instead for just that instrument, but the
        # same criteria otherwise.  That will override this one, because
        # instrument-specific rules match first, and that rule will consume
        # the Gen2 "filter" key before this rule has a chance to fire.
        self.addRule(CopyKeyHandler("band", "filter"),
                     gen2keys=("filter", "tract"),
                     consume=("filter",))

        # Copy Gen2 "tract" to Gen3 "tract".
        self.addRule(CopyKeyHandler("tract", dtype=int), gen2keys=("tract",))

        # Add valid_first, valid_last to instrument-level transmission/ datasets;
        # these are considered calibration products in Gen3.
        for datasetTypeName in ("transmission_sensor", "transmission_optics", "transmission_filter"):
            self.addRule(ConstantKeyHandler("calibration_label", "unbounded"),
                         datasetTypeName=datasetTypeName)

        # Translate Gen2 pixel_id to Gen3 skypix.
        #
        # TODO: For now, we just assume that the refcat indexer uses htm7,
        # since that's what we have generated most of our refcats at.
        # Eventually that may have to change, but it's not clear enough how to
        # do that for us to have a ticket yet.  If you found this note because
        # you've run into this limitation, please let the middleware team know
        # that it's time to make this a priority.
        self.addRule(CopyKeyHandler("htm7", gen2key="pixel_id", dtype=int), gen2keys=("pixel_id",))

    def addGenericInstrumentRules(self, instrumentName: str,
                                  calibFilterType: str = "physical_filter",
                                  detectorKey: str = "ccd",
                                  exposureKey: str = "visit"):
        """Add translation rules that depend on some properties of the
        instrument but are otherwise generic.

        Parameters
        ----------
        instrument : `str`
            The short (dimension) name of the instrument that conversion is
            going to be run on.
        calibFilterType : `str`, optional
            One of ``physical_filter`` or ``band``, indicating which
            of those the gen2 calibRegistry uses as the ``filter`` key.
        detectorKey : `str`, optional
            The gen2 key used to identify what in gen3 is `detector`.
        exposureKey : `str`, optional
            The gen2 key used to identify what in gen3 is `exposure`.
        """
        # Add instrument to Gen3 data ID if Gen2 contains exposureKey or
        # detectorKey.  (Both rules will match, so we'll actually set
        # instrument in the same dict twice).
        self.addRule(ConstantKeyHandler("instrument", instrumentName),
                     instrument=instrumentName, gen2keys=(exposureKey,), consume=False)
        self.addRule(ConstantKeyHandler("instrument", instrumentName),
                     instrument=instrumentName, gen2keys=(detectorKey,), consume=False)
        self.addRule(ConstantKeyHandler("instrument", instrumentName),
                     instrument=instrumentName, gen2keys=("calibDate",), consume=False)

        # Copy Gen2 exposureKey to Gen3 'exposure' for raw only.  Also consume
        # filter, since that's implied by 'exposure' in Gen3.
        self.addRule(CopyKeyHandler("exposure", exposureKey),
                     instrument=instrumentName, datasetTypeName="raw", gen2keys=(exposureKey,),
                     consume=(exposureKey, "filter"))

        # Copy Gen2 'visit' to Gen3 'visit' otherwise.  Also consume filter.
        self.addRule(CopyKeyHandler("visit"), instrument=instrumentName, gen2keys=("visit",),
                     consume=("visit", "filter"))

        # Copy Gen2 'ccd' to Gen3 'detector;
        self.addRule(CopyKeyHandler("detector", detectorKey),
                     instrument=instrumentName,
                     gen2keys=(detectorKey,))

        # Add instrument for transmission curve datasets (transmission_sensor is
        # already handled by the above rules).
        self.addRule(ConstantKeyHandler("instrument", instrumentName),
                     instrument=instrumentName, datasetTypeName="transmission_optics")
        self.addRule(ConstantKeyHandler("instrument", instrumentName),
                     instrument=instrumentName, datasetTypeName="transmission_atmosphere")
        self.addRule(ConstantKeyHandler("instrument", instrumentName),
                     instrument=instrumentName, datasetTypeName="transmission_filter")
        self.addRule(CopyKeyHandler("physical_filter", "filter"),
                     instrument=instrumentName, datasetTypeName="transmission_filter")

        # Add calibration mapping for filter dependent types
        for calibType in ('flat', 'sky', 'fringe'):
            self.addRule(CopyKeyHandler(calibFilterType, "filter"),
                         instrument=instrumentName, datasetTypeName=calibType)

        # Translate Gen2 calibDate and datasetType to Gen3 calibration_label.
        self.addRule(CalibKeyHandler(detectorKey), gen2keys=("calibDate",))

    def makeMatching(self, datasetTypeName: str, gen2keys: Dict[str, type], instrument: Optional[str] = None,
                     skyMap: Optional[BaseSkyMap] = None, skyMapName: Optional[str] = None):
        """Construct a Translator appropriate for instances of the given
        dataset.

        Parameters
        ----------
        datasetTypeName : `str`
            Name of the dataset type.
        gen2keys: `dict`
            Keys of a Gen2 data ID for this dataset.
        instrument: `str`, optional
            Name of the Gen3 instrument dimension for translated data IDs.
        skyMap: `~lsst.skymap.BaseSkyMap`, optional
            The skymap instance that defines any tract/patch data IDs.
            `~lsst.skymap.BaseSkyMap` instances.
        skyMapName : `str`, optional
            Gen3 SkyMap Dimension name to be associated with any tract or patch
            Dimensions.

        Returns
        -------
        translator : `Translator`
            A translator whose translate() method can be used to transform Gen2
            data IDs to Gen3 dataIds.
        """
        if instrument is not None:
            rulesForInstrument = self._rules.get(instrument, {None: []})
        else:
            rulesForInstrument = {None: []}
        rulesForAnyInstrument = self._rules[None]
        candidateRules = itertools.chain(
            rulesForInstrument.get(datasetTypeName, []),     # this instrument, this DatasetType
            rulesForInstrument[None],                         # this instrument, any DatasetType
            rulesForAnyInstrument.get(datasetTypeName, []),  # any instrument, this DatasetType
            rulesForAnyInstrument[None],                      # any instrument, any DatasetType
        )
        matchedHandlers = []
        targetKeys = set(gen2keys)
        self.log.debug("Constructing data ID translator for %s with Gen2 keys %s...",
                       datasetTypeName, gen2keys)
        for ruleKeys, ruleHandlers, consume in candidateRules:
            if ruleKeys.issubset(targetKeys):
                matchedHandlers.append(ruleHandlers)
                targetKeys -= consume
        self.log.debug("...matched %d handlers: %s, with %s unmatched.",
                       len(matchedHandlers), matchedHandlers, targetKeys)
        return Translator(matchedHandlers, skyMap=skyMap, skyMapName=skyMapName,
                          datasetTypeName=datasetTypeName, log=self.log)


class Translator:
    """Callable object that translates Gen2 Data IDs to Gen3 Data IDs for a
    particular DatasetType.

    Translators should usually be constructed via
    `TranslatorFactory.makeMatching`.

    Parameters
    ----------
    handlers : `list`
        A list of KeyHandlers this Translator should use.
    skyMap : `BaseSkyMap`, optional
        SkyMap instance used to define any tract or patch Dimensions.
    skyMapName : `str`
        Gen3 SkyMap Dimension name to be associated with any tract or patch
        Dimensions.
    datasetTypeName : `str`
        Name of the dataset type whose data IDs this translator handles.
    """
    def __init__(self, handlers: List[KeyHandler], skyMap: Optional[BaseSkyMap], skyMapName: Optional[str],
                 datasetTypeName: str, log: Log):
        self.handlers = handlers
        self.skyMap = skyMap
        self.skyMapName = skyMapName
        self.datasetTypeName = datasetTypeName
        self.log = log

    __slots__ = ("handlers", "skyMap", "skyMapName", "datasetTypeName", "log")

    def __str__(self):
        hstr = ",".join(str(h) for h in self.handlers)
        return f"{type(self).__name__}(dtype={self.datasetTypeName}, handlers=[{hstr}])"

    def __call__(self, gen2id: Dict[str, Any], *, partial: bool = False):
        """Return a Gen3 data ID that corresponds to the given Gen2 data ID.
        """
        gen3id = {}
        for handler in self.handlers:
            try:
                handler.translate(gen2id, gen3id, skyMap=self.skyMap, skyMapName=self.skyMapName,
                                  datasetTypeName=self.datasetTypeName)
            except KeyError:
                if partial:
                    self.log.debug("Failed to translate %s from %s (this may not be an error).",
                                   handler.dimension, gen2id)
                    continue
                else:
                    raise
        return gen3id

    @property
    def dimensionNames(self):
        """The names of the dimensions populated by this Translator
        (`frozenset`).
        """
        return frozenset(h.dimension for h in self.handlers)
