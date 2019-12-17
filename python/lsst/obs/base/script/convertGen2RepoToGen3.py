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

"""Convert a gen2 butler repo to gen 3. See
`lsst.obs.base.ConvertRepoConfig` for most of the config options.
"""

import argparse
import logging

import lsst.daf.butler
import lsst.log
import lsst.utils

from lsst.obs.base.gen2to3 import (ConvertRepoTask, ConvertRepoSkyMapConfig,
                                   Translator, ConstantKeyHandler, CopyKeyHandler,
                                   CalibKeyHandler)


def build_argparser():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("instrumentClass", metavar="lsst.obs.CAMERA.INSTRUMENT",
                        help=("The full import path to the gen3 Instrument class for this camera"
                              " (e.g. lsst.obs.decam.DarkEnergyCamera)."))
    parser.add_argument("--gen2root", required=True,
                        help="Root path of the gen 2 repo to be converted.")
    parser.add_argument("--gen3root", required=True,
                        help="Root path of the gen 3 repo to be produced.")
    parser.add_argument("--skymapName", default=None,
                        help="Name of the new gen3 skymap (e.g. 'discrete/ci_hsc').")
    parser.add_argument("--skymapConfig", default=None,
                        help="Path to skymap config file defining the new gen3 skymap.")
    parser.add_argument("--calibs", default=None,
                        help="Path to calibration repo; absolute, or relative to gen2root.")
    parser.add_argument("-v", "--verbose", action="store_const", dest="verbose",
                        default=lsst.log.Log.INFO, const=lsst.log.Log.DEBUG,
                        help="Set the log level to DEBUG.")
    parser.add_argument("--calibFilterType", default="physical_filter",
                        help="physical_filter or abstract_filter as the id in the gen2 calibRegistry.")
    parser.add_argument("-c", "--config", default=None,
                        help=("Path to a `ConvertRepoConfig` override to be included after "
                              "the Instrument config overrides are applied."))

    return parser


def parse_args(parser):
    args = parser.parse_args()

    skymapList = [args.skymapName, args.skymapConfig]
    if not all(x is None for x in skymapList) and not all(x is not None for x in skymapList):
        parser.error("Must specify both --skymapName and --skymapConfig, or neither.")

    return args


def configure_translators(instrument, calibFilterType, ccdKey="ccd"):
    """Configure the gen3 translators so they know the correct instrument name.

    Parameters
    ----------
    instrument : `lsst.obs.base.Instrument`
        The instrument that conversion is going to be run on.
    calibFilterType : `str`
        Whether the gen2 calibRegistry uses ``physical_filter`` or
        ``abstract_filter`` as the ``filter`` key.
    ccdKey : `str`, optional
        The gen2 key used to identify what in gen3 is `detector`.
    """
    # Add instrument to Gen3 data ID if Gen2 contains "visit" or ccdKey.
    # (Both rules will match, so we'll actually set instrument in the same dict twice).
    Translator.addRule(ConstantKeyHandler("instrument", instrument.getName()),
                       instrument=instrument.getName(), gen2keys=("visit",), consume=False)
    Translator.addRule(ConstantKeyHandler("instrument", instrument.getName()),
                       instrument=instrument.getName(), gen2keys=(ccdKey,), consume=False)

    # Copy Gen2 'visit' to Gen3 'exposure' for raw only.  Also consume filter,
    # since that's implied by 'exposure' in Gen3.
    Translator.addRule(CopyKeyHandler("exposure", "visit"),
                       instrument=instrument.getName(), datasetTypeName="raw", gen2keys=("visit",),
                       consume=("visit", "filter"))

    # Copy Gen2 'visit' to Gen3 'visit' otherwise.  Also consume filter.
    Translator.addRule(CopyKeyHandler("visit"), instrument=instrument.getName(), gen2keys=("visit",),
                       consume=("visit", "filter"))

    # Copy Gen2 'ccd' to Gen3 'detector;
    Translator.addRule(CopyKeyHandler("detector", ccdKey),
                       instrument=instrument.getName(),
                       gen2keys=(ccdKey,))

    # Add instrument for transmission curve datasets (transmission_sensor is
    # already handled by the above translators).
    Translator.addRule(ConstantKeyHandler("instrument", instrument),
                       instrument=instrument.getName(), datasetTypeName="transmission_optics")
    Translator.addRule(ConstantKeyHandler("instrument", instrument),
                       instrument=instrument.getName(), datasetTypeName="transmission_atmosphere")
    Translator.addRule(ConstantKeyHandler("instrument", instrument),
                       instrument=instrument.getName(), datasetTypeName="transmission_filter")
    Translator.addRule(CopyKeyHandler("physical_filter", "filter"),
                       instrument=instrument.getName(), datasetTypeName="transmission_filter")

    # Add calibration mapping for filter dependent types
    for calibType in ('flat', 'sky', 'fringe'):
        Translator.addRule(CopyKeyHandler(calibFilterType, "filter"),
                           instrument=instrument.getName(), datasetTypeName=calibType)

    # Translate Gen2 calibDate and datasetType to Gen3 calibration_label.
    Translator.addRule(CalibKeyHandler(ccdKey), gen2keys=("calibDate",))


def convert(gen2root, gen3root, instrumentClass, calibFilterType,
            skymapName=None, skymapConfig=None,
            calibs=None, config=None):
    """Convert the gen 2 Butler repo living at gen2root into a gen 3 repo
    living at gen3root.

    Parameters
    ----------
    gen2root : `str`
        Root path to the gen2 repo to be converted.
    gen3root : `str`
        Root path to the gen3 output repo.
    instrumentClass : `str`
        Full python path to the `lsst.obs.base.Instrument` class of the repo
        being converted.
    calibFilterType : `str`
        `abstract_filter` or `physical_filter`, depending on the type of
        ``filter`` in the gen2 calib registry.
    skymapName : `str`, optional
        Name of the skymap to be converted in the repo.
    skymapConfig : `str`, optional
        Path to the `lsst.skymap.BaseSkyMapConfig` of the gen2 skymap to be
        converted.
    calibs : `str`, optional
        Path to the gen2 calibration repository to be converted.
        If a relative path, it is assumed to be relative to ``gen2root``.
    config : `str`, optional
        Path to `lsst.obs.base.ConvertRepoConfig` configuration to load
        after all default/instrument configurations.
    """
    # instantiate the correct instrument
    instrument = lsst.utils.doImport(instrumentClass)()

    convertRepoConfig = ConvertRepoTask.ConfigClass()
    instrument.applyConfigOverrides(ConvertRepoTask._DefaultName, convertRepoConfig)
    if convertRepoConfig.raws.instrument is None:
        convertRepoConfig.raws.instrument = instrumentClass
    convertRepoConfig.raws.transfer = "symlink"
    if skymapName is not None:
        convertRepoConfig.skyMaps[skymapName] = ConvertRepoSkyMapConfig()
        convertRepoConfig.skyMaps[skymapName].load(skymapConfig)
    if config is not None:
        convertRepoConfig.load(config)

    configure_translators(instrument, calibFilterType, convertRepoConfig.ccdKey)

    butlerConfig = lsst.daf.butler.Butler.makeRepo(gen3root)
    butler = lsst.daf.butler.Butler(butlerConfig, run=instrument.getName())
    convertRepoTask = ConvertRepoTask(config=convertRepoConfig, butler3=butler)
    convertRepoTask.run(
        root=gen2root,
        # NOTE: we'd like to use `raw/NAME` and `calib/NAME` here, but if I do, I get an error about
        # `AssertionError: Multiple collections for curated calibrations is not yet supported.`
        # TODO: This is likely related to DM-23230.
        # since we also want to specify just `NAME` so that we can instantiate a gen3 Butler
        # using that as an "umbrella" collection so it knows about everything in it.
        collections=[instrument.getName()],  # ['raw/'+instrument.getName(), instrument.getName()],
        calibs=None if calibs is None else {calibs: [instrument.getName()]}
    )


def main():
    """To be run by the commandline script in `bin/`.
    """
    parser = build_argparser()
    args = parser.parse_args()

    log = lsst.log.Log.getLogger("convertRepo")
    log.setLevel(args.verbose)
    # Forward python logging to lsst logger
    logger = logging.getLogger("convertRepo")
    logger.setLevel(lsst.log.LevelTranslator.lsstLog2logging(log.getLevel()))
    logger.addHandler(lsst.log.LogHandler())

    convert(args.gen2root, args.gen3root, args.instrumentClass, args.calibFilterType,
            skymapName=args.skymapName, skymapConfig=args.skymapConfig,
            calibs=args.calibs, config=args.config)
