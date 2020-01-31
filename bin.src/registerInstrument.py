#!/usr/bin/env python

import argparse
import logging

import lsst.log
from lsst.log import Log

from lsst.daf.butler import Butler


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Add an instrument and its curated calibrations to the data repository."
    )
    parser.add_argument("root", help="Path to butler to use")
    parser.add_argument("camera", help="Camera to register.")
    parser.add_argument("-v", "--verbose", action="store_const", dest="logLevel",
                        default=Log.INFO, const=Log.DEBUG,
                        help="Set the log level to DEBUG.")

    args = parser.parse_args()
    log = Log.getLogger("lsst.daf.butler")
    log.setLevel(args.logLevel)

    # Forward python logging to lsst logger
    lgr = logging.getLogger("lsst.daf.butler")
    lgr.setLevel(logging.INFO if args.logLevel == Log.INFO else logging.DEBUG)
    lgr.addHandler(lsst.log.LogHandler())

    instrument = None
    butler = Butler(args.root, run="calib/" + args.camera)

    if args.camera in ('hsc', "HSC", 'obs_hsc', 'obs_subaru', 'subaru'):
        from lsst.obs.subaru.gen3.hsc import HyperSuprimeCam
        instrument = HyperSuprimeCam()
        args.camera = 'hsc'

        # DUPLICATED TO AVOID FILTER SHENANIGANS
        instrument.register(butler.registry)
        instrument.writeCuratedCalibrations(butler)
        # END DUPLICATION
        if args.camera == 'hsc' and False:
            instrument.ingestStrayLightData(butler,
                                            "/datasets/hsc/repo/CALIB/STRAY_LIGHT/",
                                            transfer='symlink')
    elif args.camera in ('lsst', 'obs_lsst', 'latiss', 'auxtel'):
        from lsst.obs.lsst.gen3 import LatissInstrument
        instrument = LatissInstrument()
        # DUPLICATED TO AVOID FILTER SHENANIGANS
        instrument.register(butler.registry)
        instrument.writeCuratedCalibrations(butler)
        # END DUPLICATION
