#!/usr/bin/env python

import argparse
import logging

import lsst.log
from lsst.log import Log

from lsst.daf.butler import Butler


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Ingests datasets exported by exportExternalData."
    )
    parser.add_argument(
        "root",
        help="Path to butler to ingest into (usually DATA)."
    )
    parser.add_argument(
        "filename",
        help="Path to YAML file describing external files (usually resources/external.yaml)."
    )
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

    butler = Butler(args.root, run="calib/gen2")

    butler.import_(directory="/",
                   filename=args.filename, transfer="symlink")
