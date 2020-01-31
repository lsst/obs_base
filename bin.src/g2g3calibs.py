#!/usr/bin/env python

import argparse
import sqlite3


sourceDirs = dict()

parser = argparse.ArgumentParser(description='Way to generate HSC input yaml files.')
parser.add_argument('visit', type=int, default=19712)
parser.add_argument('--ccd', '-c', default='')
parser.add_argument('--instrument', '-i', type=str, default='HSC')
parser.add_argument('--datasetBase', '-D', type=int, default=28000)
parser.add_argument('--run', '-R', default='calib/gen2')
args = parser.parse_args()


dataSql = sqlite3.connect('/datasets/hsc/repo/registry.sqlite3')
calibSql = sqlite3.connect('/datasets/hsc/repo/CALIB/calibRegistry.sqlite3')

dataC = dataSql.cursor()
calibC = calibSql.cursor()
dataC.row_factory = sqlite3.Row
calibC.row_factory = sqlite3.Row


rawQuery = 'select DISTINCT field,taiObs,filter from raw WHERE '
rawQuery += "visit = ?"

calibQueryPre = 'select DISTINCT validStart, validEnd, ccd, filter, calibDate, id from '
calibQueryPost = ' WHERE validStart <= ? AND validEnd >= ? AND filter in (?, "NONE") '

if args.ccd != '':
    rawQuery += ' AND ccd = ?'
    calibQueryPost += ' AND ccd = ?'

yL = []
yCL = []
yCT = []
yL.append("description: Butler Data Repository Export")
yL.append("version: 0")
yL.append("data:")

yCL.append("- type: dimension")
yCL.append("  element: calibration_label")
yCL.append("  records:")

yCT = ['- type: run',
       f"  id: 99",
       f"  name: {args.run}",
       '  start_time: null',
       '  end_time: null',
       '  host: null',
       '  collection: calib/gen2',
       '  pipeline: null',
       '  environment: null']

if args.ccd != '':
    dataC.execute(rawQuery, (str(args.visit), str(args.ccd)))
else:
    dataC.execute(rawQuery, (str(args.visit), ))


uniquenessCheck = dict()

for r in dataC.fetchall():
    sourceDir = f"/datasets/hsc/repo/{r['field']}/{r['taiObs']}/{r['filter']}/"
    if sourceDir not in sourceDirs:
        print(f"# Ingest: {sourceDir}")
        sourceDirs[sourceDir] = 1

    for table in ('bias', 'dark', 'flat', 'defects', 'fringe'):  # 'sky'):
        cq = calibQueryPre + table + calibQueryPost
        if args.ccd != '':
            calibC.execute(cq, (r['taiObs'], r['taiObs'], r['filter'], args.ccd))
        else:
            calibC.execute(cq, (r['taiObs'], r['taiObs'], r['filter']))

        calibRows = calibC.fetchall()
        print(table, len(calibRows))

        if len(calibRows) == 0:
            continue

        yCTupd = ['- type: dataset_type',
                  f"  name: {table}",
                  '  dimensions:',
                  '  - instrument',
                  '  - calibration_label',
                  '  - detector']
        yCT.extend(yCTupd)

        if table in ('flat', 'fringe', 'sky'):
            yCT.append('  - physical_filter')
        if table in ('bias', 'dark'):
            yCT.append('  storage_class: ImageF')
        elif table in ('flat', ):
            yCT.append('  storage_class: MaskedImageF')
        else:
            yCT.append('  storage_class: ExposureF')
        yCT.append(f"- type: dataset")
        yCT.append(f"  dataset_type: {table}")
        yCT.append(f"  run: {args.run}")
        yCT.append(f"  records:")

        for cr in calibRows:
            calibrationLabel = f"gen2/{table}_{cr['calibDate']}_{cr['id']}"

            if cr['filter'] != 'NONE':
                path = (f"/datasets/hsc/repo/CALIB/{table.upper()}/{cr['calibDate']}/{cr['filter']}"
                        f"/{table.upper()}-{cr['calibDate']}-{cr['filter']}-{cr['ccd']:03}.fits")
            else:
                path = (f"/datasets/hsc/repo/CALIB/{table.upper()}/{cr['calibDate']}/{cr['filter']}"
                        f"/{table.upper()}-{cr['calibDate']}-{cr['ccd']:03}.fits")

            if cr['filter'] in (r['filter'], 'NONE'):
                if calibrationLabel not in uniquenessCheck:
                    uniquenessCheck[calibrationLabel] = 1

                    print(f"calibration_label {calibrationLabel}")
                    yCL.append(f"  - instrument: {args.instrument}")
                    yCL.append(f"    name: {calibrationLabel}")
                    yCL.append(f"    datetime_begin: {cr['validStart']} 00:00:00")
                    yCL.append(f"    datetime_end: {cr['validEnd']} 00:00:00")
                    if table in ('flat', 'fringe', 'sky'):
                        yCL.append(f"    physical_filter: {cr['filter']}")

                print(f"dataset {path}")
                yCT.append(f"  - dataset_id: {args.datasetBase}")
                yCT.append(f"    data_id:")
                yCT.append(f"      instrument: {args.instrument}")
                yCT.append(f"      calibration_label: {calibrationLabel}")
                yCT.append(f"      detector: {cr['ccd']}")
                if table in ('flat', 'fringe', 'sky'):
                    yCT.append(f"      physical_filter: {cr['filter']}")

                yCT.append(f"    path: {path}")
                if table == 'defects':
                    yCT.append(f"    formatter: DEFECTS?")
                else:
                    yCT.append(f"    formatter: "
                               "lsst.daf.butler.formatters.fitsExposureFormatter.FitsExposureFormatter")
                args.datasetBase += 1

with open("./g2g3.yaml", "w") as f:
    outL = yL
    outL.extend(yCL)
    outL.extend(yCT)
    for l in outL:
        f.write(l + "\n")
