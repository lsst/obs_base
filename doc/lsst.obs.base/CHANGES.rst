obs_base v22.0 (2021-04-01)
===========================

New Feature
-----------

* Enhance raw data ingest such that there is no longer a need for a special subclass when ingesting DECam data.  The metadata translator can now find additional headers itself. [DM-29166]
* Add progress reporting to raw ingest, visit definition, and 2to3 conversion.
* Change raw data ingest to support remote object stores. [DM-25965]
* Raw data ingest now supports external metadata sidecar files or JSON per-directory index files. Creating these sidecar files in advance (using ``astrometadata write-index`` or ``astrometadata write-sidecar``) can significantly improve ingest performance. This is especially useful if a particular test data set is commonly re-ingested. [DM-27476]
* Raw data ingest has been modified to provide a callback feature when files fail to be ingested or are successfully ingested. This allows reporting tools to make detailed reports when doing bulk ingest. [DM-29071]
* 2to3 conversion has been significantly improved. [DM-27147]

Other
-----

* When reading exposures the formatter now checks that the filter label in the DataId is consistent with the filter label read from the file. [DM-28583]
