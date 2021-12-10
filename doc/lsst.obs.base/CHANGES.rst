obs_base v23.0.0 2021-12-10
===========================

New Features
------------

- 2to3 conversion has been improved to add a dry run facility, to defer dataId expansion when not required, and to allow templates to be overridden. (`DM-28636 <https://jira.lsstcorp.org/browse/DM-28636>`_)
- Reorganize the base ``Exposure`` and raw formatters to improve efficiency and clarify component handling. (`DM-28698 <https://jira.lsstcorp.org/browse/DM-28698>`_)
- Add ``amp`` parameter to the formatters for the Exposure StorageClass, allowing single-amplifier subimage reads. (`DM-29370 <https://jira.lsstcorp.org/browse/DM-29370>`_)
- Change raw ingest to use a reproducible UUID5 dataset ID. This means that the dataset ID for a raw ingested in one repository will be identical to that used in another.  For integer-based registries this change will have no effect. (`DM-29950 <https://jira.lsstcorp.org/browse/DM-29950>`_)
- Add support for updating exposure and visit definitions in RawIngestTask and DefineVisitsTask. (`DM-30866 <https://jira.lsstcorp.org/browse/DM-30866>`_)
- Add support for forced updates of `instrument`, `detector`, and `physical_filter` definitions during instrument registration. (`DM-31903 <https://jira.lsstcorp.org/browse/DM-31903>`_)


Bug Fixes
---------

- Not all PSFs are persistable and now if one is encountered as part of composite disassembly it will be ignored. These types of PSFs were already silently dropped when writing a full ``Exposure``. (`DM-29794 <https://jira.lsstcorp.org/browse/DM-29794>`_)
- The ``butler define-visits`` command now correctly uses the ``--collections`` option to constrain the exposures that will be processed into visits. (`DM-31079 <https://jira.lsstcorp.org/browse/DM-31079>`_)


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
