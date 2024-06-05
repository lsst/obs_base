obs_base v25.0.0 (2023-03-02)
=============================

New Features
------------

- * Removed Gen2 code from package, including 2to3 conversion code.
    Use an older release to convert any remaining Gen2 repositories to Gen3.
  * Moved support for reading curated calibrations for ``obs_*_data`` packages from ``pipe_tasks`` and added unit tests for this code. (`DM-35035 <https://rubinobs.atlassian.net/browse/DM-35035>`_)
- * Added ``focusZ`` to ``MakeRawVisitInfoViaObsInfo``. (`DM-35186 <https://rubinobs.atlassian.net/browse/DM-35186>`_)


API Changes
-----------

- The ``ingest-raws`` and ``define-visits`` subcommands no longer allow multiple config settings within a single ``--config`` option.
  We have decided that it is too dangerous to split on comma in the general case and so have removed that facility to be consistent with other commands.
  Use multiple ``--config`` options instead. (`DM-35917 <https://rubinobs.atlassian.net/browse/DM-35917>`_)


An API Removal or Deprecation
-----------------------------

- Removed deprecated ``filterLabel`` exposure component access. (`DM-27811 <https://rubinobs.atlassian.net/browse/DM-27811>`_)


obs_base v24.0.0 (2022-08-30)
=============================

New Features
------------

- * Visits will now be defined for all on-sky observations regardless of the observation type.
  * Changed ``butler define-visits`` so that it now takes a ``--where`` option.
    This can be used to restrict the visit definition to specific exposures. (`DM-33848 <https://rubinobs.atlassian.net/browse/DM-33848>`_)
- Add a ``--fail-fast`` option to ``butler ingest-raws`` (`DM-33891 <https://rubinobs.atlassian.net/browse/DM-33891>`_)
- * Modify ``ingest-raws`` to support new schema for exposure records.
  * Change ``define-visits`` to support new and old schema.
  * Visit system is now an enum rather than a configuration value.
  * Add new visit system to group by ``seq_start`` and ``seq_end`` and also allocate one-to-one visits for every exposure. (`DM-33942 <https://rubinobs.atlassian.net/browse/DM-33942>`_)
- * `lsst.obs.base.Instrument` is now a subclass of `lsst.pipe.base.Instrument`. This simplifies the dependencies of ``ctrl_mpexec`` by removing any need to understand camera geometry or curated calibrations.
  * As part of this move the ``butler register-instrument`` command has been moved to ``pipe_base``.
  * The ``PackagesFormatter`` has been moved to ``daf_butler`` and the ``PexConfigFormatter`` has been moved to ``pipe_base`` since both of these are required by ``ctrl_mpexec``.
  * ``lsst.obs.base.utils.getInstrument`` has been replaced with ``Instrument.from_string``. (`DM-34105 <https://rubinobs.atlassian.net/browse/DM-34105>`_)
- * Made choice of required ``ObservationInfo`` properties configurable
    through ``RawIngestTask.getObservationInfoSubsets``.
  * Added the concept of "dependency" records to be added to the registry before
    adding the exposure record; this makes it easier to satisfy foreign key
    constraints when the exposure relates to dimensions beyond the standard set.
  * Added ``RawIngestTask`` methods ``makeExposureRecord`` and ``makeDependencyRecords``
    to provide hooks for subclasses to provide values for additional columns. (`DM-34175 <https://rubinobs.atlassian.net/browse/DM-34175>`_)


API Changes
-----------

- Add a new option ``--track-file-attrs`` to ``butler ingest-raws``.
  This controls whether the ingested files should have file sizes and checksums tracked by the datastore.
  Use ``--no-track-files-attrs`` to disable size tracking. (`DM-33086 <https://rubinobs.atlassian.net/browse/DM-33086>`_)


An API Removal or Deprecation
-----------------------------

- `~lsst.obs.base.FilterDefinition` no longer supports ``lsst.afw.image.Filter``.
  The ``defineFilters`` and ``reset`` methods have been removed, as have all wavelength parameters to the `~lsst.obs.base.FilterDefinition` constructor.

  The old ``filter`` component for exposures has been removed, and replaced with a new ``filter`` component backed by ``lsst.afw.image.FilterLabel``.
  It functions identically to the ``filterLabel`` component, which has been deprecated. (`DM-27177 <https://rubinobs.atlassian.net/browse/DM-27177>`_)
- Remove the ``processes`` and ``pool`` arguments and the ``--processes`` command-line argument from `lsst.obs.base.DefineVisitsTask.run` and ``butler define-visits`` (respectively).
  These were already broken for ``processes > 1``, and internal parallelization here is no longer useful now that this task just does database I/O, not raw metadata reads. (`DM-33783 <https://rubinobs.atlassian.net/browse/DM-33783>`_)


obs_base v23.0.0 (2021-12-10)
=============================

New Features
------------

- 2to3 conversion has been improved to add a dry run facility, to defer dataId expansion when not required, and to allow templates to be overridden. (`DM-28636 <https://rubinobs.atlassian.net/browse/DM-28636>`_)
- Reorganize the base ``Exposure`` and raw formatters to improve efficiency and clarify component handling. (`DM-28698 <https://rubinobs.atlassian.net/browse/DM-28698>`_)
- Add ``amp`` parameter to the formatters for the ``Exposure`` `~lsst.daf.butler.StorageClass`, allowing single-amplifier subimage reads. (`DM-29370 <https://rubinobs.atlassian.net/browse/DM-29370>`_)
- Change raw ingest to use a reproducible UUID5 dataset ID. This means that the dataset ID for a raw ingested in one repository will be identical to that used in another.  For integer-based registries this change will have no effect. (`DM-29950 <https://rubinobs.atlassian.net/browse/DM-29950>`_)
- Add support for updating exposure and visit definitions in `~lsst.obs.base.RawIngestTask` and `~lsst.obs.base.DefineVisitsTask`. (`DM-30866 <https://rubinobs.atlassian.net/browse/DM-30866>`_)
- Add support for forced updates of ``instrument``, ``detector``, and ``physical_filter`` definitions during instrument registration. (`DM-31903 <https://rubinobs.atlassian.net/browse/DM-31903>`_)


Bug Fixes
---------

- Not all PSFs are persistable and now if one is encountered as part of composite disassembly it will be ignored. These types of PSFs were already silently dropped when writing a full ``Exposure``. (`DM-29794 <https://rubinobs.atlassian.net/browse/DM-29794>`_)
- The ``butler define-visits`` command now correctly uses the ``--collections`` option to constrain the exposures that will be processed into visits. (`DM-31079 <https://rubinobs.atlassian.net/browse/DM-31079>`_)


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
