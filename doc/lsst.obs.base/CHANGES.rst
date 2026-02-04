obs_base v30.0.0 (2026-01-16)
=============================

New Features
------------

- Added ``lsst.obs.base.utils.TableVStack`` for efficient Astropy Table concatenation from dataset refs.

  This utility class was made for consolidating object tables in ``pipe_tasks`` but is broadly useful.
  While ``obs_base`` isn't the most logical place to put it, it's the lowest level package it can go right now. (`DM-52169 <https://rubinobs.atlassian.net/browse/DM-52169>`_)


Bug Fixes
---------

- The ``FitsExposureFormatter`` now forces the data type returned to match the storage class, overriding the default type associated with the FITS file. (`DM-52873 <https://rubinobs.atlassian.net/browse/DM-52873>`_)


Other Changes and Additions
---------------------------

- Moved camera geometry parity flips (relative to sky) into camera geometry (between the ``FOCAL_PLANE`` and ``FIELD_ANGLE`` coordinate systems).

  Raw formatters still track whether a parity flip is needed in order to support old camera definitions.
  This may change in the future. (`DM-20746 <https://rubinobs.atlassian.net/browse/DM-20746>`_)
- Added ``transmission_filter_detector`` to curated calibration types. (`DM-46526 <https://rubinobs.atlassian.net/browse/DM-46526>`_)
- Added new ``lsst.obs.base.utils.iso_date_to_curated_calib_file_root`` to unify the logic for converting an ISO date string to a curated calibrations file name. (`DM-49900 <https://rubinobs.atlassian.net/browse/DM-49900>`_)
- Adjusted the logic for visit definition such that a visit will not be defined if it is known that the observation was not looking at the sky (such as with closed-dome tests). (`DM-50167 <https://rubinobs.atlassian.net/browse/DM-50167>`_)
- When defining visits, always update visit-detector regions if ``update_records=True``.

  Previously, we only updated visit-detector regions when the visit record itself changed, which made it impossible to fix cases where the code for the detector regions was changing but the visit record was not, as well as more mysterious cases where some visit-detector regions were not originally inserted at all. (`DM-50446 <https://rubinobs.atlassian.net/browse/DM-50446>`_)
- Dropped support for old dimension universes (and stopped requiring ``raw`` dataset existence)  in ``butler define-visits``. (`DM-50661 <https://rubinobs.atlassian.net/browse/DM-50661>`_)
- Adapt ``FitsExposureFormatter`` to the new interfaces for FITS compression. (`DM-52879 <https://rubinobs.atlassian.net/browse/DM-52879>`_)
- Added ``--skip-existing`` command line option to ``butler ingest-raws``.
  Skip existing raws is now the default (previously re-ingesting was an error). (`DM-53163 <https://rubinobs.atlassian.net/browse/DM-53163>`_)
- Capped the number of input provenance datasets recorded in Exposure FITS files to 3,000 inputs. (`DM-53326 <https://rubinobs.atlassian.net/browse/DM-53326>`_)
- Modified ``DefineVisitsTask.run`` to return a struct of visit counts. (`DM-53622 <https://rubinobs.atlassian.net/browse/DM-53622>`_)
- Added timer log messages to raw ingest to allow reporting of metadata gathering and butler ingest times separately. (`DM-53679 <https://rubinobs.atlassian.net/browse/DM-53679>`_)


obs_base v29.0.0 (2025-03-26)
=============================

New Features
------------

- Modified the FITS-based formatters to write out dataset provenance in the output FITS header. (`DM-35396 <https://rubinobs.atlassian.net/browse/DM-35396>`_)
- Added a new dataset type for defects that are manually defined. (`DM-47365 <https://rubinobs.atlassian.net/browse/DM-47365>`_)

Other Changes and Additions
---------------------------

- Fixed problem ingesting on-sky raws into a Postgres database with numpy 2. (`DM-49845 <https://rubinobs.atlassian.net/browse/DM-49845>`_)

obs_base v28.0.0 (2024-11-21)
=============================

API Changes
-----------

- The ``butler write-curated-calibrations`` command now requires at least one "label" to be included for the collection name.

  This prevents the common mistake of setting up ``<instrument>/calib`` as a ``CALIBRATION`` collection rather than a more maintainable ``CHAINED`` collection. (`DM-46297 <https://rubinobs.atlassian.net/browse/DM-46297>`_)


obs_base v27.0.0 (2024-06-05)
=============================

New Features
------------

- * Added support for ingesting raw data into a version 6 dimension universe.
    This universe includes ``day_obs`` and ``group`` as dimensions.
  * Modified visit definition to allow a visit to be defined by ``group`` dimension.
  * Added ``Instrument.translatorClass`` class property that can be used to specify the relevant ``astro_metadata_translator.MetadataTranslator`` to use.
  * Added ``Instrument.group_name_to_group_id`` method to convert a group name string to an integer suitable for use as a visit ID. (`DM-42636 <https://rubinobs.atlassian.net/browse/DM-42636>`_)
- Added a ``--update-records`` option to the ``butler ingest-raws`` command-line tool.
  This can be used if there has been a change in the metadata translator resulting in a change of definition of the ``exposure`` record.
  Only use this option if you understand why a change has occurred. (`DM-43135 <https://rubinobs.atlassian.net/browse/DM-43135>`_)


- Corrected and clarified docstrings for ``read_curated_calibs``, ``read_one_calib`` and ``read_all`` functions, ensuring variable name consistency between these two and ``check_metadata``. (`DM-22465 <https://rubinobs.atlassian.net/browse/DM-22465>`_)

obs_base v26.0.0 (2023-08-03)
=============================

New Features
------------

- Added support for defining visits incrementally as exposures are ingested.
  This allows files from the telescope to be ingested one at a time whilst redefining the existing visits.
  Additionally ``--update-records`` and ``--incremental`` have been added to the ``butler define-visits`` command-line. (`DM-36395 <https://rubinobs.atlassian.net/browse/DM-36395>`_)
- * Added ``transmission_`` curve dataset types to the set of curated calibrations.
  * Updated the curated calibration code in the ``Instrument`` definition to allow for flexibility in required dimensions.
  * Updated the read curated calibration code to allow for the same flexibility in dimensions. (`DM-36597 <https://rubinobs.atlassian.net/browse/DM-36597>`_)
- Raw ingest can now ask the ``Instrument`` class for the raw dataset type definition.
  This means it is no longer required to subclass the ``getDatasetType`` method and allows various instruments to be ingested with the base class implementation. (`DM-37950 <https://rubinobs.atlassian.net/browse/DM-37950>`_)
- ``DefineVisitsTask`` now calls ObsCore table manager to update exposure regions after visit is defined.
  New configuration field ``updateObsCoreTable`` for that task can be set to `False` to disable exposure updates. (`DM-38205 <https://rubinobs.atlassian.net/browse/DM-38205>`_)

An API Removal or Deprecation
-----------------------------

- Deprecated ``ExposureIdInfo`` in favor of ``lsst.meas.base.IdGenerator``. (`DM-31924 <https://rubinobs.atlassian.net/browse/DM-31924>`_)
- Removed deprecated ``getInstrument`` function from ``lsst.obs.base.utils``. (`DM-37534 <https://rubinobs.atlassian.net/browse/DM-37534>`_)

Bug Fixes
---------

- Fixed curated calibration reading to check parent directory if there are no sub-directories. (`DM-36598 <https://rubinobs.atlassian.net/browse/DM-36598>`_)

Other Changes and Additions
---------------------------

- Modified the raw ingest task to use resolved ``DatasetRef``. (`DM-38779 <https://rubinobs.atlassian.net/browse/DM-38779>`_)

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
