.. py:currentmodule:: lsst.obs.base

.. creating-an-obs-package:

####################################
How to create a new LSST obs package
####################################

These instructions describe how to create a package that allows the LSST Science Pipelines to read raw data from a telescope and instrument using the gen3 middleware and ingest it into an LSST repository.
"Ingestion" is the process of reading raw files, interpreting their metadata, and formatting the data in a way the :py:mod:`Data Butler <lsst.daf.butler>` can read.

This guide describes how to create an interface between the :py:mod:`Data Butler <lsst.daf.butler>` and your raw files, so that users do not have to know any details about a particular instrument.
These instructions are necessary, but not sufficient: camera geometry, config overrides for tasks, and defects and other curated calibration are other necessary components not described here.
You will have successfully made an obs package when an :ref:`ingestion test <testing>` passes.
This demonstrates that the butler can read data ingested from your instrument.


The instructions here assume that you are writing for the ``ExampleCam`` camera, with the new package called ``obs_example``.
Here we put the code in a ``lsst.obs`` module hierarchy but this is not required by the Data Butler and you can use any hierarchy that suits your needs.

.. _translator:

MetadataTranslator
==================

The `astro_metadata_translator` package abstracts away reading metadata from raw images of a given telescope+instrument.
For example, the time of the observation, temperature, and filter for a given observation might each be stored in each raw file's FITS header as ``"OBSTIME"``, ``"TEMP_W"``, ``"FILTOBS"``.
The translator knows how to get all of this relevant metadata out of the raw files, so that they can be loaded into the :py:mod:`Data Butler <lsst.daf.butler>`.

Creating a ``MetadataTranslator`` derived from `~astro_metadata_translator.translators.MetadataTranslator` or `~astro_metadata_translator.translators.FitsTranslator` is the first step in building an obs_package.
Your new metadata translator can live in a separate package as you are developing it, but you should eventually make a pull request to `astro_metadata_translator` itself, so that other people can readily use it.
You will call your new translator ``ExampleTranslator`` from here on.

See the `metadata translator <https://astro-metadata-translator.lsst.io>`_ package for details.

.. _filters:

Filters
=======

Every instrument has a particular set of photometric filters, each with their own effective wavelength (e.g., ``477 microns``), transmission function, internal name (e.g., ``HSC-G``), and filter band (e.g., ``g``).
You define these filters using a `FilterDefinitionCollection` containing multiple :py:class:`FilterDefinitions <FilterDefinition>`.

Create a file, ``python/lsst/obs/example/exampleFilters.py``, formatted like the following:

.. code-block:: python

    from lsst.obs.base import FilterDefinition, FilterDefinitionCollection

    EXAMPLE_FILTER_DEFINITIONS = FilterDefinitionCollection(
        FilterDefinition(physical_filter="example g filter",
                         band="g",
                         lambdaEff=432),
        FilterDefinition(physical_filter="example z filter",
                         band="z",
                         lambdaEff=1234),
    )

See the `FilterDefinition` docs for the various components that go into defining a filter.
Note that the ``physical_filter`` name should match the exact filter name used in that observatory's metadata.

.. _formatter:

Formatter
=========

A `~lsst.daf.butler.Formatter` defines how the Data Butler reads raw data from the original files and converts it to the LSST conventions.
At present, most astronomical data is distributed as FITS files, so you will very likely be creating a subclass of `FitsRawFormatterBase`, but this is not required.

Create a file, ``python/lsst/obs/example/rawFormatter.py``, containing an ``ExampleRawFormatter`` derived from the `FitsRawFormatterBase` base class.
At minimum, you must define a ``translatorClass`` pointing to the ``ExampleTranslator`` you made in :ref:`translator` and a ``filterDefinitions`` pointing to the filter list you created in :ref:`filters`.
This formatter can also contain specializations for your specific camera, for example to rotate the `~lsst.afw.geom.SkyWcs` for detectors in the camera that have a different orientation.

.. code-block:: python

    __all__ = ["ExampleCameraRawFormatter"]

    from astro_metadata_translator import ExampleTranslator
    from lsst.obs.base import FitsRawFormatterBase
    from .exampleFilters import EXAMPLE_FILTER_DEFINITIONS


    class ExampleCameraRawFormatter(FitsRawFormatterBase):
        translatorClass = ExampleTranslator
        filterDefinitions = EXAMPLE_FILTER_DEFINITIONS

        def getDetector(self, id):
            return ExampleCamera().getCamera()[id]

.. _instrument:

Instrument
==========

An `Instrument` defines the instrument-specific logic for the Data Butler.

First create a new file ``tests/test_instrument.py`` with a test case derived from `~lsst.obs.base.instrument_tests.InstrumentTests` and `~lsst.utils.tests.TestCase`, defining ``self.data`` and ``self.instrument`` in ``setUp``.
The `set` of ``physical_filters`` you provide here will be checked to ensure that your `FilterDefinitionCollection` is loaded correctly.

.. code-block:: python

    """Tests of the ExampleCam instrument class.
    """

    import unittest

    import lsst.utils.tests
    import lsst.obs.example
    from lsst.obs.base.instrument_tests import InstrumentTests, InstrumentTestData


    class TestExampleCam(InstrumentTests, lsst.utils.tests.TestCase):
        def setUp(self):
            physical_filters = {"example g filter",
                                "example z filter"}

            self.data = InstrumentTestData(name="Example",
                                           nDetectors=4,
                                           firstDetectorName="1_1",
                                           physical_filters=physical_filters)
            self.instrument = lsst.obs.example.ExampleCam()

    if __name__ == '__main__':
        lsst.utils.tests.init()
        unittest.main()

Run this test via

.. code-block:: bash

    pytest -sv tests/test_instrument.py

the tests should fail, as there is no Example `Instrument` yet.

Next, add a file in ``python/lsst/obs/example/_instrument.py`` containing a subclass of `Instrument`, named ```ExampleCam``, which at minimum overrides these abstract methods: `Instrument.getName`, `Instrument.getCamera`, `Instrument.register`, `Instrument.filterDefinitions`, `Instrument.getRawFormatter` and define ``self.configPaths`` in ``__init__``.
The underscore is used in the name to indicate that the class will be exported by default and referred to as ``lsst.obs.example.ExampleCam``.

Run your test again: the tests should now pass.
If they do not, you can use the test output to determine what parts of the Instrument need to be fixed.

.. _testing:

Ingest tests
============

In order to test how your new gen3 obs package works with the :py:mod:`Data Butler <lsst.daf.butler>`, you need to write a test that ingests raw data.
`~lsst.obs.base.ingest_tests.IngestTestBase` provides a base class for those tests, requiring only that you specify the input data that will be tested, and the :ref:`dataIds <lsst.daf.butler-dimensions_data_ids>` to use to check that the data was correctly ingested.
This is how our system tests that your ``Formatter`` works correctly and that the ingest process can extract the required metadata from the files.

.. code-block:: python

    """Unit tests for Gen3 ExampleCam raw data ingest.
    """

    import unittest
    import os
    import lsst.utils.tests

    from lsst.obs.base.ingest_tests import IngestTestBase
    from lsst.obs.example.hsc import ExampleCam

    testDataPackage = "testdata_example"
    try:
        testDataDirectory = lsst.utils.getPackageDir(testDataPackage)
    except lsst.pex.exceptions.NotFoundError:
        testDataDirectory = None


    @unittest.skipIf(testDataDirectory is None, "testdata_example must be set up")
    class ExampleIngestTestCase(IngestTestBase, lsst.utils.tests.TestCase):
        def setUp(self):
            self.ingestdir = os.path.dirname(__file__)
            self.instrument = Examplecam()
            self.file = os.path.join(testDataDirectory, "example", "raw", "somefile.fits.gz")
            self.dataId = dict(instrument="Example", exposure=12345, detector=123)

            super().setUp()


    def setup_module(module):
        lsst.utils.tests.init()


    if __name__ == "__main__":
        lsst.utils.tests.init()
        unittest.main()


The ingest tests do not check pixel values, so it is acceptable to run the ingest on stripped data files where the pixel values have been set to a single value and the data compressed with ``fpack``.
This can result in a very small file that can be included directly in your obs package.
