.. py:currentmodule:: lsst.obs.base

.. _lsst.obs.base:

#############
lsst.obs.base
#############

The ``lsst.obs.base`` module provides the framework and common API for telescope/camera specific "obs" packages.
New cameras will derive from the classes defined here.
`lsst.obs.base.tests` provides the tests that all "obs" packages should pass.

.. _lsst.ctrl.mpexec-changes:

Changes
=======

.. toctree::
   :maxdepth: 1

   CHANGES.rst

.. _lsst.obs.base-using:

Using lsst.obs.base
===================

.. toctree::
   :maxdepth: 1

   creating-an-obs-package

.. _lsst.obs.base-contributing:

Contributing
============

``lsst.obs.base`` is developed at https://github.com/lsst/obs_base.
You can find Jira issues for this module under the `obs_base <https://jira.lsstcorp.org/issues/?jql=project%20%3D%20DM%20AND%20component%20%3D%20obs_base>`_ component.

.. _lsst.obs.base-cli:

Command Line Interface
======================

``daf_butler`` implements a command line interface command called ``butler``. The following subcommands are
implemented by this package and available to the ``butler`` command when this package is setup.

.. click:: lsst.obs.base.cli.doc.butlerCmdDocGen:cli
   :prog: butler
   :show-nested:

.. _lsst.obs.base-pyapi:

Python API reference
====================

.. automodapi:: lsst.obs.base
   :no-main-docstr:

.. automodapi:: lsst.obs.base.tests
   :no-main-docstr:

.. automodapi:: lsst.obs.base.formatters.fitsExposure
   :no-main-docstr:

.. automodapi:: lsst.obs.base.formatters.fitsGeneric
   :no-main-docstr:
