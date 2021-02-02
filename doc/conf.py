"""Sphinx configuration file for an LSST stack package.

This configuration only affects single-package Sphinx documenation builds.
"""

from documenteer.conf.pipelinespkg import *  # noqa: F403

project = "obs_base"
html_theme_options["logotext"] = project  # noqa: F405
html_title = project
html_short_title = project
doxylink = {}
