"""Sphinx configuration file for an LSST stack package.

This configuration only affects single-package Sphinx documenation builds.
"""

from documenteer.conf.pipelinespkg import *  # noqa: F403

project = "obs_base"
html_theme_options["logotext"] = project  # noqa: F405
html_title = project
html_short_title = project
doxylink = {}
exclude_patterns = ["changes/*"]

# Try to pull in links for butler and pipe_base.
intersphinx_mapping["lsst"] = ("https://pipelines.lsst.io/v/daily/", None)  # noqa
