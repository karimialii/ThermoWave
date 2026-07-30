"""Sphinx configuration for ThermoWave's Read the Docs build."""

import os
import sys
from importlib.metadata import PackageNotFoundError, version as _pkg_version

sys.path.insert(0, os.path.abspath("../src"))

project = "ThermoWave"
copyright = "2026, ThermoWave contributors"
author = "ThermoWave contributors"
try:
    release = _pkg_version("thermowave")
except PackageNotFoundError:
    release = "0.0.0"
version = ".".join(release.split(".")[:2])

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "myst_parser",
]

myst_enable_extensions = [
    "dollarmath",
    "amsmath",
]
myst_heading_anchors = 3

autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}
autodoc_typehints = "description"
napoleon_google_docstring = True
napoleon_numpy_docstring = False

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "superpowers"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
