from datetime import datetime
import os
import sys

# Make the package importable for autodoc
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../src"))

project = "pyskyfire"
author = "Ask Haugerud Hovik"
copyright = f"{datetime.now():%Y}, {author}"

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",          # Google/NumPy docstrings
    "sphinx.ext.intersphinx",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinx_autodoc_typehints",
]

# Autodoc / Autosummary
autosummary_generate = True
autodoc_typehints = "description"   # cleaner signatures in docs
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_use_param = True
napoleon_use_rtype = True

# Cross-link to common libs
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}

# Markdown (MyST)
myst_enable_extensions = [
    "colon_fence", "deflist", "substitution", "tasklist", "attrs_inline"
]
source_suffix = {".md": "markdown", ".rst": "restructuredtext"}

html_theme = "furo"
html_title = "pyskyfire"
html_static_path = ["_static"]
templates_path = ["_templates"]
