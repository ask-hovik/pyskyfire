from datetime import datetime
import os
import sys

# Make the package importable for autodoc
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../src"))

ROOT = os.path.abspath(os.path.join(__file__, "..", ".."))
SRC  = os.path.join(ROOT, "src")
_AUTOSUMMARY_DIR = os.path.join(os.path.dirname(__file__), "api", "_generated")
os.makedirs(_AUTOSUMMARY_DIR, exist_ok=True)

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
    "autoapi.extension"
]

autoapi_type = "python"
autoapi_dirs = [SRC]              # points to your 'src' path
autoapi_add_toctree_entry = True
autoapi_root = "api"              # where to put the rendered API
autoapi_keep_files = True         # keep intermediate .rst (useful for debugging/PRs)
autoapi_python_class_content = "both"  # class + __init__ docstrings
autoapi_options = [
    "members",
    "undoc-members",
    "show-inheritance",
    "inherited-members",
    "imported-members",  # show re-exported names (works with __all__)
]

# Autodoc / Autosummary
autosummary_generate = True
autodoc_typehints = "description"   # cleaner signatures in docs
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_use_param = True
napoleon_use_rtype = True

autodoc_mock_imports = [
    "CoolProp",
    "CEA_Wrap",
    "cantera",
    "centrifugal_pump",
    "cloudpickle",
    "gmsh",
    "plotly",
    "pyvis",
]


# Type hints: in the signature + description
typehints_fully_qualified = False
always_document_param_types = True


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

# Autodoc defaults (good starting point)
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "inherited-members": True,
    "show-inheritance": True,
}
autoclass_content = "class"  # or "both" to include __init__ docstrings, too


source_suffix = {".md": "markdown", ".rst": "restructuredtext"}

html_theme = "furo"
html_title = "pyskyfire"
html_static_path = ["_static"]
templates_path = ["_templates"]

import os, sys
sys.path.insert(0, os.path.abspath(".."))  # if docs/ is next to package root
sys.path.insert(0, os.path.abspath("../src"))
