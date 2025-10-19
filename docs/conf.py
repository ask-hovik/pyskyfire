# -- Minimal Sphinx config for pyskyfire using AutoAPI + NumPyDoc ------------
from datetime import date
import os

project   = "pyskyfire"
author    = "Ask Haugerud Hovik"
copyright = f"{date.today().year}, {author}"

extensions = [
    "myst_parser",
    "autoapi.extension",
    "numpydoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.graphviz",
    "sphinx.ext.inheritance_diagram",
    # "sphinx.ext.autodoc.typehints",  # ← remove (or keep and set autodoc_typehints="none")
]

ROOT = os.path.abspath(os.path.join(__file__, "..", ".."))
SRC  = os.path.join(ROOT, "src")

autoapi_type  = "python"
autoapi_dirs  = [os.path.join(SRC, "pyskyfire")]
autoapi_root  = "autoapi"
autoapi_add_toctree_entry = True
autoapi_keep_files = True
autoapi_generate_api_docs = True
autoapi_member_order = "groupwise"
autoapi_python_class_content = "class"
autoapi_own_page_level = "class"
autoapi_python_use_implicit_namespaces = True
autoapi_options = [
    "members",
    "undoc-members",
    "private-members",
    "special-members",
    "imported-members",
    "show-inheritance",
    "show-inheritance-diagram",
    "show-module-summary",
]

# NumPyDoc
numpydoc_show_class_members = True
numpydoc_show_inherited_class_members = True
numpydoc_class_members_toctree = False
numpydoc_attributes_as_param_list = True
numpydoc_xref_param_type = True
# Optional: help short names link
numpydoc_xref_aliases = {
    "Station": "pyskyfire.common.engine_network.Station",
    "SignalBlock": "pyskyfire.common.blocks.SignalBlock",
}

# If you decided to keep the typehints extension:
# autodoc_typehints = "none"

myst_enable_extensions = ["colon_fence", "deflist", "substitution", "attrs_inline", "tasklist"]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy":  ("https://numpy.org/doc/stable/", None),
    "scipy":  ("https://docs.scipy.org/doc/scipy/", None),
}

graphviz_output_format = "svg"
graphviz_dot_args = ["-Gbgcolor=transparent", "-Gtransparent=true"]

source_suffix = {".md": "markdown", ".rst": "restructuredtext"}
html_theme   = "furo"
html_title   = "pyskyfire"
templates_path = ["_templates"]
html_static_path = ["_static"]

# Python changes
python_use_unqualified_type_names = True
# Render modern short generics (list[int] not typing.List[int])
python_display_short_literal_types = True

def skip_autoapi_members(app, what, name, obj, skip, options):
    if what == "attribute":
        return True
    if what == "method" and name.endswith(".__init__"):
        return True
    return skip

def setup(app):
    app.connect("autoapi-skip-member", skip_autoapi_members)
