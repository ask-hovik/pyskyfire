# -- Minimal Sphinx config for pyskyfire using AutoAPI + NumPyDoc ------------
from datetime import date
import os
from pathlib import Path
import shutil
from sphinx.errors import SphinxError

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
    "sphinx.ext.mathjax",
    # "sphinx.ext.autodoc.typehints",  # ← remove (or keep and set autodoc_typehints="none")
]

ROOT = os.path.abspath(os.path.join(__file__, "..", ".."))
SRC  = os.path.join(ROOT, "src")

autoapi_type  = "python"
autoapi_dirs  = [os.path.join(SRC, "pyskyfire")]
autoapi_root  = "autoapi"
autoapi_add_toctree_entry = True
autoapi_keep_files = False
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
    #"imported-members",
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

autosummary_generate = False
# If you decided to keep the typehints extension:
# autodoc_typehints = "none"

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "substitution",
    "attrs_inline",
    "tasklist",
    "dollarmath", 
    "amsmath",         
]

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
html_css_files = [
    "custom.css",
]

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

def copy_example_html_files(
    source: Path,
    destination: Path,
    required_files: set[str],
) -> None:
    """Copy flat, pre-generated example HTML files into the docs build."""

    missing = sorted(
        filename for filename in required_files
        if not (source / filename).is_file()
    )
    if missing:
        raise SphinxError(
            f"Missing generated documentation files in {source}:\n"
            + "\n".join(f"  - {filename}" for filename in missing)
            + "\n\nRun the corresponding example manually before building docs."
        )

    if destination.exists():
        shutil.rmtree(destination)

    destination.mkdir(parents=True)

    for html_file in source.glob("*.html"):
        shutil.copy2(html_file, destination / html_file.name)


def copy_docs_artifacts(app) -> None:
    """Copy pre-generated HTML tutorial artifacts for HTML builds."""

    if app.builder.format != "html":
        return

    repository_root = Path(__file__).resolve().parent.parent
    output_static = Path(app.outdir) / "_static"

    copy_example_html_files(
        source=repository_root / "examples" / "minimal",
        destination=(
            output_static
            / "tutorial-artifacts"
            / "minimal-simulation"
        ),
        required_files={
            "minimal-report.html",
            "engine-3d.html",
            "contour.html",
            "heat-flux.html",
        },
    )

    copy_example_html_files(
        source=repository_root / "examples" / "advanced",
        destination=(
            output_static
            / "tutorial-artifacts"
            / "advanced-cycle-simulation"
        ),
        required_files={
            "methane_engine_report.html",
        },
    )

    copy_example_html_files(
        source=repository_root / "examples" / "MR_optimisation",
        destination=(
            output_static
            / "howto-artifacts"
            / "mixture-ratio-optimization"
        ),
        required_files={
            "mixture-ratio-optimisation.html",
            "optimized-contour.html",
        },
    )

def setup(app):
    app.connect("autoapi-skip-member", skip_autoapi_members)
    app.connect("builder-inited", copy_docs_artifacts)

# Logo
html_logo = "_static/pyskyfire_header.png"   # or .png

# Optional but recommended with Furo
html_theme_options = {
    "sidebar_hide_name": True,      # Hides the text "pyskyfire" when logo is shown
    # "light_logo": "_static/pyskyfire-logo-light.svg",   # if you have separate light/dark versions
    # "dark_logo": "_static/pyskyfire-logo-dark.svg",
}