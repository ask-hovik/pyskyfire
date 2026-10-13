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
    "sphinx_codeautolink",
    "sphinx_tippy",
    # "sphinx.ext.autodoc.typehints",  # ← remove (or keep and set autodoc_typehints="none")
]

ROOT = os.path.abspath(os.path.join(__file__, "..", ".."))
SRC  = os.path.join(ROOT, "src")

autoapi_type  = "python"
autoapi_dirs  = [os.path.join(SRC, "pyskyfire")]
autoapi_root  = "autoapi"
# The index is Markdown, so AutoAPI cannot insert its own toctree entry; the
# Reference toctree in index.md lists the subpackages instead. Add new
# subpackages there, or their pages will be missing from the sidebar.
autoapi_add_toctree_entry = False
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
    "Station": "pyskyfire.network.engine_network.Station",
    "SignalBlock": "pyskyfire.network.blocks.SignalBlock",
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
# Generate heading ids so Markdown links like ``page.md#some-heading`` resolve.
myst_heading_anchors = 3

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy":  ("https://numpy.org/doc/stable/", None),
    "scipy":  ("https://docs.scipy.org/doc/scipy/", None),
}

# Code-block cross-references: names in example code link to their API page.
# Concatenating blocks lets later fragments of a page reuse earlier imports, and
# the preface covers literalinclude fragments that start below the import line.
codeautolink_concat_default = True
codeautolink_global_preface = "import pyskyfire as psf"
# Tippy shows the hover preview; the browser's own title tooltip would double it.
codeautolink_disable_link_titles = True

# Hover previews for internal links, including the ones codeautolink adds.
tippy_props = {"placement": "auto-start", "maxWidth": 600, "interactive": False}
tippy_enable_mathjax = True
# Wikipedia/DOI tips call external APIs during the build; keep builds offline.
tippy_enable_wikitips = False
tippy_enable_doitips = False

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

def orphan_autoapi_index(app, docname, source):
    """Mark AutoAPI's top-level pages as reachable by link only.

    The Reference toctree in index.md lists the subpackages directly, so the
    generated root and ``pyskyfire`` package pages sit outside every toctree.
    """
    if docname in (f"{autoapi_root}/index", f"{autoapi_root}/pyskyfire/index"):
        source[0] = ":orphan:\n\n" + source[0]

def copy_example_html_files(
    source: Path,
    destination: Path,
    required_files: set[str],
    exclude: set[str] = frozenset(),
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
        if html_file.name in exclude:
            continue
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
            / "advanced-simulation"
        ),
        required_files={
            "methane_engine_report.html",
            "contour.html",
            "network.html",
            "engine-3d.html"
        },
    )

    # Single-figure exports embedded by the capabilities page. Each one is
    # written by the case's own post-processor, so the page always shows the
    # current result rather than a checked-in snapshot of an older run.
    copy_example_html_files(
        source=repository_root / "examples" / "film_cooling" / "standalone",
        destination=output_static / "standalone" / "film-cooling",
        required_files={
            "film-regen-wall-temperature.html",
            "film-regen-heat-flux.html",
        },
    )

    copy_example_html_files(
        source=repository_root / "validation" / "RL10" / "standalone",
        destination=output_static / "standalone" / "rl10",
        # engine-3d.html is deliberately not listed: it is 30 MB and is
        # published once under /validation/ by publish_validation_reports.
        required_files={
            "network.html",
            "wall-temperature.html",
            "coolant-temperature.html",
            "coolant-pressure.html",
            "heat-flux.html",
            "temperature-profile-throat.html",
            "combustion-viscosity.html",
            "station-fuel-pressure.html",
            "station-fuel-temperature.html",
            "station-fuel-mass-flow.html",
        },
        exclude={"engine-3d.html"},
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

    copy_example_html_files(
        source=repository_root / "examples" / "pipe" / "output",
        destination=output_static / "howto-artifacts" / "pipes",
        required_files={
            "rounded-bend-path.html",
            "rounded-bend-pipe.html",
            "bend-comparison.html",
            "flared-bend-pipe.html",
            "pressure-loss-comparison.html",
            "flared-pressure.html",
            "flared-pressure-3d.html",
            "spline-curvature.html",
            "spline-pipe.html",
            "manual-path.html",
            "manual-curvature.html",
            "tapered-stations.html",
        },
    )

    copy_example_html_files(
        source=repository_root / "examples" / "area_preserving_helix",
        destination=(
            output_static
            / "howto-artifacts"
            / "area-preserving-helix"
        ),
        required_files={"engine-3d.html"},
    )


def publish_validation_reports(app, exception) -> None:
    """Replace validation stubs with their standalone HTML reports."""

    if exception is not None or app.builder.format != "html":
        return

    repository_root = Path(__file__).resolve().parent.parent
    rl10_dir = repository_root / "validation" / "RL10"
    source = rl10_dir / "RL10A-3-3A_Report.html"
    output_validation = Path(app.outdir) / "validation"
    destination = output_validation / "rl10a-3-3a.html"

    if not source.is_file():
        raise SphinxError(
            f"Missing generated validation report: {source}\n\n"
            "Run validation/RL10/post_process.py before building docs."
        )

    output_validation.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)

    # The report links to these as siblings, and the 3D viewer is also linked
    # directly from the README and the capabilities page. Publishing them next
    # to the report keeps both the in-report links and those references working.
    standalone_dir = rl10_dir / "standalone"
    companions = {"engine-3d.html", "network.html"}
    missing = sorted(
        name for name in companions if not (standalone_dir / name).is_file()
    )
    if missing:
        raise SphinxError(
            "Missing generated RL10 validation artifacts in "
            f"{standalone_dir}:\n"
            + "\n".join(f"  - {name}" for name in missing)
            + "\n\nRun validation/RL10/post_process.py before building docs."
        )

    for name in sorted(companions):
        shutil.copy2(standalone_dir / name, output_validation / name)


def add_code_link_tips(app, exception) -> None:
    """Give codeautolink's code-block links the same hover tips as prose links.

    sphinx-tippy collects links while pages render, but codeautolink only
    writes its links into the finished HTML at build-finished. This runs after
    both, adds the code links to tippy's per-page data and rewrites the tips.
    """

    if exception is not None or app.builder.format != "html":
        return

    import posixpath
    from bs4 import BeautifulSoup
    import sphinx_tippy

    pages = sphinx_tippy.get_tippy_data(app)["pages"]

    # A link to a whole AutoAPI object page would otherwise only show the
    # page title; show the object's signature and docstring summary instead.
    for pagename, data in pages.items():
        if pagename.startswith("autoapi/") and not pagename.endswith("/index"):
            object_id = pagename.removeprefix("autoapi/").replace("/", ".")
            if object_id in data["id_to_html"]:
                data["id_to_html"][None] = data["id_to_html"][object_id]

    for pagename, data in pages.items():
        html_file = Path(app.outdir, pagename + ".html")
        if not html_file.is_file():
            continue

        body = BeautifulSoup(html_file.read_text(encoding="utf8"), "html.parser")
        page_parent = posixpath.dirname(pagename)
        for anchor in body.select("a.sphinx-codeautolink-a[href]"):
            path, _, target = anchor["href"].partition("#")
            if not path.endswith(".html"):
                continue
            refpage = posixpath.normpath(posixpath.join(page_parent, path[:-5]))
            if refpage in app.env.all_docs:
                data["refs_in_page"].add((refpage, target or None))

        sphinx_tippy.write_tippy_props_page(app, pagename, {}, {}, {})


def setup(app):
    app.connect("autoapi-skip-member", skip_autoapi_members)
    app.connect("source-read", orphan_autoapi_index)
    app.connect("builder-inited", copy_docs_artifacts)
    app.connect("build-finished", publish_validation_reports)
    # Default priority is 500; codeautolink and sphinx-tippy must run first.
    app.connect("build-finished", add_code_link_tips, priority=900)
    # Hover tips typeset their math with MathJax, so every page needs it, not
    # only the pages that contain equations themselves.
    app.set_html_assets_policy("always")

# Logo
html_logo = "_static/pyskyfire_header.png"   # or .png

# Optional but recommended with Furo
html_theme_options = {
    "sidebar_hide_name": True,      # Hides the text "pyskyfire" when logo is shown
    # "light_logo": "_static/pyskyfire-logo-light.svg",   # if you have separate light/dark versions
    # "dark_logo": "_static/pyskyfire-logo-dark.svg",
}
