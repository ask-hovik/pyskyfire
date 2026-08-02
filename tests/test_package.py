"""Smoke tests for the installed Pyskyfire distribution."""

from importlib import metadata, resources
import subprocess
import sys

import pyskyfire
import pytest


PACKAGED_RESOURCES = (
    "regen/data/theta_e.json",
    "regen/data/theta_n.json",
    "viz/assets/GLTFLoader.js",
    "viz/assets/OrbitControls.js",
    "viz/assets/three.min.js",
    "viz/assets/x6.min.js",
)


def test_installed_package_version() -> None:
    """The import must come from an installed distribution with metadata."""

    installed_version = metadata.version("pyskyfire")

    assert pyskyfire.__version__ == installed_version
    assert installed_version != "0+unknown"


def test_core_import_does_not_load_visualisation() -> None:
    """The base package must remain importable without loading viz libraries."""

    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; import pyskyfire; "
            "assert 'pyskyfire.viz' not in sys.modules; "
            "assert 'plotly' not in sys.modules; "
            "assert 'pyvista' not in sys.modules",
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr


def test_missing_visualisation_dependencies_have_install_hint() -> None:
    """Accessing optional visualization gives users an actionable error."""

    script = """
import importlib.util

real_find_spec = importlib.util.find_spec
importlib.util.find_spec = lambda name, *args, **kwargs: (
    None if name in {"plotly", "pyvista"} else real_find_spec(name, *args, **kwargs)
)

import pyskyfire

try:
    pyskyfire.viz
except ImportError as error:
    message = str(error)
    assert "pyskyfire[viz]" in message
    assert "plotly" in message
    assert "pyvista" in message
else:
    raise AssertionError("Importing pyskyfire.viz should have failed")
"""
    result = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("relative_path", PACKAGED_RESOURCES)
def test_packaged_resource_exists(relative_path: str) -> None:
    """Runtime data and browser assets must be present in the built wheel."""

    resource = resources.files("pyskyfire").joinpath(relative_path)

    assert resource.is_file(), f"Missing packaged resource: {relative_path}"
