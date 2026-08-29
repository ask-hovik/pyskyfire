"""Run the public examples as isolated subprocess smoke tests."""

from collections.abc import Iterable
import os
from pathlib import Path
import subprocess
import sys

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TIMEOUT_SECONDS = 30 * 60


def run_example(
    script: str,
    output_dir: Path,
    expected_files: Iterable[str],
    *,
    extra_args: Iterable[str | Path] = (),
    timeout: int = DEFAULT_TIMEOUT_SECONDS,
) -> None:
    """Run one example and verify that it produced non-empty artifacts."""

    output_dir.mkdir(parents=True, exist_ok=True)

    environment = os.environ.copy()
    environment["MPLBACKEND"] = "Agg"
    environment["MPLCONFIGDIR"] = str(output_dir / ".matplotlib")
    environment["PYVISTA_OFF_SCREEN"] = "true"
    environment["PYTHONUNBUFFERED"] = "1"

    command = [
        sys.executable,
        str(PROJECT_ROOT / script),
        "--output-dir",
        str(output_dir),
        *(str(argument) for argument in extra_args),
    ]

    try:
        result = subprocess.run(
            command,
            cwd=output_dir,
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as error:
        output = error.stdout or ""
        pytest.fail(
            f"{script} exceeded its {timeout}-second timeout.\n{output}",
            pytrace=False,
        )

    if result.returncode != 0:
        pytest.fail(
            f"{script} exited with status {result.returncode}.\n{result.stdout}",
            pytrace=False,
        )

    for relative_path in expected_files:
        artifact = output_dir / relative_path
        assert artifact.is_file(), f"{script} did not create {relative_path}"
        assert artifact.stat().st_size > 0, f"{relative_path} is empty"


def test_minimal_example(tmp_path: Path) -> None:
    """The minimal example is the per-push simulation smoke test."""

    run_example(
        "examples/minimal/minimal_sim.py",
        tmp_path / "minimal",
        (
            "minimal-report.html",
            "engine-3d.html",
            "contour.html",
            "heat-flux.html",
        ),
        timeout=15 * 60,
    )


def test_mixture_ratio_optimisation_example(tmp_path: Path) -> None:
    run_example(
        "examples/MR_optimisation/MR_opt.py",
        tmp_path / "mixture-ratio",
        ("mixture-ratio-optimisation.html", "optimized-contour.html"),
    )


def test_boiling_example(tmp_path: Path) -> None:
    run_example(
        "examples/boiling/boiling_sim.py",
        tmp_path / "boiling",
        ("boiling-report.html",),
    )


def test_curvature_correction_example(tmp_path: Path) -> None:
    run_example(
        "examples/curvature_correction/curvature_correction.py",
        tmp_path / "curvature-correction",
        ("curvature-correction-comparison.html",),
        timeout=15 * 60,
    )


def test_minimal_film_cooling_example(tmp_path: Path) -> None:
    run_example(
        "examples/film_cooling/minimal_film_cooling.py",
        tmp_path / "minimal-film-cooling",
        ("film_engine.html", "film_report.html"),
    )


def test_coupled_film_regen_example(tmp_path: Path) -> None:
    run_example(
        "examples/film_cooling/coupled_film_regen.py",
        tmp_path / "coupled-film-regen",
        ("coupled_film_report.html",),
    )


'''def test_advanced_simulation_and_post_processing(tmp_path: Path) -> None:
    """Post-process the result produced in the same isolated test run."""

    output_dir = tmp_path / "advanced"
    results_path = output_dir / "results.pkl"

    run_example(
        "examples/advanced/sizer_sim.py",
        output_dir,
        ("results.pkl",),
        timeout=60 * 60,
    )
    run_example(
        "examples/advanced/post_process.py",
        output_dir,
        ("methane_engine_report.html", "contour.html", "network.html"),
        extra_args=("--input", results_path),
    )'''
