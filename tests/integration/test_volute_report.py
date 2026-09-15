"""Exercise the user-facing volute example end to end."""
import math
import os
from pathlib import Path
import subprocess
import sys

import pytest


def test_volute_example_creates_a_scene_and_a_report(tmp_path):
    pytest.importorskip("plotly")
    pytest.importorskip("pyvista")
    script = Path(__file__).resolve().parents[2]/"examples/volute/volute_report.py"
    result = subprocess.run(
        [sys.executable, str(script), "--output-dir", str(tmp_path)],
        cwd=tmp_path, capture_output=True, text=True, timeout=300,
        env={**os.environ, "MPLBACKEND":"Agg", "MPLCONFIGDIR":str(tmp_path/"mpl")},
    )
    assert result.returncode == 0, result.stdout+result.stderr
    assert (tmp_path/"volute-report.html").stat().st_size > 1000
    # One standalone scene per volute, rather than a single shared one.
    scenes = sorted(tmp_path.glob("volute-scene-*.html"))
    assert len(scenes) == 3
    assert all(path.stat().st_size > 1000 for path in scenes)
    assert not (tmp_path/"volute-scene.html").exists()
    report = (tmp_path/"volute-report.html").read_text()
    # One embedded scene per volute, not a single combined one.
    assert report.count("<iframe") == 3
    assert "Coolant inlet manifold" in report
    assert "Sizing parameter SP" in report


def test_the_example_designs_two_radial_volutes_and_one_inclined():
    """The three cases are distinct duties, and all realise their design SP."""
    pytest.importorskip("plotly")
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from examples.volute.volute_report import build_cases

    cases = build_cases()
    assert [round(math.degrees(case.volute.tilt)) for case in cases] == [0, 0, 45]
    assert len({case.volute.inlet.alpha for case in cases}) == 3
    for case in cases:
        assert case.result.SP == pytest.approx(1.0, rel=2e-4)
        assert case.result.omega_tangential == pytest.approx(0.0, abs=1e-6)
        assert case.result.pressure_loss > 0.0
    # The manifold's inlet is the steepest, so it needs the largest passage.
    assert cases[2].result.A_6 > cases[1].result.A_6 > cases[0].result.A_6


def test_the_example_gives_each_case_a_different_outlet_shape():
    """Tangential, radial, and an intermediate angle with elevation."""
    pytest.importorskip("plotly")
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from examples.volute.volute_report import build_cases

    lox, methane, manifold = build_cases()
    assert methane.volute.exit_geometry.is_straight
    assert not lox.volute.exit_geometry.is_straight
    assert not manifold.volute.exit_geometry.is_straight

    # LOX asked for "radial", which resolves to exactly the radial bound.
    assert lox.volute.exit_geometry.first_angle == pytest.approx(
        lox.volute.exit_geometry.radial_angle
    )
    assert lox.volute.exit_geometry.second_angle == pytest.approx(0.0)

    # The manifold is the only one that also lifts out of its own plane.
    assert manifold.volute.exit_geometry.second_angle > 0.0
    assert lox.volute.exit_geometry.second_angle == pytest.approx(0.0)
    assert methane.volute.exit_geometry.second_angle == pytest.approx(0.0)

    # Each duct reaches its own flange diameter exactly.
    for case in (lox, methane, manifold):
        assert case.volute.exit_geometry.final_diameter == pytest.approx(
            case.exit_diameter
        )
