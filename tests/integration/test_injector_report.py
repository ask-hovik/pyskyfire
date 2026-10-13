"""Exercise the user-facing injector example end to end."""
import os
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]


def test_injector_example_creates_scene_and_report(tmp_path):
    pytest.importorskip("plotly")
    pytest.importorskip("pyvista")
    script = ROOT/"examples/injector/injector_report.py"
    result = subprocess.run(
        [sys.executable, str(script), "--output-dir", str(tmp_path)],
        cwd=tmp_path, capture_output=True, text=True, timeout=300,
        env={**os.environ, "MPLBACKEND": "Agg", "MPLCONFIGDIR": str(tmp_path/"mpl")},
    )
    assert result.returncode == 0, result.stdout + result.stderr
    for name in ("injector-scene.html", "injector-report.html"):
        assert (tmp_path/name).stat().st_size > 1000
    report = (tmp_path/"injector-report.html").read_text()
    assert "<iframe" in report
    assert "Rupe ratio" in report and "Film ring" in report
    assert "Plenum for design flow" in report


def test_the_drilled_injector_meets_its_design_where_it_should():
    pytest.importorskip("plotly")
    sys.path.insert(0, str(ROOT))
    from examples.injector.injector_report import (DESIGN_FLOW, FILM, PLENUM_RATIO, P_C,
                                                   at_flows, at_pressures, build_design,
                                                   film_sweep)

    design = build_design()
    # Undrilled holes pass exactly the design flows at the design pressures.
    exact = at_pressures(design.designed)
    for feed, flow in DESIGN_FLOW.items():
        assert exact.feeds[feed].mass_flow == pytest.approx(flow, rel=1e-8)
    # The inclinations send the doublet's resultant down the axis.
    assert exact.sets["Main doublets"].diagnostics.resultant_angle == pytest.approx(0.0, abs=1e-9)

    # Drilled holes need their own plenum pressures to pass the same flows.
    needed = at_flows(design.drilled)
    for feed in DESIGN_FLOW:
        drilled, designed = design.diameters[feed][1], design.diameters[feed][0]
        moved = needed.feeds[feed].plenum_pressure/P_C - PLENUM_RATIO[feed]
        assert (moved > 0) == (drilled < designed)

    # A higher film ring pressure always carries more of the fuel.
    shares = film_sweep(design.drilled, [1.1, 1.3, 1.5])
    assert shares[0] < shares[1] < shares[2]
    assert design.drilled.checks(chamber_radius=0.065).overlaps == 0
    assert FILM in design.drilled.feeds
