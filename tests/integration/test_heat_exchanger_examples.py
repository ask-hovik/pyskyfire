"""Run the two public exchanger arrangements, checking budgets and reports."""

import importlib.util
from pathlib import Path
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]


def load_example(name):
    spec = importlib.util.spec_from_file_location(
        name, ROOT / "examples" / "heat_exchanger" / f"{name}.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_two_fluid_example_heats_oxygen_and_reports_both_wall_interfaces(tmp_path):
    example = load_example("two_fluid")
    result = example.main(tmp_path, nodes=8)
    assert result.circuits["oxygen"].T[-1] > result.circuits["oxygen"].T[0]
    assert result.circuits["hydrogen"].T[-1] < result.circuits["hydrogen"].T[0]
    assert abs(result.energy_balance_error) < 0.01
    text = (tmp_path / "two-fluid-report.html").read_text()
    assert "Wall temperatures by exposed surface" in text
    assert "opposite interface" in text
    assert "mesh3d" in text


def test_baffles_are_paired_v_rings_and_straight_spokes_with_series_flow(tmp_path):
    from pyskyfire.viz import heat_exchanger_report

    example = load_example("baffles")
    chamber, installation, inlets = example.build_baffles(nodes=8)
    exchanger = installation.exchanger
    inner = exchanger.circuits["inner ring down"]
    outer = exchanger.circuits["outer ring return"]
    assert inner.path.centerline[-1, 0] > 0 > outer.path.centerline[-1, 0]
    # Radius increases along flow on both paths: the return runs toward the injector.
    assert inner.path.centerline[-1, 1] > 0
    assert outer.path.centerline[-1, 1] > 0
    assert inner.cross_section.floor == inner.cross_section.roof == "arc"
    assert exchanger.circuits["spokes down"].cross_section.floor == "line"
    assert inner.n_channels == 48
    assert exchanger.circuits["spokes down"].n_channels == 30
    np.testing.assert_allclose(
        exchanger.circuits["spokes down"].path.centerline[:, 1:], 0
    )
    result = installation.solve(inlets)
    for edge in exchanger.flow_connections:
        a, b = result.circuits[edge.upstream], result.circuits[edge.downstream]
        assert a.outlet.T == pytest.approx(b.T[0])
        assert a.outlet.p == pytest.approx(b.p[0])
        assert a.mass_flow == b.mass_flow
    assert abs(result.energy_balance_error) < 0.1
    assert installation.volume_change.outside == pytest.approx(0, abs=1e-10)
    assert installation.volume_change.chamber_displaced > 0
    assert all(
        np.isfinite([w.temperature_a, w.temperature_b]).all() for w in result.walls
    )
    report = heat_exchanger_report(exchanger, result, installation=installation)
    report.save_html(tmp_path / "baffles.html")
    assert (
        "Wall temperatures by exposed surface"
        in (tmp_path / "baffles.html").read_text()
    )
