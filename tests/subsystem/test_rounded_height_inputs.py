"""Equivalent height conventions through chamber assembly and reporting."""

import numpy as np
import pytest

import pyskyfire as psf


class _Contour:
    xs = np.array([0.0, 1.0])
    x_t = 0.0

    def r(self, x):
        return 0.1 + 0.02 * np.asarray(x)

    def dr_dx(self, x):
        return np.full_like(np.asarray(x, dtype=float), 0.02)


class _Material:
    def get_k(self, temperature):
        return 350.0


def _base(x):
    return 0.002 + 0.001 * np.asarray(x)


def _gamma(x):
    return 0.3 + 0.1 * np.asarray(x)


def _expected_full_height(circuit):
    """Full internal height of the section the construction actually places.

    The lobes are sized where they end up, so the inner radius follows from the
    solved lobe-centre distance rather than from the contour plus a liner.
    """
    prof = psf.regen.build_channel_construction(circuit).profiles
    half_theta = prof.theta / 2
    minimum = 2 * (prof.sector_radius * np.sin(half_theta) - prof.t_wall)
    return minimum + _base(circuit.x_domain) * (1 + np.sin(half_theta))


def _chamber(mode, section_type=psf.regen.CrossSectionRounded, height=None):
    contour = _Contour()
    circuit = psf.regen.CoolingCircuit(
        name="Rounded",
        contour=contour,
        cross_section=section_type(height_input=mode),
        placement=psf.regen.SurfacePlacement(60, helix_angle=_gamma),
        channel_height=_base if height is None else height,
        walls=[
            psf.regen.Wall(_Material(), 0.0003, role="liner"),
            psf.regen.Wall(_Material(), 0.001),
        ],
        span=[0.0, 1.0], coolant_transport=None, roughness=0.0,
    )
    chamber = psf.regen.ThrustChamber(
        contour=contour, cooling_circuits=[circuit], combustion_transport=None,
        n_nodes=21, compute_gas=False,
    )
    return chamber, circuit


def _equivalent_pair(section_type=psf.regen.CrossSectionRounded):
    """Return the same channel described as a base height and as a full one."""
    _, base = _chamber("base", section_type)
    heights = _expected_full_height(base)
    full_chamber, full = _chamber(
        "full",
        section_type,
        height=lambda x: float(np.interp(x, base.x_domain, heights)),
    )
    return base, full_chamber, full


@pytest.mark.parametrize("section_type", [psf.regen.CrossSectionRounded, psf.regen.CrossSectionRoundedInternal])
def test_equivalent_inputs_produce_identical_solver_paths_and_construction(section_type):
    base, _, full = _equivalent_pair(section_type)
    np.testing.assert_allclose(base.channel_full_heights, _expected_full_height(base))
    np.testing.assert_allclose(base.channel_full_heights, full.channel_full_heights)
    np.testing.assert_allclose(base.channel_base_heights, full.channel_base_heights)
    np.testing.assert_allclose(base.channel_heights, _base(base.x_domain))
    assert not np.allclose(base.channel_heights, base.channel_full_heights)
    for attr in ("reference_path", "cooling_channel_centroid_path", "A_coolant_vals", "Dh_coolant_vals", "ds_dx_vals", "signed_curvature_vals"):
        np.testing.assert_allclose(getattr(base, attr), getattr(full, attr), rtol=1e-11, atol=1e-13)
    assert base.volume == pytest.approx(full.volume, rel=1e-13)
    base_geom = psf.regen.build_channel_construction(base)
    full_geom = psf.regen.build_channel_construction(full)
    for index in range(len(base.x_domain)):
        np.testing.assert_allclose(base_geom.polygon(index), full_geom.polygon(index), atol=1e-15)


def test_reports_default_to_full_height_and_can_select_base_or_both():
    chamber, circuit = _chamber("base")
    full_plot = psf.viz.PlotChannelHeight(chamber)
    base_plot = psf.viz.PlotChannelHeight(chamber, height="base")
    both = psf.viz.PlotChannelHeight(chamber, circuit_index=0, height="both")
    np.testing.assert_allclose(full_plot.fig.data[0].y, _expected_full_height(circuit))
    np.testing.assert_allclose(base_plot.fig.data[0].y, _base(circuit.x_domain))
    assert [trace.name for trace in both.fig.data] == ["Rounded (full)", "Rounded (base)"]
    assert both.fig.data[1].line.dash == "dash"
    np.testing.assert_allclose(both.fig.data[0].y, full_plot.fig.data[0].y)
    np.testing.assert_allclose(both.fig.data[1].y, base_plot.fig.data[0].y)
    report = psf.viz.Report(title="Channel heights")
    report.add_tab("Geometry").add_figure(both)


def test_wall_layer_plot_uses_full_height_for_base_input():
    base_chamber, _ = _chamber("base")
    _, full_chamber, _ = _equivalent_pair()
    with pytest.warns(RuntimeWarning, match="helical"):
        base_plot = psf.viz.PlotWallLayers(base_chamber, interpolation_points=21)
    with pytest.warns(RuntimeWarning, match="helical"):
        full_plot = psf.viz.PlotWallLayers(full_chamber, interpolation_points=21)
    for base_trace, full_trace in zip(base_plot.fig.data, full_plot.fig.data, strict=True):
        np.testing.assert_allclose(base_trace.y, full_trace.y, rtol=1e-7, atol=1e-8)
