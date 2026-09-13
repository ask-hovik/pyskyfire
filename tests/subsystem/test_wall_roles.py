"""Wall roles through geometry, conduction, reporting and serialization."""

from types import SimpleNamespace

import cloudpickle
import numpy as np
import pytest

import pyskyfire as psf
from pyskyfire.regen.coupled_solver import CoupledHeatExchangerPhysics


class _Material:
    def get_k(self, temperature):
        return 350.0


class _Contour:
    xs = np.array([0.0, 1.0])
    x_t = 0.0

    def r(self, x):
        return np.full_like(np.asarray(x, dtype=float), 0.1)

    def dr_dx(self, x):
        return np.zeros_like(np.asarray(x, dtype=float))


def _wall(thickness, **kwargs):
    return psf.regen.Wall(_Material(), thickness, **kwargs)


def _circuit(walls, section=None, **kwargs):
    return psf.regen.CoolingCircuit(
        name="wall roles", contour=_Contour(),
        cross_section=section or psf.regen.CrossSectionRounded(height_input="base"),
        placement=psf.regen.SurfacePlacement(60),
        channel_height=0.004, walls=walls, span=[0.0, 1.0],
        coolant_transport=None, roughness=0.0,
        **kwargs,
    )


def _chamber(circuit):
    return psf.regen.ThrustChamber(
        contour=circuit.contour, cooling_circuits=[circuit],
        combustion_transport=None, n_nodes=11, compute_gas=False,
    )


def test_default_and_legacy_walls_are_main_even_when_named_liner():
    wall = _wall(0.001, name="liner")
    assert wall.role == "main"
    state = wall.__dict__.copy()
    for field in ("role", "_evaluator", "_contour"):
        state.pop(field)
    restored = psf.regen.Wall.__new__(psf.regen.Wall)
    restored.__setstate__(state)
    assert restored.role == "main"
    assert restored.thickness(0.0) == pytest.approx(0.001)
    liner = _wall(0.002, role="liner")
    assert cloudpickle.loads(cloudpickle.dumps(liner)).role == "liner"


def test_wall_role_and_hot_to_cold_order_are_validated():
    with pytest.raises(ValueError, match="Wall.role"):
        _wall(0.001, role="lining")
    with pytest.raises(ValueError, match="Liner walls must precede"):
        _circuit([_wall(0.001), _wall(0.002, role="liner")])
    with pytest.raises(ValueError, match="at least one main wall"):
        _circuit([_wall(0.002, role="liner")])


@pytest.mark.parametrize("x", [0.5, np.array([0.0, 0.5, 1.0]), np.array([[0.0, 1.0]])])
def test_thickness_sums_use_roles_and_preserve_input_shape(x):
    circuit = _circuit([
        _wall(lambda x: 0.001 + 0.001 * x, role="liner"),
        _wall(0.0002, role="liner"),
        _wall(0.0003),
        _wall(lambda x: 0.0005 + 0.0001 * x),
    ])
    main = 0.0008 + 0.0001 * np.asarray(x)
    liner = 0.0012 + 0.001 * np.asarray(x)
    for actual, expected in [
        (circuit.channel_wall_thickness(x), main),
        (circuit.hot_side_layer_thickness(x), liner),
        (circuit.total_thickness(x), main + liner),
    ]:
        assert np.shape(actual) == np.shape(x)
        np.testing.assert_allclose(actual, expected)


def test_unmarked_wall_stubs_are_main_and_absent_liners_return_shaped_zero():
    circuit = _circuit([
        SimpleNamespace(thickness=lambda x: 0.0003),
        _wall(0.0005),
    ])
    x = np.array([0.0, 0.5, 1.0])
    np.testing.assert_allclose(circuit.channel_wall_thickness(x), [0.0008] * 3)
    np.testing.assert_array_equal(circuit.hot_side_layer_thickness(x), np.zeros(3))
    assert circuit.hot_side_layer_thickness(0.5) == 0.0


@pytest.mark.parametrize("section", [
    psf.regen.CrossSectionSquared(0.001),
    psf.regen.CrossSectionRounded(height_input="base"),
])
def test_multiple_main_layers_match_one_combined_wall_in_simulation_and_cad(section):
    split = _circuit([
        _wall(0.002, role="liner"), _wall(0.0003), _wall(0.0005),
    ], section)
    combined = _circuit([_wall(0.002, role="liner"), _wall(0.0008)], section)
    _chamber(split)
    _chamber(combined)

    # The reference path is the ideal contour; every layer is applied when the
    # section is placed, so how the main stack is split cannot matter.
    np.testing.assert_allclose(split.reference_path[:, 1], 0.1)
    np.testing.assert_allclose(split.cross_section_base_point_path,
                               combined.cross_section_base_point_path)
    prof = split._prof(split.cross_section_base_point_path)
    np.testing.assert_allclose(prof.t_wall, 0.0008)
    np.testing.assert_allclose(prof.t_hot_layer, 0.002)
    np.testing.assert_allclose(prof.t_wall_total, 0.0028)
    np.testing.assert_allclose(split.section_profiles_at(0.5).t_wall, 0.0008)
    for attribute in ["A_coolant_vals", "Dh_coolant_vals", "cooling_channel_centroid_path"]:
        np.testing.assert_allclose(getattr(split, attribute), getattr(combined, attribute),
                                   atol=1e-15)
    split_cad = psf.regen.build_channel_construction(split)
    combined_cad = psf.regen.build_channel_construction(combined)
    for i in range(len(split.x_domain)):
        np.testing.assert_allclose(split_cad.polygon(i), combined_cad.polygon(i))


def test_liners_remain_in_series_conduction_and_interface_temperatures():
    circuit = _circuit([_wall(0.002, role="liner"), _wall(0.0003), _wall(0.0005)])
    chamber = _chamber(circuit)
    physics = CoupledHeatExchangerPhysics.__new__(CoupledHeatExchangerPhysics)
    physics.thrust_chamber = chamber
    physics.circuit_index = 0
    area = circuit.dA_dx_thermal_exhaust(0.5)
    assert physics.dQ_cond_dx(0.5, 600.0, 300.0) == pytest.approx(
        300.0 * 350.0 * area / 0.0028
    )
    np.testing.assert_allclose(
        physics.interface_temperatures(0.5, 600.0, 300.0),
        600.0 - 300.0 * np.array([0.0, 0.002, 0.0023, 0.0028]) / 0.0028,
    )


def test_wall_plot_repeats_all_main_layers_in_closeout_but_no_liners():
    chamber = _chamber(_circuit([
        _wall(0.002, role="liner", name="silver"),
        _wall(0.0003, name="outer tube"),
        _wall(0.0005, name="inner tube"),
    ]))
    plot = psf.viz.PlotWallLayers(chamber, mirror=False, show_contour=False)
    closeouts = [trace for trace in plot.fig.data if "closeout" in trace.name]
    assert [trace.name for trace in closeouts] == [
        "inner tube (closeout, assumed)", "outer tube (closeout, assumed)",
    ]
    assert max(closeouts[1].y) - min(closeouts[0].y) == pytest.approx(0.0008)


class _RibMaterial:
    def __init__(self):
        self.temperatures = []

    def get_k(self, temperature):
        self.temperatures.append(temperature)
        return 10.0 * temperature


@pytest.mark.parametrize("section", [
    psf.regen.CrossSectionSquared(0.001),
    psf.regen.CrossSectionRounded(height_input="base"),
])
def test_circuits_keep_independent_fin_settings_and_materials(section):
    walls = [_wall(0.002, role="liner"), _wall(0.0008)]
    material = _RibMaterial()
    disabled_material = _RibMaterial()
    disabled = _circuit(walls, section, enable_fin=False, rib_material=disabled_material)
    default = _circuit(walls, section)
    explicit = _circuit(walls, section, rib_material=material)
    assert default.enable_fin is True
    assert default.rib_material is None
    chamber = psf.regen.ThrustChamber(
        contour=default.contour, cooling_circuits=[disabled, default, explicit],
        combustion_transport=None, n_nodes=11, compute_gas=False,
    )
    assert [c.enable_fin for c in chamber.cooling_circuits] == [False, True, True]
    assert not hasattr(chamber, "enable_fin")

    resistances = [
        c.R_coolant_per_len(0.5, h_c=1000.0, T_wall_rep=500.0)
        for c in chamber.cooling_circuits
    ]
    assert resistances[0] == pytest.approx(
        1.0 / (1000.0 * disabled.dA_dx_thermal_coolant(0.5))
    )
    assert resistances[2] < resistances[1] < resistances[0]
    assert material.temperatures == [500.0]
    assert disabled_material.temperatures == []
    for circuit in [default, explicit]:
        np.testing.assert_allclose(circuit.A_coolant_vals, disabled.A_coolant_vals)

    # The explicit rib material does not replace any conduction-layer material.
    physics = CoupledHeatExchangerPhysics.__new__(CoupledHeatExchangerPhysics)
    physics.thrust_chamber = chamber
    conduction = []
    for i in range(3):
        physics.circuit_index = i
        conduction.append(physics.dQ_cond_dx(0.5, 600.0, 300.0))
    np.testing.assert_allclose(conduction, conduction[0])
    assert material.temperatures == [500.0]


def test_fin_settings_survive_serialization_and_legacy_circuits_keep_disabled_state():
    circuit = _circuit([_wall(0.0008)], enable_fin=False, rib_material=_RibMaterial())
    restored = cloudpickle.loads(cloudpickle.dumps(circuit))
    assert restored.enable_fin is False
    assert restored.rib_material.get_k(500.0) == 5000.0

    state = circuit.__dict__.copy()
    state.pop("rib_material")
    legacy = psf.regen.CoolingCircuit.__new__(psf.regen.CoolingCircuit)
    legacy.__setstate__(state)
    assert legacy.enable_fin is False
    assert legacy.rib_material is None
