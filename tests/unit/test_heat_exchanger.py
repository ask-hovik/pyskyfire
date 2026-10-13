"""Conservation, counterflow, shared geometry, and prescribed-gas contracts."""

from types import SimpleNamespace
import numpy as np
import pytest
from pyskyfire import common, ducts, regen, thermo, heat_exchanger as hx
from pyskyfire.heat_exchanger import solver
from pyskyfire.common.physics import (
    gnielinski_temperature_correction,
    h_internal,
    nusselt_gnielinski,
)


class ConstantTransport:
    fluid = "manufactured-fluid"

    def get_mu(self, T, p):
        return 0.001

    def get_k(self, T, p):
        return 0.6

    def get_Pr(self, T, p):
        return 7.0


class ConstantResolver:
    def __init__(self, *args):
        self.capability = SimpleNamespace(enthalpy_march=True)

    def enthalpy(self, T, p):
        return 4000 * T

    def from_hp(self, h, p):
        return thermo.FluidState(
            h / 4000, p, h, np.nan, 1000.0, thermo.PHASE_LIQUID, np.nan
        )


@pytest.fixture
def constant_fluid(monkeypatch):
    monkeypatch.setattr(solver, "FluidStateResolver", ConstantResolver)


def circuit(
    name="a",
    length=1.0,
    count=40,
    copies=1,
    transport=None,
    area_width=0.01,
    exposure=None,
):
    path = ducts.Path.from_nodes(
        [ducts.PathNode((0, 0, 0)), ducts.PathNode((length, 0, 0))]
    ).refined(length / count)
    wall = regen.Wall(common.Material("constant", k=common.ConstantModel(20.0)), 0.001)
    section = regen.CrossSectionSquared(
        0.001, floor="line", roof="line", exposure=exposure
    )
    return hx.HeatExchangerCircuit(
        name=name,
        path=path,
        cross_section=section,
        wall=wall,
        transport=transport or ConstantTransport(),
        dimensions=dict(height=0.01, floor_width=area_width, roof_width=area_width),
        placements=[ducts.PathPlacement((0, i * 0.02, 0)) for i in range(copies)],
        heat_transfer=lambda **kw: 1000.0,
    )


def test_gnielinski_reference_and_laminar_branch():
    assert nusselt_gnielinski(1e5, 1.2, 0.0185) == pytest.approx(254.62682749359632)
    assert h_internal(100, 7, 0.6, 0.01, 0.64) == pytest.approx(3.66 * 0.6 / 0.01)
    with pytest.raises(ValueError, match="Gnielinski"):
        nusselt_gnielinski(100, 1, 0.02)


def test_gnielinski_temperature_correction_heats_only():
    assert gnielinski_temperature_correction(200.0, 800.0) == pytest.approx(0.25**0.45)
    assert gnielinski_temperature_correction(800.0, 200.0) == 1.0
    bulk = h_internal(1e5, 1.2, 0.1, 0.01, 0.0185)
    hot = h_internal(1e5, 1.2, 0.1, 0.01, 0.0185, bulk_temperature=200, wall_temperature=800)
    assert hot == pytest.approx(bulk * 0.25**0.45)
    # Laminar Nu is left alone.
    assert h_internal(
        100, 7, 0.6, 0.01, 0.64, bulk_temperature=200, wall_temperature=800
    ) == pytest.approx(3.66 * 0.6 / 0.01)
    with pytest.raises(ValueError, match="both"):
        h_internal(1e5, 1.2, 0.1, 0.01, 0.0185, wall_temperature=800)


def test_default_coefficient_uses_converged_wetted_wall_temperature(constant_fluid):
    c = circuit(count=4)
    c.heat_transfer = None
    result = hx.HeatExchanger([c]).solve(
        {"a": hx.Inlet(300, 2e6, 0.1)}, boundaries={"a": hx.PrescribedGas(900, 5000)}
    )
    r = result.circuits["a"]
    walls = solver._wetted_wall_temperatures("a", len(c.s) - 1, result.walls)
    bulk = solver.midpoint(r.T)
    assert np.all(walls > bulk)
    expected = [
        h_internal(Re, 7.0, 0.6, Dh, fd, bulk_temperature=Tb, wall_temperature=Tw)
        for Re, Dh, fd, Tb, Tw in zip(
            r.reynolds, solver.midpoint(c.diameter), r.darcy_factor, bulk, walls
        )
    ]
    assert r.heat_transfer_coefficient == pytest.approx(expected, rel=1e-6)


def test_bulk_properties_can_be_requested(constant_fluid):
    c = circuit(count=4)
    c.heat_transfer, c.property_correction = None, False
    r = (
        hx.HeatExchanger([c])
        .solve(
            {"a": hx.Inlet(300, 2e6, 0.1)},
            boundaries={"a": hx.PrescribedGas(900, 5000)},
        )
        .circuits["a"]
    )
    expected = [
        h_internal(Re, 7.0, 0.6, Dh, fd)
        for Re, Dh, fd in zip(r.reynolds, solver.midpoint(c.diameter), r.darcy_factor)
    ]
    assert r.heat_transfer_coefficient == pytest.approx(expected, rel=1e-6)


def test_pipe_and_exchanger_share_exact_path_and_placement():
    c = circuit()
    p = ducts.Pipe(c.path, diameter=0.01)
    assert p.path is c.path
    assert p.stations is c.path.stations
    assert c.displaced_volume == pytest.approx(0.012**2)
    points = c.section_points(samples=4)
    assert np.ptp(points[:, :, 1]) == pytest.approx(0.012)
    assert np.ptp(points[:, :, 2]) == pytest.approx(0.012)


def test_transported_frame_straight_and_generated_arc():
    path = ducts.Path.from_nodes(
        [
            ducts.PathNode((0, 0, 0)),
            ducts.PathNode((1, 0, 0), bend=ducts.Continuous(0.2)),
            ducts.PathNode((1, 1, 0)),
        ]
    )
    p = ducts.PathPlacement(
        (2, 3, 4), ducts.rotation_about_axis((1, 0, 0), 0.7), (0, 0, 1)
    )
    frames = p.frames(path)
    np.testing.assert_allclose(
        frames.transpose(0, 2, 1) @ frames,
        np.broadcast_to(np.eye(3), frames.shape),
        atol=1e-12,
    )
    np.testing.assert_allclose(frames[:, :, 2], p.place(path).tangents, atol=1e-12)
    np.testing.assert_allclose(p.place(path).centerline[0], [2, 3, 4])
    np.testing.assert_allclose(p.place(path).lengths, path.lengths)
    np.testing.assert_allclose(path.refined(0.01).s[-1], path.s[-1])
    with pytest.raises(ValueError, match="parallel"):
        ducts.PathPlacement(initial_floor_direction=(1, 0, 0)).frames(path)
    with pytest.raises(ValueError, match="orthonormal"):
        ducts.PathPlacement(rotation=-np.eye(3))


def test_adiabatic_march_conserves_energy_and_darcy_loss(constant_fluid):
    c = circuit(count=4)
    result = hx.HeatExchanger([c]).solve({"a": hx.Inlet(300, 2e6, 0.1)})
    r = result.circuits["a"]
    np.testing.assert_allclose(r.T, 300)
    assert r.energy_change == pytest.approx(0, abs=1e-8)
    expected = (
        common.darcy_friction_factor(10000, 0.01, 0, model="fluids")
        * 1
        / 0.01
        * 1000
        / 2
    )
    assert r.p[0] - r.p[-1] == pytest.approx(expected, rel=1e-6)


@pytest.mark.parametrize("counterflow", [False, True])
def test_two_fluid_matches_effectiveness_and_conserves_heat(
    constant_fluid, counterflow
):
    a, b = circuit("cold"), circuit("hot")
    n = len(a.s) - 1
    pairs = tuple((i, n - 1 - i if counterflow else i) for i in range(n))
    exchanger = hx.HeatExchanger(
        [a, b], thermal_links=[hx.ThermalLink("cold", "hot", pairs, a.wall)]
    )
    result = exchanger.solve(
        {"cold": hx.Inlet(300, 2e6, 0.01), "hot": hx.Inlet(500, 2e6, 0.01)}
    )
    UA = 0.01 / (1 / 1000 + 0.001 / 20 + 1 / 1000)
    NTU = UA / (0.01 * 4000)
    effectiveness = NTU / (1 + NTU) if counterflow else (1 - np.exp(-2 * NTU)) / 2
    assert result.circuits["cold"].total_heat == pytest.approx(
        effectiveness * 0.01 * 4000 * 200, rel=3e-4
    )
    assert abs(result.energy_balance_error) < 0.002
    for w in result.walls:
        assert 300 < w.temperature_a < w.temperature_b < 500
    assert result.circuits["cold"].total_heat == -result.circuits["hot"].total_heat


def test_equal_split_and_series_flow(constant_fluid):
    a, b = circuit("down", count=6, copies=3), circuit("up", count=6, copies=3)
    hxobj = hx.HeatExchanger([a, b], flow_connections=[hx.FlowConnection("down", "up")])
    result = hxobj.solve(
        {"down": hx.Inlet(300, 2e6, 0.3)},
        boundaries={name: hx.PrescribedGas(500, 1000) for name in ("down", "up")},
    )
    assert result.circuits["down"].velocity[0] == pytest.approx(1.0)
    assert result.circuits["up"].T[0] == pytest.approx(result.circuits["down"].T[-1])
    assert result.circuits["up"].p[0] == pytest.approx(result.circuits["down"].p[-1])
    assert abs(result.energy_balance_error) < 0.01
    assert result.circuits["up"].T[-1] > result.circuits["down"].T[-1] > 300


def test_thermal_area_budget_and_boundary_equivalence(constant_fluid):
    a, b = circuit("a", count=2, copies=2), circuit("b", count=2, copies=2)
    with pytest.raises(ValueError, match="area"):
        hx.HeatExchanger(
            [a, b], thermal_links=[hx.ThermalLink("a", "b", ((0, 0), (0, 0)), a.wall)]
        )
    with pytest.raises(ValueError, match="different gas conditions"):
        hx.HeatExchanger([a]).solve(
            {"a": hx.Inlet(300, 2e6, 0.1)},
            boundaries={"a": hx.PrescribedGas(lambda p: 500 + p[1] * 100, 1000)},
        )
    with pytest.raises(ValueError, match="indices"):
        hx.HeatExchanger(
            [a, b], thermal_links=[hx.ThermalLink("a", "b", ((-1, 0),), a.wall)]
        )


def test_nonconvergence_is_explicit(constant_fluid):
    c = circuit(count=2)
    with pytest.raises(RuntimeError, match="did not converge"):
        hx.HeatExchanger([c]).solve(
            {"a": hx.Inlet(300, 2e6, 0.1)},
            boundaries={"a": hx.PrescribedGas(500, 1000)},
            max_iterations=1,
        )


def test_coolprop_helium_heating_and_acceleration():
    tr = thermo.FluidTransport(
        thermo.Fluid(type="fuel", propellants=["Helium"], fractions=[1.0])
    )
    c = circuit(count=4, transport=tr)
    r = (
        hx.HeatExchanger([c])
        .solve(
            {"a": hx.Inlet(200, 5e6, 0.01)},
            boundaries={"a": hx.PrescribedGas(500, 500)},
        )
        .circuits["a"]
    )
    assert r.T[-1] > r.T[0]
    assert np.all(r.acceleration_drop > 0)
    assert r.p[0] - r.p[-1] == pytest.approx(
        sum(r.friction_loss + r.acceleration_drop), rel=1e-5
    )
    assert r.energy_change == pytest.approx(r.total_heat, rel=1e-6)


def test_chamber_volume_clips_envelope_at_injector():
    c = circuit(length=0.2, count=4)
    contour = regen.Contour([-0.1, -0.05, 0.0, 0.05], [0.1, 0.1, 0.05, 0.07])
    chamber = SimpleNamespace(contour=contour)
    install = hx.ChamberInstallation(
        chamber, hx.HeatExchanger([c]), position=(-0.2, 0, 0), coefficient=1000
    )
    assert install.volume_change.displaced == pytest.approx(
        0.5 * c.displaced_volume, rel=1e-3
    )
    assert install.volume_change.outside == pytest.approx(
        0.5 * c.displaced_volume, rel=1e-3
    )
    assert install.volume_change.chamber_displaced == pytest.approx(
        install.volume_change.displaced
    )


def test_partial_exposure_scales_heat_and_thick_roof_changes_wall_temperature(
    constant_fluid,
):
    full = circuit("full", count=8)
    partial = circuit("partial", count=8, exposure={"floor": 0.5})
    inlet = hx.Inlet(300, 2e6, 0.1)
    results = hx.HeatExchanger([full, partial]).solve(
        {"full": inlet, "partial": inlet},
        boundaries={n: hx.PrescribedGas(500, 1000) for n in ("full", "partial")},
    )
    ratio = results.circuits["partial"].total_heat / results.circuits["full"].total_heat
    assert 0.5 < ratio < 0.51  # half area, with the colder partial-flow bulk state
    c = circuit(exposure={"floor": 1.0, "roof": 1.0}, count=8)
    section = regen.CrossSectionSquared(
        0.001,
        roof_thickness=0.01,
        floor="line",
        roof="line",
        exposure={"floor": 1.0, "roof": 1.0},
    )
    thick = hx.HeatExchangerCircuit(
        name="thick",
        path=c.path,
        cross_section=section,
        wall=c.wall,
        transport=c.transport,
        dimensions=dict(height=0.01, floor_width=0.01, roof_width=0.01),
        heat_transfer=lambda **kw: 1000.0,
    )
    result = hx.HeatExchanger([thick]).solve(
        {"thick": inlet}, boundaries={"thick": hx.PrescribedGas(500, 1000)}
    )
    walls = {w.surface_a: w for w in result.walls if w.segment_a == 0}
    assert (
        walls["roof"].temperature_b - walls["roof"].temperature_a
        > walls["floor"].temperature_b - walls["floor"].temperature_a
    )


def test_two_phase_heating_is_rejected():
    tr = thermo.FluidTransport(
        thermo.Fluid(type="fuel", propellants=["Water"], fractions=[1.0])
    )
    c = circuit(count=2, transport=tr)
    with pytest.raises(ValueError, match="two-phase"):
        hx.HeatExchanger([c]).solve(
            {"a": hx.Inlet(400, 1e6, 0.001)},
            boundaries={"a": hx.PrescribedGas(700, 1000)},
        )


def test_chamber_adapter_uses_transformed_axial_position_and_shared_bartz():
    calls = []

    def state(x, T=None):
        calls.append((x, T))
        return SimpleNamespace(
            T=3000.0 if T is None else T,
            gamma=1.2,
            M=0.1,
            Pr=0.7,
            k=0.2,
            cp=3000.0,
            mu=1e-4,
        )

    chamber = SimpleNamespace(
        contour=regen.Contour([-0.2, -0.1, 0.0, 0.1], [0.15, 0.15, 0.08, 0.1]),
        combustion_transport=SimpleNamespace(get_state=state, mdot=10.0),
    )
    gas = hx.ChamberGas(chamber, position=(-0.15, 0, 0))
    T, h = gas.sample(np.array([0.02, 0, 0]), 500.0)
    assert T > 3000 and h > 0
    assert calls[0][0] == pytest.approx(-0.13)
    assert calls[1][1] == pytest.approx(1750.0)


def test_fin_conductivity_follows_wall_temperature(constant_fluid, monkeypatch):
    c = circuit(count=4)
    c.cross_section.enable_fin = True
    c.wall.material = common.Material(
        "temperature dependent", k=common.PolynomialModel([0.0, 0.01])
    )
    observed = []
    original = c.cross_section.R_coolant_per_len

    def resistance(prof, h, k):
        observed.append(np.asarray(k).copy())
        return original(prof, h, k)

    monkeypatch.setattr(c.cross_section, "R_coolant_per_len", resistance)
    result = hx.HeatExchanger([c]).solve(
        {"a": hx.Inlet(300, 2e6, 0.1)}, boundaries={"a": hx.PrescribedGas(500, 1000)}
    )
    assert np.min(observed[-1]) > 3.1
    assert abs(result.energy_balance_error) < 0.01


def _linear_segment(J, c):
    J, c = np.asarray(J, float), np.asarray(c, float)
    return lambda x: (J @ np.asarray(x, float) + c, None)


def test_segment_outlet_converges_in_a_few_evaluations_on_a_contractive_map():
    J, c = [[0.05, 0.02], [-0.03, 0.04]], [1.0, 2.0]
    x, _, evaluations = solver._segment_outlet(_linear_segment(J, c), (0.0, 0.0), np.ones(2))
    np.testing.assert_allclose(x, np.linalg.solve(np.eye(2) - np.asarray(J), c), rtol=1e-9)
    # The half-step iteration this replaces needed ~30 here.
    assert evaluations <= 6


def test_segment_outlet_falls_back_to_half_steps_where_a_trial_cannot_be_evaluated():
    # Substitution overshoots by 150 % and lands where the map refuses to go;
    # the half step does not.
    def G(x):
        if x[0] > 1.3:
            raise ValueError("outside the fluid's range")
        return np.array([2.5 - 1.5 * x[0], 0.0]), None

    x, _, _ = solver._segment_outlet(G, (0.0, 0.0), np.ones(2))
    np.testing.assert_allclose(x, (1.0, 0.0), atol=1e-9)


def test_segment_outlet_reports_no_fixed_point():
    x, _, evaluations = solver._segment_outlet(
        lambda x: (np.asarray(x, float) + 1.0, None), (0.0, 0.0), np.ones(2))
    assert x is None
    assert evaluations >= solver.SEGMENT_EVALUATIONS
