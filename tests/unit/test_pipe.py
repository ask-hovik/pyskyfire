"""Geometry and pressure-loss contracts for explicit pipes."""

from dataclasses import replace
import math
import warnings
from types import SimpleNamespace

import fluids
import numpy as np
import pytest

from pyskyfire import ducts, viz
from pyskyfire.common.friction import annulus_laminar_factor, darcy_friction_factor
from pyskyfire.ducts import pipe as pipe_module
from pyskyfire.thermo import FluidState, PHASE_UNKNOWN


class ConstantTransport:
    def __init__(self, rho=1000.0, mu=1.0e-3):
        self._rho = rho
        self._mu = mu

    def get_rho(self, T, p):
        return self._rho

    def get_mu(self, T, p):
        return self._mu

    def get_cp(self, T, p):
        return 4180.0


def inlet(p=2.0e6, T=300.0):
    return FluidState(
        T=T,
        p=p,
        h=np.nan,
        quality=np.nan,
        rho=np.nan,
        phase=PHASE_UNKNOWN,
        T_sat=np.nan,
    )


def test_generator_replaces_corner_vertex_with_forward_tangent_stations():
    nodes = [
        ducts.PathNode((0.0, 0.0, 0.0)),
        ducts.PathNode((0.0, 0.0, 0.5), ducts.RoundedBend(0.05)),
        ducts.PathNode((0.0, 0.5, 0.75)),
    ]

    stations, connections = ducts.generate_path_stations(nodes)
    bend = connections[1]
    incoming = np.array([0.0, 0.0, 1.0])
    outgoing = np.array([0.0, 0.5, 0.25])
    outgoing /= np.linalg.norm(outgoing)
    angle = math.acos(float(np.dot(incoming, outgoing)))
    tangent = 0.05 * math.tan(angle / 2.0)

    assert [type(connection) for connection in connections] == [
        ducts.Straight,
        ducts.RoundedBend,
        ducts.Straight,
    ]
    np.testing.assert_allclose(stations[1].position, np.array(nodes[1].position) - tangent * incoming)
    np.testing.assert_allclose(stations[2].position, np.array(nodes[1].position) + tangent * outgoing)
    assert bend.angle == pytest.approx(angle)
    assert np.dot(np.array(stations[1].position), incoming) > 0.0
    assert np.dot(np.array(nodes[-1].position) - np.array(stations[2].position), outgoing) > 0.0


def test_continuous_corner_is_one_exact_arc_connection():
    radius = 0.2
    nodes = [
        ducts.PathNode((0, 0, 0)),
        ducts.PathNode((1, 0, 0), ducts.Continuous(radius)),
        ducts.PathNode((1, 1, 0)),
    ]

    stations, connections = ducts.generate_path_stations(nodes)

    assert [type(item) for item in connections] == [
        ducts.Straight, ducts.Continuous, ducts.Straight]
    arc = connections[1]
    assert arc.radius == radius
    assert arc.length == pytest.approx(radius * math.pi / 2.0)
    center = np.array([0.8, 0.2, 0.0])
    np.testing.assert_allclose(arc.center, center, atol=1e-12)
    arc_points = np.asarray([station.position for station in stations[1:3]])
    np.testing.assert_allclose(np.linalg.norm(arc_points - center, axis=1), radius)


def test_constant_bore_arc_loss_needs_no_extra_stations():
    """One exact arc connection loses what the same arc split into many does."""
    def loss(path):
        pipe = ducts.Pipe(path, diameter=0.02, roughness=1e-5)
        return pipe.solve(inlet_state=inlet(), mass_flow=0.5,
                          transport=ConstantTransport()).total_pressure_loss

    path = ducts.Path.from_nodes([
        ducts.PathNode((0, 0, 0)),
        ducts.PathNode((1, 0, 0), ducts.Continuous(0.2)),
        ducts.PathNode((1, 1, 0)),
    ])
    assert len(path.stations) == 4
    assert loss(path) == pytest.approx(loss(path.refined(0.005)), rel=1e-9)


@pytest.mark.parametrize(
    "nodes, message",
    [
        ([ducts.PathNode((0, 0, 0), ducts.RoundedBend(0.1)), ducts.PathNode((1, 0, 0))], "first and last"),
        ([ducts.PathNode((0, 0, 0)), ducts.PathNode((1, 0, 0)), ducts.PathNode((1, 1, 0))], "must define"),
        ([ducts.PathNode((0, 0, 0)), ducts.PathNode((1, 0, 0), ducts.RoundedBend(0.1)), ducts.PathNode((0, 0, 0))], "reversing"),
        ([ducts.PathNode((0, 0, 0)), ducts.PathNode((0.1, 0, 0), ducts.RoundedBend(1.0)), ducts.PathNode((0.1, 1, 0))], "incoming leg"),
    ],
)
def test_generator_rejects_ambiguous_or_impossible_paths(nodes, message):
    with pytest.raises(ValueError, match=message):
        ducts.generate_path_stations(nodes)


def test_generator_rejects_overlapping_neighboring_bends():
    nodes = [
        ducts.PathNode((0, 0, 0)),
        ducts.PathNode((1, 0, 0), ducts.RoundedBend(0.6)),
        ducts.PathNode((1, 1, 0), ducts.RoundedBend(0.6)),
        ducts.PathNode((2, 1, 0)),
    ]
    with pytest.raises(ValueError, match="overlap or touch"):
        ducts.generate_path_stations(nodes)


def test_collinear_control_node_is_preserved_as_path_station():
    stations, connections = ducts.generate_path_stations([
        ducts.PathNode((0, 0, 0)),
        ducts.PathNode((1, 0, 0)),
        ducts.PathNode((2, 0, 0)),
    ])
    assert [station.position[0] for station in stations] == [0, 1, 2]
    assert all(isinstance(connection, ducts.Straight) for connection in connections)


def test_pipe_accepts_manual_public_representation_and_exposes_visual_geometry():
    stations = (
        ducts.PathStation((0, 0, 0)),
        ducts.PathStation((1, 0, 0)),
    )
    pipe = ducts.Pipe(stations, (ducts.Straight(),), roughness=1.0e-5, diameter=0.02)
    np.testing.assert_allclose(pipe.centerline, [[0, 0, 0], [1, 0, 0]])
    np.testing.assert_allclose(pipe.radius, [0.01, 0.01])
    assert pipe.stations is stations


def test_pipe_validates_connection_count_and_lumped_bend_diameter():
    stations = [
        ducts.PathStation((0, 0, 0)),
        ducts.PathStation((1, 0, 0)),
    ]
    with pytest.raises(ValueError, match="exactly one connection"):
        ducts.Pipe(stations, (), diameter=[0.02, 0.03])
    with pytest.raises(ValueError, match="equal endpoint diameters"):
        ducts.Pipe(stations, (ducts.RoundedBend(0.1, angle=math.pi / 2),), diameter=[0.02, 0.03])


def test_straight_loss_matches_darcy_weisbach_for_constant_properties():
    diameter = 0.02
    length = 1.5
    mass_flow = 0.2
    rho = 1000.0
    mu = 1.0e-3
    roughness = 2.0e-5
    pipe = ducts.Pipe(
        [ducts.PathStation((0, 0, 0)), ducts.PathStation((length, 0, 0))],
        [ducts.Straight()],
        roughness=roughness,
    diameter=diameter)

    result = pipe.solve(
        inlet_state=inlet(),
        mass_flow=mass_flow,
        transport=ConstantTransport(rho, mu),
        friction_model="fluids",
    )
    velocity = mass_flow / (rho * math.pi * diameter**2 / 4.0)
    reynolds = rho * velocity * diameter / mu
    fd = ducts.darcy_friction_factor(reynolds, diameter, roughness, model="fluids")
    expected = fd * length / diameter * rho * velocity**2 / 2.0

    assert result.total_pressure_loss == pytest.approx(expected)
    assert result.pressure_loss[0] == pytest.approx(expected)
    assert result.darcy_factor[0] == pytest.approx(fd)
    assert result.loss_coefficient[0] == pytest.approx(fd * length / diameter)


def test_variable_diameter_changes_station_velocity_but_adds_no_acceleration_loss():
    pipe = ducts.Pipe(
        [ducts.PathStation((0, 0, 0)), ducts.PathStation((1, 0, 0))],
        [ducts.Straight()],
    diameter=[0.01, 0.02])
    result = pipe.solve(inlet_state=inlet(), mass_flow=0.05, transport=ConstantTransport())

    assert result.velocity[0] / result.velocity[-1] == pytest.approx(4.0)
    np.testing.assert_allclose(result.p_stagnation[1:],
                               result.p_stagnation[0] - np.cumsum(result.pressure_loss))
    assert len(result.states) == len(pipe.stations)
    assert np.all(result.T_stagnation == inlet().T)


def test_static_state_is_stagnation_less_local_dynamic_head():
    pipe = ducts.Pipe(
        [ducts.PathStation((0, 0, 0)), ducts.PathStation((1, 0, 0))],
        [ducts.Straight()],
    diameter=[0.01, 0.02])
    result = pipe.solve(inlet_state=inlet(), mass_flow=0.05, transport=ConstantTransport())

    assert result.p_stagnation[0] == inlet().p
    np.testing.assert_allclose(result.p_static,
                               result.p_stagnation - 0.5*result.density*result.velocity**2)
    np.testing.assert_allclose(result.T_static, inlet().T - result.velocity**2/(2*4180.0))
    # The narrow inlet carries sixteen times the outlet's dynamic pressure.
    assert result.p_stagnation[0] - result.p_static[0] > result.p_stagnation[1] - result.p_static[1]
    assert result.T_static[0] < result.T_static[1] < result.T_stagnation[1]


def test_continuous_connection_applies_explicit_curvature(monkeypatch):
    captured = []
    monkeypatch.setattr(pipe_module, "darcy_friction_factor", lambda *args, **kwargs: 0.02)

    def multiplier(reynolds, diameter, radius, **kwargs):
        captured.append(radius)
        return 3.0

    monkeypatch.setattr(pipe_module, "curvature_multiplier", multiplier)
    pipe = ducts.Pipe(
        [ducts.PathStation((0, 0, 0)), ducts.PathStation((0.1, 0.1, 0))],
        [ducts.Continuous(radius=0.1, length=math.pi * 0.1 / 2.0)],
    diameter=0.02)
    result = pipe.solve(inlet_state=inlet(), mass_flow=0.1, transport=ConstantTransport())

    assert captured and all(radius == 0.1 for radius in captured)
    assert result.darcy_factor[0] == pytest.approx(0.06)


def test_manual_continuous_connections_estimate_menger_curvature(monkeypatch):
    radii = []
    monkeypatch.setattr(pipe_module, "darcy_friction_factor", lambda *args, **kwargs: 0.02)
    monkeypatch.setattr(
        pipe_module,
        "curvature_multiplier",
        lambda reynolds, diameter, radius, **kwargs: radii.append(radius) or 1.0,
    )
    pipe = ducts.Pipe(
        [
            ducts.PathStation((1, 0, 0)),
            ducts.PathStation((0, 1, 0)),
            ducts.PathStation((-1, 0, 0)),
        ],
        [ducts.Continuous(), ducts.Continuous()],
    diameter=0.02)
    pipe.solve(inlet_state=inlet(), mass_flow=0.1, transport=ConstantTransport())
    assert radii
    np.testing.assert_allclose(radii, 1.0)


def test_rounded_bend_matches_complete_fluids_K_without_arc_double_counting():
    diameter = 0.02
    radius = 0.06
    angle = math.pi / 2.0
    mass_flow = 0.1
    rho = 1000.0
    mu = 1.0e-3
    pipe = ducts.Pipe(
        [ducts.PathStation((0, 0, 0)), ducts.PathStation((radius, radius, 0))],
        [ducts.RoundedBend(radius, angle=angle)],
        roughness=1.0e-5,
    diameter=diameter)
    result = pipe.solve(
        inlet_state=inlet(), mass_flow=mass_flow, transport=ConstantTransport(rho, mu)
    )
    velocity = mass_flow / (rho * math.pi * diameter**2 / 4.0)
    reynolds = rho * velocity * diameter / mu
    fd = ducts.darcy_friction_factor(reynolds, diameter, 1.0e-5, model="fluids")
    expected_K = fluids.fittings.bend_rounded(
        Di=diameter,
        angle=90.0,
        rc=radius,
        fd=fd,
        method="Rennels",
    )

    assert result.loss_coefficient[0] == pytest.approx(expected_K)
    assert result.total_pressure_loss == pytest.approx(expected_K * rho * velocity**2 / 2.0)
    assert np.isnan(result.darcy_factor[0])


def test_pipe_rejects_invalid_mass_flow_and_nonpositive_pressure():
    pipe = ducts.Pipe(
        [ducts.PathStation((0, 0, 0)), ducts.PathStation((100, 0, 0))],
        [ducts.Straight()],
    diameter=0.001)
    with pytest.raises(ValueError, match="mass_flow"):
        pipe.solve(inlet_state=inlet(), mass_flow=0.0, transport=ConstantTransport())
    with pytest.raises(RuntimeError, match="pressure fell"):
        pipe.solve(inlet_state=inlet(1.0e5), mass_flow=1.0, transport=ConstantTransport())


def test_tapered_pipe_adds_stations_until_each_connection_varies_little():
    stations = [ducts.PathStation((0, 0, 0)), ducts.PathStation((1, 0, 0))]
    pipe = ducts.Pipe(stations, [ducts.Straight()], diameter=[0.01, 0.02])

    change = pipe.diameter[1:] / pipe.diameter[:-1] - 1.0
    assert np.all(change > 0) and change.max() <= 0.02 + 1e-12
    # ln(2)/ln(1.02) pieces, equidistributed in ln D.
    assert len(pipe.connections) == math.ceil(math.log(2) / math.log1p(0.02))
    np.testing.assert_allclose(change, change[0], rtol=0.05)
    np.testing.assert_allclose(pipe.centerline[[0, -1]], [[0, 0, 0], [1, 0, 0]])
    assert pipe.path.s[-1] == pytest.approx(1.0)

    coarse = ducts.Pipe(stations, [ducts.Straight()], diameter=[0.01, 0.02],
                        max_diameter_change=0.5)
    assert len(coarse.connections) == 2


def test_tapered_loss_converges_with_the_diameter_tolerance():
    stations = [ducts.PathStation((0, 0, 0)), ducts.PathStation((1, 0, 0))]

    def loss(tolerance):
        pipe = ducts.Pipe(stations, [ducts.Straight()], diameter=[0.01, 0.02],
                          max_diameter_change=tolerance)
        return pipe.solve(inlet_state=inlet(), mass_flow=0.05,
                          transport=ConstantTransport()).total_pressure_loss

    assert loss(0.02) == pytest.approx(loss(0.001), rel=5e-4)
    assert loss(2.0) < 0.8 * loss(0.001)  # one connection: far too little


def test_cone_angle_is_signed_full_angle_along_the_centerline():
    path = ducts.Path.from_nodes([
        ducts.PathNode((0, 0, 0)),
        ducts.PathNode((1, 0, 0), ducts.Continuous(0.2)),
        ducts.PathNode((1, 1, 0)),
    ])
    widening = ducts.Pipe(path, diameter=lambda s: 0.02 + 0.02 * s / path.s[-1])
    narrowing = ducts.Pipe(path, diameter=lambda s: 0.04 - 0.02 * s / path.s[-1])

    expected = 2.0 * math.atan(0.01 / path.s[-1])
    np.testing.assert_allclose(widening.cone_angle, expected, rtol=1e-9)
    np.testing.assert_allclose(narrowing.cone_angle, -expected, rtol=1e-9)


def straight_pipe(s, diameter, **kwargs):
    """A straight pipe along x with one station at each ``s`` [m]."""
    stations = [ducts.PathStation((x, 0, 0)) for x in s]
    return ducts.Pipe(stations, [ducts.Straight()] * (len(s) - 1), diameter=diameter, **kwargs)


def area_change_K(pipe, reference):
    """Total area-change loss coefficient on the velocity head at ``reference`` diameter."""
    mean = 0.5 * (pipe.diameter[:-1] + pipe.diameter[1:])
    return float(np.sum(pipe.area_change_loss_coefficient * (reference / mean)**4))


@pytest.mark.parametrize("tolerance", [0.02, 0.002, math.inf])
def test_cone_loses_rennels_expansion_term_however_finely_it_is_split(tolerance):
    length = 0.0125 / math.tan(math.radians(5.0))  # 25 to 50 mm at 10 degrees
    pipe = straight_pipe([0.0, length], [0.025, 0.05], max_diameter_change=tolerance)
    expected = fluids.fittings.diffuser_conical(0.025, 0.05, l=length, fd=0.0)

    assert area_change_K(pipe, 0.025) == pytest.approx(expected, rel=1e-9)


def test_step_reaches_the_sudden_expansion_limit():
    pipe = straight_pipe([0.0, 1.0e-4], [0.025, 0.05])

    assert area_change_K(pipe, 0.025) == pytest.approx((1 - 0.25)**2, rel=0.02)


def test_solve_adds_the_expansion_loss_to_wall_friction():
    pipe = straight_pipe([0.0, 0.1], [0.025, 0.045])
    result = pipe.solve(inlet_state=inlet(), mass_flow=1.0, transport=ConstantTransport())
    mean = 0.5 * (pipe.diameter[:-1] + pipe.diameter[1:])
    friction = result.darcy_factor * np.diff(result.s) / mean

    assert np.all(pipe.area_change_loss_coefficient > 0)
    np.testing.assert_allclose(result.loss_coefficient - friction,
                               pipe.area_change_loss_coefficient, rtol=1e-12)


def test_nearly_parallel_tail_does_not_dilute_a_steep_cone():
    steep = 0.0125 / math.tan(math.radians(5.0))       # 25 to 50 mm at 10 degrees
    tail = 0.0005 / math.tan(math.radians(0.1))        # 50 to 51 mm at 0.2 degrees
    pipe = straight_pipe([0.0, steep, steep + tail], [0.025, 0.05, 0.051])
    cone = fluids.fittings.diffuser_conical(0.025, 0.05, l=steep, fd=0.0)
    mean_angle = fluids.fittings.diffuser_conical(0.025, 0.051, l=steep + tail, fd=0.0)

    assert area_change_K(pipe, 0.025) == pytest.approx(cone, rel=0.05)
    assert mean_angle < 0.2 * cone


def test_short_constant_stretch_joins_two_cones_and_a_long_one_separates_them():
    first, second = 0.05, 0.08                  # 25 to 35 mm, then 35 to 50 mm

    def pipe(gap):
        return straight_pipe([0.0, first, first + gap, first + gap + second],
                             [0.025, 0.035, 0.035, 0.05])

    joined = straight_pipe([0.0, first, first + second], [0.025, 0.035, 0.05])
    apart = (fluids.fittings.diffuser_conical(0.025, 0.035, l=first, fd=0.0)
             + fluids.fittings.diffuser_conical(0.035, 0.05, l=second, fd=0.0)
             * (0.025 / 0.035)**4)

    assert area_change_K(pipe(2 * 0.035), 0.025) == pytest.approx(area_change_K(joined, 0.025), rel=1e-12)
    assert area_change_K(pipe(10 * 0.035), 0.025) == pytest.approx(apart, rel=1e-9)
    assert apart < area_change_K(joined, 0.025)


@pytest.mark.parametrize("tolerance", [0.02, 0.002, math.inf])
def test_contraction_loses_rennels_contraction_term_however_finely_it_is_split(tolerance):
    length = 0.0125 / math.tan(math.radians(15.0))  # 50 to 25 mm at 30 degrees
    pipe = straight_pipe([0.0, length], [0.05, 0.025], max_diameter_change=tolerance)
    expected = fluids.fittings.contraction_conical(0.05, 0.025, l=length, fd=0.0)

    assert area_change_K(pipe, 0.025) == pytest.approx(expected, rel=1e-9)


def test_step_down_reaches_the_sudden_contraction_limit():
    pipe = straight_pipe([0.0, 1.0e-4], [0.05, 0.025])

    assert area_change_K(pipe, 0.025) == pytest.approx(
        fluids.fittings.contraction_sharp(0.05, 0.025), rel=0.02)


def test_narrowing_ends_a_diffuser():
    bulge = straight_pipe([0.0, 0.1, 0.11, 0.21], [0.025, 0.05, 0.025, 0.05])
    cone = fluids.fittings.diffuser_conical(0.025, 0.05, l=0.1, fd=0.0)
    neck = fluids.fittings.contraction_conical(0.05, 0.025, l=0.01, fd=0.0)

    assert area_change_K(bulge, 0.025) == pytest.approx(2 * cone + neck, rel=1e-9)


def test_bore_variation_inside_a_connection_is_resolved():
    """A bulge between equal end diameters still gets its stations."""
    def bore(s):
        return 0.02 * (1.0 + 0.5 * np.sin(np.pi * s))

    pipe = ducts.Pipe(
        [ducts.PathStation((0, 0, 0)), ducts.PathStation((1, 0, 0))],
        [ducts.Straight()],
        diameter=bore,
    )
    assert len(pipe.stations) > 2
    s = pipe.path.s
    for start, end in zip(s[:-1], s[1:]):
        inside = bore(np.linspace(start, end, 50))
        assert inside.max() / inside.min() <= 1.02 + 1e-3


def test_flaring_arc_is_split_along_the_exact_arc():
    path = ducts.Path.from_nodes([
        ducts.PathNode((0, 0, 0)),
        ducts.PathNode((1, 0, 0), ducts.Continuous(0.2)),
        ducts.PathNode((1, 1, 0)),
    ])
    pipe = ducts.Pipe(path, diameter=lambda s: 0.02 + 0.02 * s / path.s[-1])
    arc = [i for i, c in enumerate(pipe.connections) if isinstance(c, ducts.Continuous)]

    assert len(arc) > 1
    center = np.asarray(pipe.connections[arc[0]].center)
    points = pipe.centerline[arc[0]:arc[-1] + 2]
    np.testing.assert_allclose(np.linalg.norm(points - center, axis=1), 0.2, atol=1e-12)
    assert pipe.path.s[-1] == pytest.approx(path.s[-1])
    assert all(c.radius == 0.2 for c in pipe.connections if isinstance(c, ducts.Continuous))


def test_split_manual_continuous_keeps_the_curvature_of_its_original_stations(monkeypatch):
    radii = []
    monkeypatch.setattr(pipe_module, "darcy_friction_factor", lambda *args, **kwargs: 0.02)
    monkeypatch.setattr(
        pipe_module, "curvature_multiplier",
        lambda reynolds, diameter, radius, **kwargs: radii.append(radius) or 1.0,
    )
    pipe = ducts.Pipe(
        [ducts.PathStation((1, 0, 0)), ducts.PathStation((0, 1, 0)), ducts.PathStation((-1, 0, 0))],
        [ducts.Continuous(), ducts.Continuous()],
        diameter=[0.02, 0.025, 0.03],
    )
    pipe.solve(inlet_state=inlet(), mass_flow=0.1, transport=ConstantTransport())
    assert len(pipe.connections) > 2
    np.testing.assert_allclose(radii, 1.0)


def test_max_diameter_change_is_validated():
    stations = [ducts.PathStation((0, 0, 0)), ducts.PathStation((1, 0, 0))]
    for bad in (0.0, -0.1, math.nan):
        with pytest.raises(ValueError, match="max_diameter_change"):
            ducts.Pipe(stations, [ducts.Straight()], diameter=0.01, max_diameter_change=bad)
    kept = ducts.Pipe(stations, [ducts.Straight()], diameter=[0.01, 0.02],
                      max_diameter_change=math.inf)
    assert kept.stations == tuple(stations)


def bent_pipe(radius=0.05, diameter=0.01):
    stations, connections = ducts.generate_path_stations([
        ducts.PathNode((0, 0, 0)),
        ducts.PathNode((1, 0, 0), bend=ducts.RoundedBend(radius=radius)),
        ducts.PathNode((1, 1, 0)),
    ])
    return ducts.Pipe(stations, connections, diameter=diameter)


def test_display_path_reopens_lumped_bends_without_touching_the_solver_model():
    radius = 0.05
    pipe = bent_pipe(radius)
    bend = pipe.connections[1]
    path = pipe.display_path(segments=12)

    # Three connections in, and only the bend is refined.
    assert len(pipe.stations) == 4
    assert len(path.points) == 3 + 12
    assert path.station_count == 4 and path.connection_count == 3

    arc = path.points[1:][path.connection == 1]
    np.testing.assert_allclose(
        np.linalg.norm(arc - np.asarray(bend.center), axis=1), radius, atol=1e-12
    )
    # Endpoints and the unrefined stations are the solver's own.
    np.testing.assert_allclose(path.points[0], pipe.centerline[0])
    np.testing.assert_allclose(path.points[-1], pipe.centerline[-1])
    np.testing.assert_allclose(arc[-1], pipe.centerline[2])
    np.testing.assert_allclose(path.radius, pipe.radius[0])
    assert pipe.centerline.shape == (4, 3)


def test_display_path_opens_generated_continuous_arcs_per_quarter_turn():
    path = ducts.Path.from_nodes([
        ducts.PathNode((0, 0, 0)),
        ducts.PathNode((1, 0, 0), ducts.Continuous(0.2)),
        ducts.PathNode((1, 1, 0)),
    ])
    pipe = ducts.Pipe(path, diameter=0.01)
    display = pipe.display_path(segments=6)

    assert len(display.points) == 4 + 6 - 1
    arc = display.points[1:][display.connection == 1]
    np.testing.assert_allclose(np.linalg.norm(arc - [0.8, 0.2, 0], axis=1), 0.2, atol=1e-12)


def test_display_path_station_and_connection_indices_track_their_source():
    path = bent_pipe().display_path(segments=4)
    np.testing.assert_allclose(path.station, [0, 1, 1.25, 1.5, 1.75, 2, 3])
    np.testing.assert_array_equal(path.connection, [0, 1, 1, 1, 1, 2])
    assert np.all(np.diff(path.station) > 0)


def test_display_path_keeps_the_chord_when_a_bend_has_no_resolved_plane():
    pipe = ducts.Pipe(
        [ducts.PathStation((0, 0, 0)), ducts.PathStation((1, 0, 0))],
        [ducts.RoundedBend(0.1, angle=math.pi / 2)],
    diameter=0.01)
    path = pipe.display_path(segments=8)
    np.testing.assert_allclose(path.points, pipe.centerline)
    np.testing.assert_array_equal(path.connection, [0])


def test_display_path_rejects_a_bend_plane_that_contradicts_its_stations():
    stations, connections = ducts.generate_path_stations([
        ducts.PathNode((0, 0, 0)),
        ducts.PathNode((1, 0, 0), bend=ducts.RoundedBend(radius=0.05)),
        ducts.PathNode((1, 1, 0)),
    ])
    broken = replace(connections[1], center=(5.0, 5.0, 5.0))
    pipe = ducts.Pipe(stations, (connections[0], broken, connections[2]), diameter=0.01)
    with pytest.raises(ValueError, match="inconsistent with its tangent stations"):
        pipe.display_path()


def test_display_path_validates_segments_and_interpolates_changing_radius():
    pipe = bent_pipe()
    for bad in (0, 1.5, True):
        with pytest.raises(ValueError, match="segments must be a positive integer"):
            pipe.display_path(segments=bad)

    tapered = ducts.Pipe(
        [ducts.PathStation((0, 0, 0)), ducts.PathStation((1, 0, 0))],
        [ducts.Straight()],
    diameter=[0.02, 0.04])
    np.testing.assert_allclose(tapered.display_path().radius, tapered.radius)
    np.testing.assert_allclose(tapered.radius[[0, -1]], [0.01, 0.02])


def test_pipe_pressure_plots_follow_arc_length() -> None:
    pytest.importorskip("plotly")
    from pyskyfire.viz import PlotPipePressure, PlotPipePressureLoss

    result = SimpleNamespace(
        s=np.array([0.0, 0.5, 1.0]),
        p_stagnation=np.array([2.0e6, 1.99e6, 1.98e6]),
        p_static=np.array([1.99e6, 1.985e6, 1.975e6]),
    )

    pressure = PlotPipePressure(result, labels=["duct"])
    assert [trace.name for trace in pressure.fig.data] == [
        "Stagnation pressure", "Static pressure"]
    np.testing.assert_allclose(pressure.fig.data[0].x, result.s)
    np.testing.assert_allclose(pressure.fig.data[1].y, result.p_static)

    loss = PlotPipePressureLoss(result, result, labels=["a", "b"])
    assert [trace.name for trace in loss.fig.data] == ["a", "b"]
    np.testing.assert_allclose(loss.fig.data[0].y, [0.0, 1.0e4, 2.0e4])

    with pytest.raises(ValueError, match="one entry per result"):
        PlotPipePressureLoss(result, labels=["a", "b"])


def test_from_curve_recovers_a_circle_exactly():
    radius = 0.3

    def circle(t):
        return np.column_stack((radius * np.cos(t), radius * np.sin(t), np.zeros_like(t)))

    path = ducts.Path.from_curve(circle, (0.0, math.pi / 2), max_turn=math.radians(5))

    assert all(isinstance(c, ducts.Continuous) for c in path.connections)
    np.testing.assert_allclose([c.radius for c in path.connections], radius, rtol=1e-9)
    np.testing.assert_allclose([c.center for c in path.connections], 0.0, atol=1e-9)
    assert path.s[-1] == pytest.approx(radius * math.pi / 2, rel=1e-12)
    assert len(path.connections) == 32  # 90 degrees halved to at most 5 per connection
    # Arcs are exact, so drawing reopens them consistently.
    ducts.Pipe(path, diameter=0.01).display_path()


def test_from_curve_crowds_stations_into_tight_bends_and_keeps_straights_sparse():
    def hook(t):
        # Straight along x for t < 1, then a tight quarter circle of radius 0.05.
        t = np.asarray(t, dtype=float)
        angle = np.clip(t - 1.0, 0.0, 1.0) * math.pi / 2
        return np.column_stack((np.minimum(t, 1.0) + 0.05 * np.sin(angle),
                                0.05 * (1 - np.cos(angle)), np.zeros_like(t)))

    path = ducts.Path.from_curve(hook, (0.0, 2.0))
    straight = [c for c in path.connections if isinstance(c, ducts.Straight)]
    # The straight half keeps its initial intervals, plus one halving where it
    # meets the arc and the curvature jumps.
    assert len(straight) <= 9
    assert len(path.connections) - len(straight) >= 18  # 90 degrees at <= 5 each
    arc = [c for c in path.connections if isinstance(c, ducts.Continuous)]
    np.testing.assert_allclose([c.radius for c in arc], 0.05, rtol=1e-9)


def test_from_curve_spline_length_and_loss_converge():
    from scipy.interpolate import CubicSpline

    control = np.array([[0, 0, 0], [0, 0, 0.3], [0.15, 0.1, 0.5], [0.4, 0.25, 0.55]])
    t = np.concatenate(([0.0], np.cumsum(np.linalg.norm(np.diff(control, axis=0), axis=1))))
    spline = CubicSpline(t, control, bc_type="natural")
    dense = np.linspace(0.0, t[-1], 200001)
    length = np.trapezoid(np.linalg.norm(spline(dense, 1), axis=1), dense)

    def loss(path):
        return ducts.Pipe(path, diameter=0.025, roughness=1.5e-5).solve(
            inlet_state=inlet(), mass_flow=2.0, transport=ConstantTransport()
        ).total_pressure_loss

    default = ducts.Path.from_curve(spline, (0.0, t[-1]))
    fine = ducts.Path.from_curve(spline, (0.0, t[-1]), max_turn=math.radians(0.25))
    assert default.s[-1] == pytest.approx(length, rel=1e-5)
    assert len(fine.stations) > 10 * len(default.stations)
    assert loss(default) == pytest.approx(loss(fine), rel=1e-3)


def test_from_curve_accepts_scalar_and_planar_curves():
    def scalar_only(t):
        return [math.cos(t), math.sin(t)]  # one planar point per scalar

    path = ducts.Path.from_curve(scalar_only, (0.0, 1.0))
    np.testing.assert_allclose(path.centerline[:, 2], 0.0)
    np.testing.assert_allclose(path.centerline[-1], [math.cos(1), math.sin(1), 0], atol=1e-12)


def test_from_curve_rejects_corners_stalls_and_bad_input():
    def corner(t):
        t = np.asarray(t, dtype=float)
        return np.column_stack((np.minimum(t, 0.5), np.maximum(t - 0.5, 0.0), 0 * t))

    with pytest.raises(ValueError, match="corner"):
        ducts.Path.from_curve(corner, (0.0, 1.0))
    with pytest.raises(ValueError, match="does not advance"):
        ducts.Path.from_curve(lambda t: np.zeros((np.size(t), 3)), (0.0, 1.0))
    with pytest.raises(ValueError, match="shape"):
        ducts.Path.from_curve(lambda t: np.zeros((np.size(t), 4)), (0.0, 1.0))
    with pytest.raises(ValueError, match="t_span"):
        ducts.Path.from_curve(lambda t: t, (1.0, 1.0))
    with pytest.raises(ValueError, match="max_turn"):
        ducts.Path.from_curve(lambda t: t, (0.0, 1.0), max_turn=0.0)


def test_flaring_pipe_on_a_sampled_curve_splits_along_its_arcs():
    def arc(t):
        return np.column_stack((0.2 * np.sin(t), 0.2 * (1 - np.cos(t)), np.zeros_like(t)))

    path = ducts.Path.from_curve(arc, (0.0, math.pi / 2))
    pipe = ducts.Pipe(path, diameter=lambda s: 0.01 + 0.02 * s / path.s[-1])
    assert len(pipe.stations) > len(path.stations)
    np.testing.assert_allclose(np.linalg.norm(pipe.centerline - [0, 0.2, 0], axis=1), 0.2,
                               atol=1e-9)
    assert pipe.path.s[-1] == pytest.approx(path.s[-1], rel=1e-12)


# -- a core inside the bore -------------------------------------------------

def straight_path(length=1.0):
    return ducts.Path([ducts.PathStation((0, 0, 0)), ducts.PathStation((length, 0, 0))],
                      [ducts.Straight()])


def bent_path(radius=0.1):
    return ducts.Path.from_nodes([
        ducts.PathNode((0, 0, 0)),
        ducts.PathNode((0, 0, 0.5), ducts.RoundedBend(radius)),
        ducts.PathNode((0, 0.5, 0.5)),
    ])


def test_a_pipe_without_a_core_is_its_bore():
    pipe = ducts.Pipe(straight_path(), diameter=0.04)

    assert not pipe.has_core
    np.testing.assert_allclose(pipe.hydraulic_diameter, 0.04)
    np.testing.assert_allclose(pipe.area, math.pi/4*0.04**2)
    np.testing.assert_allclose(pipe.r_mean, 0.01)


def test_straight_annulus_matches_darcy_weisbach():
    pipe = ducts.Pipe(straight_path(), diameter=0.04, core_diameter=0.02)
    result = pipe.solve(inlet_state=inlet(), mass_flow=2.0, transport=ConstantTransport())
    area = math.pi/4*(0.04**2 - 0.02**2)
    velocity = 2.0/(1000.0*area)
    reynolds = 1000.0*velocity*0.02/1e-3
    expected = fluids.friction.Clamond(reynolds, 0.0)/0.02*1000.0*velocity**2/2

    assert pipe.has_core
    assert result.total_pressure_loss == pytest.approx(expected, rel=1e-9)
    np.testing.assert_allclose(result.area, area)


def test_laminar_annulus_uses_the_exact_annulus_solution():
    pipe = ducts.Pipe(straight_path(), diameter=0.04, core_diameter=0.02)
    result = pipe.solve(inlet_state=inlet(), mass_flow=0.01, transport=ConstantTransport(mu=0.1))
    velocity = 0.01/(1000.0*math.pi/4*(0.04**2 - 0.02**2))
    reynolds = 1000.0*velocity*0.02/0.1
    expected = 64*annulus_laminar_factor(0.5)/reynolds/0.02*1000.0*velocity**2/2

    assert result.total_pressure_loss == pytest.approx(expected, rel=1e-9)


def test_rounded_bend_round_a_core_uses_the_pipe_correlation_on_the_chosen_scale_and_warns():
    pipe = ducts.Pipe(bent_path(), diameter=0.04, core_diameter=0.02)
    kwargs = dict(inlet_state=inlet(), mass_flow=2.0, transport=ConstantTransport())
    with pytest.warns(ducts.PipeCorrelationWarning, match="circular-pipe correlation"):
        hydraulic = pipe.solve(**kwargs)
    with pytest.warns(ducts.PipeCorrelationWarning):
        outer = pipe.solve(bend_length_scale="outer", **kwargs)
    bend = pipe.connections[1]
    velocity = hydraulic.meridional_velocity[1]
    reynolds = 1000.0*velocity*0.02/1e-3
    fd = darcy_friction_factor(reynolds, 0.02, model="fluids")
    expected = fluids.fittings.bend_rounded(Di=0.02, angle=math.degrees(bend.angle), rc=bend.radius,
                                            fd=fd, method="Rennels")

    assert hydraulic.loss_coefficient[1] == pytest.approx(expected, rel=1e-9)
    assert np.isnan(hydraulic.darcy_factor[1])
    assert outer.loss_coefficient[1] < hydraulic.loss_coefficient[1]


def test_swirl_into_a_rounded_bend_warns_with_or_without_a_core():
    kwargs = dict(inlet_state=inlet(), mass_flow=2.0, transport=ConstantTransport(), inlet_swirl=5.0)
    with pytest.warns(ducts.PipeCorrelationWarning, match="without swirl"):
        ducts.Pipe(bent_path(), diameter=0.04).solve(**kwargs)
    with warnings.catch_warnings():
        warnings.simplefilter("error", ducts.PipeCorrelationWarning)
        ducts.Pipe(bent_path(), diameter=0.04).solve(**{**kwargs, "inlet_swirl": 0.0})


def test_core_validation():
    with pytest.raises(ValueError, match="smaller than diameter"):
        ducts.Pipe(straight_path(), diameter=0.04, core_diameter=0.04)
    with pytest.raises(ValueError, match="non-negative"):
        ducts.Pipe(straight_path(), diameter=0.04, core_diameter=-0.01)
    with pytest.raises(ValueError, match="bore radius"):
        ducts.Pipe(bent_path(radius=0.015), diameter=0.04, core_diameter=0.01)
    with pytest.raises(ValueError, match="equal endpoint diameters"):
        ducts.Pipe(bent_path(), diameter=0.04, core_diameter=lambda s: 0.01 + 0.01*s)
    with pytest.raises(ValueError, match="bend_length_scale"):
        ducts.Pipe(straight_path(), diameter=0.04, core_diameter=0.01).solve(
            inlet_state=inlet(), mass_flow=1.0, transport=ConstantTransport(), bend_length_scale="mean")


def test_swirl_lengthens_the_friction_path_and_decays_in_a_plain_pipe():
    pipe = ducts.Pipe(straight_path(), diameter=0.04)
    plain = pipe.solve(inlet_state=inlet(), mass_flow=2.0, transport=ConstantTransport())
    velocity = plain.meridional_velocity[0]
    swirling = pipe.solve(inlet_state=inlet(), mass_flow=2.0, transport=ConstantTransport(),
                          inlet_swirl=velocity)

    assert np.all(plain.swirl_velocity == 0.0)
    np.testing.assert_allclose(plain.velocity, plain.meridional_velocity)
    assert swirling.total_pressure_loss > 2.0*plain.total_pressure_loss
    assert swirling.swirl_velocity[0] == pytest.approx(velocity)
    assert np.all(np.diff(swirling.swirl_velocity) < 0.0)
    assert swirling.swirl_angle[0] == pytest.approx(np.pi/4)


def test_a_core_growing_from_a_point_is_resolved_finitely_and_contracts_the_flow():
    def nose(s):
        return 0.03*min(s, 0.1)/0.1    # a cone from a point to 30 mm

    pipe = ducts.Pipe(straight_path(0.2), diameter=0.04, core_diameter=nose)
    hydraulic = pipe.hydraulic_diameter
    equivalent = np.sqrt(4*pipe.area/math.pi)

    assert pipe.core_diameter[0] == 0.0 and len(pipe.stations) < 400
    for quantity in (hydraulic, equivalent):
        assert np.max(np.abs(np.diff(np.log(quantity)))) <= math.log1p(0.02) + 1e-3
    assert np.all(pipe.cone_angle <= 1e-12)
    assert np.sum(pipe.area_change_loss_coefficient) > 0
    result = pipe.solve(inlet_state=inlet(), mass_flow=1.0, transport=ConstantTransport())
    assert result.meridional_velocity[-1] > result.meridional_velocity[0]


def test_a_blunt_core_end_is_a_sudden_expansion():
    s = [0.0, 0.5, 0.5001, 1.0]
    pipe = straight_pipe(s, 0.04, core_diameter=[0.02, 0.02, 0.0, 0.0])
    area = pipe.area
    mid_D = 0.5*(pipe.diameter[:-1] + pipe.diameter[1:])
    mid_d = 0.5*(pipe.core_diameter[:-1] + pipe.core_diameter[1:])
    mid = math.pi/4*(mid_D**2 - mid_d**2)
    K = float(np.sum(pipe.area_change_loss_coefficient*(area[0]/mid)**2))

    assert K == pytest.approx((1 - area[0]/area[-1])**2, rel=0.02)


def test_component_draws_the_core_inside_a_translucent_bore():
    pipe = ducts.Pipe(bent_path(), diameter=0.04, core_diameter=0.02)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ducts.PipeCorrelationWarning)
        result = pipe.solve(inlet_state=inlet(), mass_flow=2.0, transport=ConstantTransport())
    node = viz.PipeComponent(pipe, result, n_theta=12).build()
    core, bore = node.children

    assert [c.name for c in node.children] == ["core", "bore"]
    assert bore.opacity == pytest.approx(0.35) and core.opacity == 1.0
    assert set(node.anchors) == {"inlet", "outlet"}
    for wall in (core, bore):
        assert len(wall.fields["p_static"].values) == len(wall.mesh.points)
        assert len(wall.fields["pressure_loss"].values) == len(wall.mesh.triangles)
        assert "swirl_velocity" in wall.fields

    plain = viz.PipeComponent(ducts.Pipe(bent_path(), diameter=0.04)).build()
    assert plain.mesh is not None and not plain.children


def test_component_draws_a_core_only_where_it_runs():
    pipe = straight_pipe([0.0, 0.5, 0.5001, 1.0], 0.04, core_diameter=[0.02, 0.02, 0.0, 0.0])
    result = pipe.solve(inlet_state=inlet(), mass_flow=1.0, transport=ConstantTransport())
    node = viz.PipeComponent(pipe, result, n_theta=8).build()
    core, bore = node.children
    reach = np.ptp(core.mesh.points[:, 0])

    assert reach == pytest.approx(0.5001, abs=1e-3)
    assert len(core.fields["p_static"].values) == len(core.mesh.points)
    assert len(core.fields["reynolds"].values) == len(core.mesh.triangles)
    assert len(core.mesh.points) < len(bore.mesh.points)

