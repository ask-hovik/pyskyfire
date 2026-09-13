"""Geometry and pressure-loss contracts for explicit pipes."""

from dataclasses import replace
import math

import fluids
import numpy as np
import pytest

from pyskyfire import ducts
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
        ducts.PathNode((0.0, 0.0, 0.0), 0.01),
        ducts.PathNode((0.0, 0.0, 0.5), 0.01, ducts.RoundedBend(0.05)),
        ducts.PathNode((0.0, 0.5, 0.75), 0.01),
    ]

    stations, connections = ducts.generate_pipe_stations(nodes)
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


def test_continuous_corner_is_discretized_with_exact_arc_metadata():
    radius = 0.2
    segments = 8
    nodes = [
        ducts.PathNode((0, 0, 0), 0.02),
        ducts.PathNode((1, 0, 0), 0.02, ducts.Continuous(radius, segments)),
        ducts.PathNode((1, 1, 0), 0.02),
    ]

    stations, connections = ducts.generate_pipe_stations(nodes)
    continuous = [item for item in connections if isinstance(item, ducts.Continuous)]

    assert len(continuous) == segments
    assert all(item.radius == radius for item in continuous)
    assert sum(item.length for item in continuous) == pytest.approx(radius * math.pi / 2.0)
    center = np.array([0.8, 0.2, 0.0])
    arc_points = np.asarray([station.position for station in stations[1 : segments + 2]])
    np.testing.assert_allclose(np.linalg.norm(arc_points - center, axis=1), radius)


@pytest.mark.parametrize(
    "nodes, message",
    [
        ([ducts.PathNode((0, 0, 0), 0.01, ducts.RoundedBend(0.1)), ducts.PathNode((1, 0, 0), 0.01)], "first and last"),
        ([ducts.PathNode((0, 0, 0), 0.01), ducts.PathNode((1, 0, 0), 0.01), ducts.PathNode((1, 1, 0), 0.01)], "must define"),
        ([ducts.PathNode((0, 0, 0), 0.01), ducts.PathNode((1, 0, 0), 0.01, ducts.RoundedBend(0.1)), ducts.PathNode((0, 0, 0), 0.01)], "reversing"),
        ([ducts.PathNode((0, 0, 0), 0.01), ducts.PathNode((0.1, 0, 0), 0.01, ducts.RoundedBend(1.0)), ducts.PathNode((0.1, 1, 0), 0.01)], "incoming leg"),
    ],
)
def test_generator_rejects_ambiguous_or_impossible_paths(nodes, message):
    with pytest.raises(ValueError, match=message):
        ducts.generate_pipe_stations(nodes)


def test_generator_rejects_overlapping_neighboring_bends():
    nodes = [
        ducts.PathNode((0, 0, 0), 0.01),
        ducts.PathNode((1, 0, 0), 0.01, ducts.RoundedBend(0.6)),
        ducts.PathNode((1, 1, 0), 0.01, ducts.RoundedBend(0.6)),
        ducts.PathNode((2, 1, 0), 0.01),
    ]
    with pytest.raises(ValueError, match="overlap or touch"):
        ducts.generate_pipe_stations(nodes)


def test_collinear_control_node_is_preserved_as_diameter_station():
    stations, connections = ducts.generate_pipe_stations([
        ducts.PathNode((0, 0, 0), 0.01),
        ducts.PathNode((1, 0, 0), 0.02),
        ducts.PathNode((2, 0, 0), 0.03),
    ])
    assert [station.diameter for station in stations] == [0.01, 0.02, 0.03]
    assert all(isinstance(connection, ducts.Straight) for connection in connections)


def test_pipe_accepts_manual_public_representation_and_exposes_visual_geometry():
    stations = (
        ducts.PipeStation((0, 0, 0), 0.02),
        ducts.PipeStation((1, 0, 0), 0.04),
    )
    pipe = ducts.Pipe(stations, (ducts.Straight(),), roughness=1.0e-5)
    np.testing.assert_allclose(pipe.centerline, [[0, 0, 0], [1, 0, 0]])
    np.testing.assert_allclose(pipe.radius, [0.01, 0.02])
    assert pipe.stations is stations


def test_pipe_validates_connection_count_and_lumped_bend_diameter():
    stations = [
        ducts.PipeStation((0, 0, 0), 0.02),
        ducts.PipeStation((1, 0, 0), 0.03),
    ]
    with pytest.raises(ValueError, match="exactly one connection"):
        ducts.Pipe(stations, ())
    with pytest.raises(ValueError, match="equal endpoint diameters"):
        ducts.Pipe(stations, (ducts.RoundedBend(0.1, angle=math.pi / 2),))


def test_straight_loss_matches_darcy_weisbach_for_constant_properties():
    diameter = 0.02
    length = 1.5
    mass_flow = 0.2
    rho = 1000.0
    mu = 1.0e-3
    roughness = 2.0e-5
    pipe = ducts.Pipe(
        [ducts.PipeStation((0, 0, 0), diameter), ducts.PipeStation((length, 0, 0), diameter)],
        [ducts.Straight()],
        roughness=roughness,
    )

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
        [ducts.PipeStation((0, 0, 0), 0.01), ducts.PipeStation((1, 0, 0), 0.02)],
        [ducts.Straight()],
    )
    result = pipe.solve(inlet_state=inlet(), mass_flow=0.05, transport=ConstantTransport())

    assert result.velocity[0] / result.velocity[1] == pytest.approx(4.0)
    assert result.p[1] == pytest.approx(result.p[0] - result.pressure_loss[0])
    assert len(result.states) == 2
    assert np.all(result.T == inlet().T)


def test_continuous_connection_applies_explicit_curvature(monkeypatch):
    captured = []
    monkeypatch.setattr(pipe_module, "darcy_friction_factor", lambda *args, **kwargs: 0.02)

    def multiplier(reynolds, diameter, radius, **kwargs):
        captured.append(radius)
        return 3.0

    monkeypatch.setattr(pipe_module, "curvature_multiplier", multiplier)
    pipe = ducts.Pipe(
        [ducts.PipeStation((0, 0, 0), 0.02), ducts.PipeStation((0.1, 0.1, 0), 0.02)],
        [ducts.Continuous(radius=0.1, length=math.pi * 0.1 / 2.0)],
    )
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
            ducts.PipeStation((1, 0, 0), 0.02),
            ducts.PipeStation((0, 1, 0), 0.02),
            ducts.PipeStation((-1, 0, 0), 0.02),
        ],
        [ducts.Continuous(), ducts.Continuous()],
    )
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
        [ducts.PipeStation((0, 0, 0), diameter), ducts.PipeStation((radius, radius, 0), diameter)],
        [ducts.RoundedBend(radius, angle=angle)],
        roughness=1.0e-5,
    )
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
        [ducts.PipeStation((0, 0, 0), 0.001), ducts.PipeStation((100, 0, 0), 0.001)],
        [ducts.Straight()],
    )
    with pytest.raises(ValueError, match="mass_flow"):
        pipe.solve(inlet_state=inlet(), mass_flow=0.0, transport=ConstantTransport())
    with pytest.raises(RuntimeError, match="pressure fell"):
        pipe.solve(inlet_state=inlet(1.0e5), mass_flow=1.0, transport=ConstantTransport())


def bent_pipe(radius=0.05, diameter=0.01):
    stations, connections = ducts.generate_pipe_stations([
        ducts.PathNode((0, 0, 0), diameter),
        ducts.PathNode((1, 0, 0), diameter, bend=ducts.RoundedBend(radius=radius)),
        ducts.PathNode((1, 1, 0), diameter),
    ])
    return ducts.Pipe(stations, connections)


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


def test_display_path_station_and_connection_indices_track_their_source():
    path = bent_pipe().display_path(segments=4)
    np.testing.assert_allclose(path.station, [0, 1, 1.25, 1.5, 1.75, 2, 3])
    np.testing.assert_array_equal(path.connection, [0, 1, 1, 1, 1, 2])
    assert np.all(np.diff(path.station) > 0)


def test_display_path_keeps_the_chord_when_a_bend_has_no_resolved_plane():
    pipe = ducts.Pipe(
        [ducts.PipeStation((0, 0, 0), 0.01), ducts.PipeStation((1, 0, 0), 0.01)],
        [ducts.RoundedBend(0.1, angle=math.pi / 2)],
    )
    path = pipe.display_path(segments=8)
    np.testing.assert_allclose(path.points, pipe.centerline)
    np.testing.assert_array_equal(path.connection, [0])


def test_display_path_rejects_a_bend_plane_that_contradicts_its_stations():
    stations, connections = ducts.generate_pipe_stations([
        ducts.PathNode((0, 0, 0), 0.01),
        ducts.PathNode((1, 0, 0), 0.01, bend=ducts.RoundedBend(radius=0.05)),
        ducts.PathNode((1, 1, 0), 0.01),
    ])
    broken = replace(connections[1], center=(5.0, 5.0, 5.0))
    pipe = ducts.Pipe(stations, (connections[0], broken, connections[2]))
    with pytest.raises(ValueError, match="inconsistent with its tangent stations"):
        pipe.display_path()


def test_display_path_validates_segments_and_interpolates_changing_radius():
    pipe = bent_pipe()
    for bad in (0, 1.5, True):
        with pytest.raises(ValueError, match="segments must be a positive integer"):
            pipe.display_path(segments=bad)

    tapered = ducts.Pipe(
        [ducts.PipeStation((0, 0, 0), 0.02), ducts.PipeStation((1, 0, 0), 0.04)],
        [ducts.Straight()],
    )
    np.testing.assert_allclose(tapered.display_path().radius, [0.01, 0.02])
