"""Contracts for lazy cooling-channel section planes."""

import numpy as np

from pyskyfire.regen.channel_helpers import build_section_planes
from pyskyfire.regen.thrust_chamber import CoolingCircuit


def test_axial_centerline_produces_radial_section_frames() -> None:
    centerline = np.array(
        [
            [0.0, 2.0, 0.0],
            [1.0, 2.0, 0.0],
            [2.0, 2.0, 0.0],
        ]
    )

    planes = build_section_planes(centerline)

    np.testing.assert_allclose(
        planes.origins,
        [[0.0, 2.0, 0.0], [1.0, 2.0, 0.0], [2.0, 2.0, 0.0]],
    )
    np.testing.assert_allclose(planes.normals, [[1.0, 0.0, 0.0]] * 3)
    np.testing.assert_allclose(planes.x_axes, [[0.0, 1.0, 0.0]] * 3)
    np.testing.assert_allclose(planes.y_axes, [[0.0, 0.0, 1.0]] * 3)


def test_helical_section_frames_are_orthonormal() -> None:
    axial = np.linspace(0.0, 1.0, 11)
    centerline = np.column_stack((axial, np.full_like(axial, 0.2), axial * 0.4))

    planes = build_section_planes(centerline)

    np.testing.assert_allclose(np.linalg.norm(planes.normals, axis=1), 1.0)
    np.testing.assert_allclose(np.linalg.norm(planes.x_axes, axis=1), 1.0)
    np.testing.assert_allclose(np.linalg.norm(planes.y_axes, axis=1), 1.0)
    np.testing.assert_allclose(
        np.sum(planes.normals * planes.x_axes, axis=1),
        0.0,
        atol=1.0e-14,
    )


class _Cylinder:
    """Constant-radius contour, so the section offsets are purely radial."""

    def __init__(self, radius=0.1, length=0.2):
        self.radius = float(radius)
        self.xs = np.array([0.0, length])
        self.x_t = 0.0

    def r(self, x):
        return np.full_like(np.asarray(x, dtype=float), self.radius)

    def dr_dx(self, x):
        return np.zeros_like(np.asarray(x, dtype=float))


class _Wall:
    def __init__(self, thickness):
        self._thickness = float(thickness)
        self.role = "main"

    def thickness(self, x):
        return np.full_like(np.asarray(x, dtype=float), self._thickness)


def _built_circuit(*, radius=0.1, thickness=0.001, height=0.004, helix_angle=0.0):
    from pyskyfire.regen.channel_placement import SurfacePlacement
    from pyskyfire.regen.cross_section import CrossSectionSquared
    from pyskyfire.regen.thrust_chamber import ThrustChamber

    contour = _Cylinder(radius=radius)
    circuit = CoolingCircuit(
        name="planes",
        contour=contour,
        cross_section=CrossSectionSquared(0.001),
        span=[0.0, 1.0],
        placement=SurfacePlacement(n_channel_positions=60, helix_angle=helix_angle),
        channel_height=lambda x: height,
        walls=[_Wall(thickness)],
        coolant_transport=None,
        roughness=0.0,
    )
    ThrustChamber(
        contour=contour,
        cooling_circuits=[circuit],
        combustion_transport=None,
        n_nodes=12,
        compute_gas=False,
    )
    return circuit


def test_circuit_planes_are_anchored_on_the_cross_section_base_point() -> None:
    radius, thickness, height = 0.1, 0.001, 0.004
    circuit = _built_circuit(radius=radius, thickness=thickness, height=height)

    planes = circuit.section_planes()

    # The base point sits on the coolant-facing wall: outside the contour the
    # reference path follows, and inside the centroid the coolant travels on.
    np.testing.assert_allclose(np.linalg.norm(planes.origins[:, 1:], axis=1),
                               radius + thickness)
    np.testing.assert_allclose(circuit.reference_path[:, 1], radius)
    centroid_radius = circuit.cooling_channel_centroid_path[:, 1]
    assert np.all(centroid_radius > radius + thickness)
    assert np.all(centroid_radius < radius + thickness + height)
    # Close to mid-height, but measured from the section rather than assumed:
    # the curved floor and roof pull the area centroid slightly inward.
    np.testing.assert_allclose(centroid_radius, radius + thickness + height / 2, rtol=1e-3)
    assert np.all(centroid_radius < radius + thickness + height / 2)


def test_circuit_plane_calculation_is_not_stored_on_circuit() -> None:
    circuit = _built_circuit()
    before = set(circuit.__dict__)

    circuit.section_planes()

    assert set(circuit.__dict__) == before


def test_circuit_axis_normal_planes_are_transverse_to_the_engine_axis() -> None:
    from pyskyfire.regen.channel_helpers import build_axis_normal_section_planes

    circuit = _built_circuit(helix_angle=0.4)
    planes = circuit.section_planes(axis_normal=True)
    expected = build_axis_normal_section_planes(circuit.cross_section_base_point_path)

    np.testing.assert_allclose(planes.normals, [[1.0, 0.0, 0.0]] * len(planes.normals))
    np.testing.assert_allclose(planes.x_axes, expected.x_axes)
    np.testing.assert_allclose(planes.origins, expected.origins)
    # A helical circuit still has flow-normal planes that lean out of transverse.
    assert not np.allclose(circuit.section_planes().normals, planes.normals)


def test_axis_normal_frames_are_radial_even_on_a_sloping_helix() -> None:
    from pyskyfire.regen.channel_helpers import build_axis_normal_section_planes

    centerline = np.array([[0.0, 0.1, 0.0], [0.1, 0.15, 0.5], [0.2, 0.2, 1.0]])
    planes = build_axis_normal_section_planes(centerline)
    np.testing.assert_allclose(planes.normals, [[1, 0, 0]] * 3)
    # The chamber axis is exactly (-r, 0) in the local plane.
    axis_points = planes.origins - centerline[:, 1, None] * planes.x_axes
    np.testing.assert_allclose(axis_points[:, 1:], 0.0, atol=1e-16)
    np.testing.assert_allclose(np.cross(planes.x_axes, planes.y_axes), planes.normals)
