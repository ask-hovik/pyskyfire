"""Physical section projection, separate from the one-dimensional model."""

from types import SimpleNamespace

import numpy as np
import pytest

from pyskyfire.regen import (
    CrossSectionRounded,
    CrossSectionSquared,
    build_channel_construction,
)


def circuit_geometry(*, slope=0.0, gamma=0.0, section=None, hot_layer=0.0, n_stations=9):
    x = np.linspace(0.0, 0.04, n_stations)
    # The reference path follows the ideal contour; wall layers are added by
    # the construction, not baked into the path.
    radius = 0.1 + slope * x
    # Wind the path at exactly the rate the stored gamma law asks for, the way
    # ThrustChamber.build_reference_paths integrates it.
    rate = np.tan(gamma) * np.sqrt(1 + slope**2)
    angle = (rate * np.log(radius / radius[0]) / slope if slope
             else rate * x / radius[0])
    reference = np.column_stack((x, radius, angle))
    return SimpleNamespace(
        cross_section=CrossSectionSquared((0.001) / 2) if section is None else section,
        reference_path=reference,
        x_domain=x,
        channel_heights=np.full_like(x, 0.012),
        n_stations=n_stations,
        channel_local_sector=np.full_like(x, 0.1 * np.cos(gamma)),
        helix_angle_vals=np.full_like(x, gamma),
        channel_wall_thickness=lambda x: np.full_like(x, 0.001),
        hot_side_layer_thickness=lambda x: np.full_like(x, hot_layer),
        total_thickness=lambda x: np.full_like(x, hot_layer + 0.001),
    )


@pytest.mark.parametrize("slope", [-0.6, 0.0, 0.6])
@pytest.mark.parametrize("gamma", [-0.7, 0.0, 0.7])
def test_squared_arcs_remain_on_transverse_circles_with_projected_height(slope, gamma):
    circuit = circuit_geometry(slope=slope, gamma=gamma, hot_layer=0.002)
    construction = build_channel_construction(circuit)
    expected_scale = np.sqrt(1 + slope**2)
    np.testing.assert_allclose(construction.pitch, 0.1)
    np.testing.assert_allclose(construction.alpha, np.arctan(slope), atol=1e-14)
    np.testing.assert_allclose(construction.gamma, gamma, atol=1e-14)
    np.testing.assert_allclose(construction.planes.normals, [[1, 0, 0]] * 9)
    for i in range(9):
        polygon = construction.polygon(i, 31)
        np.testing.assert_allclose(polygon[:, 0], circuit.x_domain[i])
        inner = 0.1 + slope * circuit.x_domain[i] + 0.003 * expected_scale
        outer = inner + 0.012 * expected_scale
        np.testing.assert_allclose(np.linalg.norm(polygon[:31, 1:], axis=1), inner)
        np.testing.assert_allclose(np.linalg.norm(polygon[31:, 1:], axis=1), outer)


@pytest.mark.parametrize("slope", [-0.6, 0.0, 0.6])
@pytest.mark.parametrize("gamma", [-0.7, 0.0, 0.7])
def test_rib_faces_point_at_the_engine_axis(slope, gamma):
    circuit = circuit_geometry(slope=slope, gamma=gamma)
    construction = build_channel_construction(circuit)
    for i in range(9):
        polygon = construction.polygon(i, 8)
        # Floor-to-roof along each of the two rib faces, each against the rib
        # it belongs to. Neither leans, however the wall slopes and however
        # steeply the channel winds.
        for side, face in ((1, polygon[8] - polygon[7]), (-1, polygon[15] - polygon[0])):
            rib_angle = circuit.reference_path[i, 2] + side * 0.05
            radial = np.array([0, np.cos(rib_angle), np.sin(rib_angle)])
            np.testing.assert_allclose(np.cross(face, radial), 0.0, atol=1e-17)


@pytest.mark.parametrize("slope", [0.0, 0.6])
@pytest.mark.parametrize("gamma", [0.0, 0.7])
def test_neighbouring_channels_leave_the_requested_material_between_them(slope, gamma):
    """The rib is what is left in the solid, not the gap across one section.

    Both faces of a rib are helicoids, so their closest approach is reached at
    a different station, and at the channel roof, where the same angular rate
    is a steeper helix. That full width is twice ``sidewall_half_thickness``.
    """
    circuit = circuit_geometry(slope=slope, gamma=gamma, n_stations=801)
    first = build_channel_construction(circuit)
    neighbour_path = circuit.reference_path.copy()
    neighbour_path[:, 2] += 0.1
    second = build_channel_construction(circuit, neighbour_path)

    # Sample the neighbour's facing rib wall over every station and height.
    edges = np.array([second.polygon(j, 8)[[0, -1]] for j in range(len(neighbour_path))])
    fraction = np.linspace(0, 1, 201)[:, None]
    wall = (edges[:, None, 0]*(1 - fraction) + edges[:, None, 1]*fraction).reshape(-1, 3)

    floor, roof = first.polygon(len(neighbour_path)//2, 8)[[7, 8]]
    material = np.array([
        np.min(np.linalg.norm(wall - ((1 - f)*floor + f*roof), axis=1))
        for f in np.linspace(0, 1, 9)
    ])

    # Never thinner than asked for, and no more than a per cent over: on a
    # taper the neighbouring roof edge has moved outward by the time the two
    # faces close on each other, which only leaves extra material.
    assert material.min() >= 0.001 * (1.0 - 1e-9)
    assert material.min() == pytest.approx(0.001, rel=1e-2)
    if gamma:
        # A straight channel has no twist to thin the rib; a helical one does,
        # and it thins towards the roof.
        assert material.argmin() >= len(material) // 2
        assert material[0] > material.min() * 1.005
    else:
        np.testing.assert_allclose(material, 0.001, rtol=1e-6)


def transverse_point(construction, index, u, v):
    """Follow the flow direction from a normal-plane point to x = x[index].

    A locally prismatic channel meets the axis-normal plane in exactly this
    swept image, so the gas area at a station is what these images leave open.
    """
    planes = construction.flow_planes
    point = (planes.origins[index] + u * planes.x_axes[index]
             + v * planes.y_axes[index])
    normal = planes.normals[index]
    point = point + (construction.reference_path[index, 0] - point[0]) / normal[0] * normal
    return point[1:]


def ray_offsets(construction, index, side):
    """Return the affine map from normal-plane (u, v) to sector-ray offset."""
    angle = construction.reference_path[index, 2] + side * construction.pitch[index] / 2
    normal = np.array([-np.sin(angle), np.cos(angle)])
    origin = transverse_point(construction, index, 0.0, 0.0) @ normal
    return origin, np.array([
        transverse_point(construction, index, 1.0, 0.0) @ normal - origin,
        transverse_point(construction, index, 0.0, 1.0) @ normal - origin,
    ])


@pytest.mark.parametrize("section_type", [CrossSectionRounded])
def test_rounded_lobes_keep_the_gas_area_of_the_sector_they_replace(section_type):
    circuit = circuit_geometry(slope=0.3, gamma=0.5, section=section_type())
    construction = build_channel_construction(circuit)
    assert construction.section_type == "rounded"
    radius = circuit.reference_path[:, 1]
    np.testing.assert_allclose(
        construction.gas_area, construction.pitch * radius**2 / 2, rtol=1e-12,
    )
    # The lobe centre is where the section is anchored, and it sits further out
    # than the contour it replaces because the tube has to clear the sector.
    base = construction.profiles.cross_section_base_point_path
    np.testing.assert_allclose(base[:, 0], circuit.x_domain)
    assert np.all(base[:, 1] > radius)
    for index in range(len(radius)):
        polygon = construction.polygon(index, 17)
        offset = polygon - construction.planes.origins[index]
        np.testing.assert_allclose(offset @ construction.planes.normals[index], 0, atol=1e-16)


@pytest.mark.parametrize("slope", [0.0, 0.6])
@pytest.mark.parametrize("gamma", [0.0, 0.7])
def test_rounded_gas_side_walls_stop_exactly_on_the_sector_boundary(slope, gamma):
    # The lobes are sized at the radius they end up at, so both tubes of the
    # section touch the sector rays that separate them from their neighbours.
    circuit = circuit_geometry(slope=slope, gamma=gamma, section=CrossSectionRounded())
    circuit.channel_heights[:] = 0.011
    construction = build_channel_construction(circuit)
    prof = construction.profiles
    tube = prof.sector_radius * np.sin(prof.theta / 2)
    base = circuit.cross_section.base_height(prof)
    for index in range(len(prof.h)):
        for side in (-1.0, 1.0):
            origin, gradient = ray_offsets(construction, index, side)
            for centre, gas_radius in (
                (np.zeros(2), tube[index]),
                (np.array([base[index], 0.0]), tube[index] + base[index] * np.sin(prof.theta[index] / 2)),
            ):
                reach = origin + gradient @ centre + side * gas_radius * np.linalg.norm(gradient)
                assert reach == pytest.approx(0.0, abs=1e-15)


def test_minimum_rounded_height_forms_a_circle_without_side_segments():
    circuit = circuit_geometry(section=CrossSectionRounded())
    # Without slope or helix the sector is the plain transverse wedge, so the
    # tangent tube's centre distance follows from the gas-area constraint.
    half = 0.05
    factor = np.sin(half) * np.cos(half) - (np.pi / 2 - half) * np.sin(half)**2
    distance = np.sqrt(0.1**2 * half / factor)
    # Lobe sizing does not depend on the height, so the first pass reports the
    # minimum that the second pass is then built at.
    sized = build_channel_construction(circuit)
    np.testing.assert_allclose(sized.profiles.sector_radius, distance)
    circuit.channel_heights[:] = circuit.cross_section.minimum_channel_height(sized.profiles)
    np.testing.assert_allclose(circuit.channel_heights, 2 * (distance * np.sin(half) - 0.001))
    construction = build_channel_construction(circuit)
    assert np.all(construction.circular)
    np.testing.assert_allclose(construction.inner_arc.centers, construction.outer_arc.centers)
    np.testing.assert_allclose(construction.inner_arc.radii, construction.outer_arc.radii)
    np.testing.assert_allclose(
        construction.inner_arc.points(0, 2)[::-1], construction.outer_arc.points(0, 2), atol=1e-16,
    )


def test_construction_does_not_mutate_simulation_dimensions_or_rib_law_inputs():
    def rib(h, gamma, t_wall_total, t_hot_layer, pitch_mid):
        assert h == pytest.approx(0.012)
        assert gamma == pytest.approx(0.5)
        assert t_wall_total == pytest.approx(0.003)
        assert t_hot_layer == pytest.approx(0.002)
        return 0.05 * pitch_mid

    circuit = circuit_geometry(slope=0.6, gamma=0.5, hot_layer=0.002, section=CrossSectionSquared(rib))
    before = {key: value.copy() for key, value in vars(circuit).items() if isinstance(value, np.ndarray)}
    centerline = circuit.reference_path.copy()
    construction = build_channel_construction(circuit)
    for key, value in before.items():
        np.testing.assert_array_equal(getattr(circuit, key), value)
    np.testing.assert_array_equal(circuit.reference_path, centerline)
    np.testing.assert_allclose(construction.profiles.h, 0.012)
    np.testing.assert_allclose(construction.profiles.theta, 0.1 * np.cos(0.5))


def test_transverse_rib_overlap_is_rejected_even_if_analytic_pitch_is_open():
    circuit = circuit_geometry(section=CrossSectionSquared((0.010098) / 2))
    # 0.010098 < analytic floor pitch 0.0101, but exceeds the chord available.
    with pytest.raises(ValueError, match="no open channel"):
        build_channel_construction(circuit)


def test_numerical_layout_helix_change_is_used_for_rib_projection():
    circuit = circuit_geometry(gamma=0.2)
    path = circuit.reference_path.copy()
    path[:, 2] += 2 * path[:, 0]
    construction = build_channel_construction(circuit, path)
    np.testing.assert_allclose(construction.gamma, np.arctan(np.tan(0.2) + 0.2))
    np.testing.assert_allclose(construction.pitch, 0.1)


@pytest.mark.parametrize("gamma", [np.nan, np.pi / 2, -np.pi / 2])
def test_singular_or_nonfinite_helix_is_rejected(gamma):
    circuit = circuit_geometry()
    circuit.helix_angle_vals[:] = gamma
    with pytest.raises(ValueError, match="helix angle"):
        build_channel_construction(circuit)


def test_viewer_stride_retains_the_final_station_and_lofts_a_closed_mesh():
    from pyskyfire.viz.engine_viz import build_sections_for_centerlines, loft_sections_to_polydata

    circuit = circuit_geometry(slope=0.3, gamma=0.5)
    sections = build_sections_for_centerlines(
        [circuit.reference_path], [0], SimpleNamespace(cooling_circuits=[circuit]), stride=3,
    )[0]
    assert len(sections) == 4  # 0, 3, 6, 8, including the final partial stride.
    np.testing.assert_allclose(sections[-1][:, 0], circuit.x_domain[-1])
    mesh = loft_sections_to_polydata(sections, cap_ends=True)
    assert mesh.n_open_edges == 0
    assert mesh.volume > 0
