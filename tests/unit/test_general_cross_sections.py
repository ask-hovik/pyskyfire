"""Geometric and thermal contracts shared by chamber and exchanger sections."""

from dataclasses import replace
import pickle

import numpy as np
import pytest

from pyskyfire.regen import CrossSectionRounded, CrossSectionSquared, SectionProfiles
from pyskyfire.common.section_geometry import SectionArc, SectionLine


def rectangle(**kw):
    return SectionProfiles.from_dimensions(
        floor_width=0.006, roof_width=0.006, height=0.01,
        wall_thickness=0.001, **kw,
    )


def test_rectangle_geometry_and_independent_wall_thicknesses():
    prof = rectangle()
    section = CrossSectionSquared(0.002, roof_thickness=0.003)
    assert all(isinstance(curve, SectionLine) for curve in section.boundary(prof))
    np.testing.assert_allclose(section.widths(prof), [[0.006], [0.006]])
    np.testing.assert_allclose(section.A_coolant(prof), 0.006*0.01)
    np.testing.assert_allclose(section.P_wetted(prof), 2*(0.006+0.01))
    np.testing.assert_allclose(section.Dh_coolant(prof), 4*0.006*0.01/(2*(0.006+0.01)))
    np.testing.assert_allclose(section.A_exterior(prof), (0.006+2*0.002)*(0.01+0.001+0.003))
    for name, value in {"floor": 0.001, "roof": 0.003, "side_left": 0.002, "side_right": 0.002}.items():
        np.testing.assert_allclose(section.surface_thicknesses(prof)[name], value)


@pytest.mark.parametrize("height", [0.006, 0.02])
def test_parallel_rounded_section_is_a_capsule_or_circle(height):
    width, wall = 0.006, 0.0005
    prof = SectionProfiles.from_dimensions(floor_width=width, roof_width=width,
                                          height=height, wall_thickness=wall)
    section = CrossSectionRounded()
    distance = height-width
    expected_area = np.pi*(width/2)**2 + width*distance
    expected_perimeter = np.pi*width + 2*distance
    np.testing.assert_allclose(section.A_coolant(prof), expected_area)
    np.testing.assert_allclose(section.P_wetted(prof), expected_perimeter)
    np.testing.assert_allclose(section.full_height(prof), height)
    np.testing.assert_allclose(section.base_height(prof), distance)
    np.testing.assert_allclose(section.widths(prof), [[width], [width]])
    np.testing.assert_allclose(section.A_exterior(prof),
                               np.pi*(width/2+wall)**2 + (width+2*wall)*distance)


@pytest.mark.parametrize("angle", [0.06, 0.2, 0.6])
def test_rounded_sector_and_endpoint_width_inputs_describe_the_same_curves(angle):
    section = CrossSectionRounded(height_input="base")
    prof = SectionProfiles.from_dimensions(radius=0.1, angle=angle, height=0.04,
                                          wall_thickness=0.0005)
    lower, upper = section.widths(prof)
    explicit = SectionProfiles.from_dimensions(floor_width=lower, roof_width=upper,
        height=section.full_height(prof), wall_thickness=prof.t_wall)
    direct = CrossSectionRounded()
    for actual, expected in zip(direct.boundary(explicit), section.boundary(prof)):
        np.testing.assert_allclose(actual.points(0, 31), expected.points(0, 31), atol=1e-13)
    np.testing.assert_allclose(direct.A_coolant(explicit), section.A_coolant(prof))


@pytest.mark.parametrize("lower,upper", [(0.004, 0.007), (0.007, 0.004)])
def test_unequal_rounded_chords_are_exact_and_sidewalls_are_tangent(lower, upper):
    prof = SectionProfiles.from_dimensions(floor_width=lower, roof_width=upper,
                                          height=0.02, wall_thickness=0.001)
    section = CrossSectionRounded()
    np.testing.assert_allclose(section.widths(prof), [[lower], [upper]])
    np.testing.assert_allclose(section.full_height(prof), 0.02)
    surfaces = section.surfaces(prof)
    ordered = list(surfaces.values())
    for first, second in zip(ordered, ordered[1:]+ordered[:1]):
        np.testing.assert_allclose(first.evaluate([1]), second.evaluate([0]), atol=1e-16)
        a, b = first.derivative([1])[0, 0], second.derivative([0])[0, 0]
        np.testing.assert_allclose(np.cross(np.r_[a, 0], np.r_[b, 0]), 0, atol=1e-15)
        assert a@b > 0
    for name, inside in surfaces.items():
        outside = section.surfaces(prof, exterior=True)[name]
        distance = np.linalg.norm(outside.evaluate([0.5])-inside.evaluate([0.5]), axis=-1)
        np.testing.assert_allclose(distance, 0.001)


@pytest.mark.parametrize("floor", ["arc", "line"])
@pytest.mark.parametrize("roof", ["arc", "line"])
def test_squared_chords_do_not_depend_on_the_end_curve_kind(floor, roof):
    prof = rectangle(radius=0.05)
    section = CrossSectionSquared(0.001, floor=floor, roof=roof)
    np.testing.assert_allclose(section.widths(prof), [[0.006], [0.006]])
    for curve, kind in zip(section.boundary(prof), (floor, roof)):
        assert isinstance(curve, SectionArc if kind == "arc" else SectionLine)
    if floor == roof == "line":
        np.testing.assert_allclose(section.A_coolant(prof), 0.00006)
    assert section.A_exterior(prof)[0] > section.A_coolant(prof)[0]


def test_sector_line_selection_preserves_sidewall_endpoints_and_rib_thickness():
    prof = SectionProfiles.from_dimensions(radius=0.1, angle=0.15, height=0.01,
                                          wall_thickness=0.001)
    curved = CrossSectionSquared(0.0004)
    straight = CrossSectionSquared(0.0004, floor="line", roof="line")
    np.testing.assert_allclose(curved.widths(prof), straight.widths(prof))
    for a, b in zip(curved.boundary(prof), straight.boundary(prof)):
        np.testing.assert_allclose(a.points(0, 2), b.points(0, 2))
    right = straight.surfaces(prof)["side_right"]
    normal = np.array([-np.sin(0.15/2), np.cos(0.15/2)])
    points = right.points(0, 2) + [0.1, 0]
    np.testing.assert_allclose(points@normal, -0.0004)
    floor, roof = straight.boundary(prof)
    np.testing.assert_allclose(straight.full_height(prof),
        roof.evaluate([0.5])[:, 0, 0]-floor.evaluate([0.5])[:, 0, 0])


def test_tapered_squared_sidewall_thickness_is_normal_to_its_face():
    prof = SectionProfiles.from_dimensions(floor_width=0.004, roof_width=0.009,
        height=0.01, wall_thickness=0.001)
    section = CrossSectionSquared(0.0007, roof_thickness=0.002)
    inside, outside = section.surfaces(prof), section.surfaces(prof, exterior=True)
    for name in ("side_left", "side_right"):
        side = inside[name]
        tangent = side.ends-side.starts
        normal = np.column_stack((-tangent[:, 1], tangent[:, 0]))/np.linalg.norm(tangent, axis=1)[:, None]
        np.testing.assert_allclose(np.sum((outside[name].starts-side.starts)*normal, axis=1), 0.0007)


def test_exposure_only_changes_thermal_lengths():
    prof = rectangle()
    section = CrossSectionSquared(0.001, exposure={
        "floor": 0.5, "roof": 0.25, "side_left": 0.2, "side_right": 0.0})
    plain = CrossSectionSquared(0.001)
    np.testing.assert_array_equal(section.A_coolant(prof), plain.A_coolant(prof))
    np.testing.assert_array_equal(section.Dh_coolant(prof), plain.Dh_coolant(prof))
    np.testing.assert_allclose(section.P_coolant(prof), 0.5*0.006+0.25*0.006+0.2*0.01)
    np.testing.assert_allclose(section.R_coolant_per_len(prof, 1000, 20),
        1/(1000*(0.5*0.006+0.25*0.006+0.2*0.01)))
    inactive = CrossSectionSquared(0.001, exposure={})
    assert inactive.P_coolant(prof)[0] == 0
    assert np.isinf(inactive.R_coolant_per_len(prof, 1000, 20)[0])


def test_angular_exposure_reproduces_the_former_112_degree_contact():
    angle = 112*np.pi/180
    prof = SectionProfiles.from_dimensions(radius=0.1, angle=0.2, height=0.025,
                                          wall_thickness=0.001)
    section = CrossSectionRounded(exposure_angles={"floor": angle})
    q = 0.1*np.sin(0.1)-0.001
    np.testing.assert_allclose(section.P_coolant(prof), q*angle)
    np.testing.assert_allclose(section.P_thermal(prof), (q+0.001)*angle)
    fractional = CrossSectionRounded(exposure={"floor": angle/(np.pi-0.2)})
    np.testing.assert_allclose(fractional.P_coolant(prof), section.P_coolant(prof))
    np.testing.assert_allclose(fractional.P_thermal(prof), section.P_thermal(prof))


@pytest.mark.parametrize("side", ["side_left", "side_right"])
def test_any_direct_sidewall_exposure_conflicts_with_symmetric_fins(side):
    with pytest.raises(ValueError, match="mutually exclusive"):
        CrossSectionSquared(0.001, exposure={side: 0.001}, enable_fin=True)
    section = CrossSectionRounded(exposure={side: lambda x: x}, enable_fin=True)
    prof = SectionProfiles.from_dimensions(floor_width=0.006, roof_width=0.006,
        height=0.01, wall_thickness=0.001, stations=[0, 0.5])
    with pytest.raises(ValueError, match="mutually exclusive"):
        section.R_coolant_per_len(prof, 1000, 20)


def test_symmetric_fins_work_on_a_standalone_section():
    prof = rectangle()
    section = CrossSectionSquared(0.001, enable_fin=True)
    bare = section.R_coolant_per_len(prof, 1000, 20, enable_fin=False)
    enhanced = section.R_coolant_per_len(prof, 1000, 20)
    assert enhanced[0] < bare[0]
    m = np.sqrt(2*1000/(20*0.002))
    efficiency = np.tanh(m*0.01)/(m*0.01)
    np.testing.assert_allclose(enhanced, 1/(1000*(0.006+efficiency*0.02)))


def test_local_station_callables_and_pickle_roundtrip():
    prof = rectangle(stations=[0, 0.5, 1])
    section = CrossSectionSquared(lambda x: 0.001+0.001*x,
        roof_thickness=lambda t_wall: 2*t_wall, exposure={"floor": lambda x: x})
    np.testing.assert_allclose(section.sidewall_half_thickness(prof), [0.001, 0.0015, 0.002])
    np.testing.assert_allclose(section.roof_thickness(prof), 0.002)
    np.testing.assert_allclose(section.P_coolant(prof), [0, 0.003, 0.006])
    simple = CrossSectionSquared(0.001, roof_thickness=0.003, exposure={"roof": 0.5})
    restored = pickle.loads(pickle.dumps(simple))
    np.testing.assert_allclose(restored.P_thermal(prof), simple.P_thermal(prof))


@pytest.mark.parametrize("fraction", [-0.1, 1.1, np.nan, np.inf])
def test_invalid_exposure_fractions_raise(fraction):
    with pytest.raises(ValueError, match="fraction"):
        CrossSectionRounded(exposure={"floor": fraction})


def test_invalid_angles_and_input_combinations_raise():
    with pytest.raises(ValueError, match="both exposure"):
        CrossSectionRounded(exposure={"floor": 0.5}, exposure_angles={"floor": 1})
    with pytest.raises(ValueError, match="only supported"):
        CrossSectionRounded(exposure_angles={"side_left": 1})
    with pytest.raises(ValueError, match="exceeds"):
        CrossSectionRounded(exposure_angles={"floor": 4}).P_coolant(rectangle())
    with pytest.raises(ValueError, match="both"):
        SectionProfiles.from_dimensions(floor_width=0.004, height=0.01)
    with pytest.raises(ValueError, match="without angle"):
        SectionProfiles.from_dimensions(floor_width=0.004, roof_width=0.004, height=0.01, angle=0.1)
    with pytest.raises(ValueError, match="use widths"):
        SectionProfiles.from_dimensions(radius=np.inf, angle=0, height=0.01)
    invalid = replace(rectangle(), h=np.array([0.001]))
    with pytest.raises(ValueError, match="geometric minimum"):
        CrossSectionRounded().A_coolant(invalid)
    with pytest.raises(ValueError, match="diameter"):
        CrossSectionSquared(0.001).A_coolant(replace(rectangle(), section_radius=np.array([0.001])))


def test_enabling_another_surface_does_not_change_the_floor_exposure():
    prof = SectionProfiles.from_dimensions(radius=0.1, angle=0.1, height=0.01,
                                          wall_thickness=0.001)
    floor = CrossSectionSquared(0.001, exposure={"floor": 0.5})
    both = CrossSectionSquared(0.001, exposure={"floor": 0.5, "roof": 0.5})
    np.testing.assert_array_equal(floor.thermal_perimeters(prof, exterior=True)["floor"],
                                 both.thermal_perimeters(prof, exterior=True)["floor"])


def test_minimum_height_for_unequal_chords_is_a_real_geometric_bound():
    section = CrossSectionRounded()
    prof = SectionProfiles.from_dimensions(floor_width=0.004, roof_width=0.007,
                                          height=0.03, wall_thickness=0.001)
    minimum = section.minimum_channel_height(prof)
    assert minimum[0] > max(prof.floor_width[0], prof.roof_width[0])
    limiting = replace(prof, h=minimum)
    np.testing.assert_allclose(section.full_height(limiting), minimum)
    with pytest.raises(ValueError, match="minimum"):
        section.A_coolant(replace(prof, h=minimum*0.999))


def test_straight_end_curves_reach_the_chamber_construction_and_export():
    from types import SimpleNamespace
    from pyskyfire.regen import build_channel_construction
    from pyskyfire.regen.cad_geometry import transverse_construction

    x = np.linspace(0, 0.1, 4)
    circuit = SimpleNamespace(
        cross_section=CrossSectionSquared(0.0005, floor="line", roof="line"),
        reference_path=np.column_stack((x, np.full_like(x, 0.1), np.zeros_like(x))),
        x_domain=x, channel_heights=np.full_like(x, 0.01),
        channel_local_sector=np.full_like(x, 0.1),
        channel_wall_thickness=lambda x: np.full_like(x, 0.001),
    )
    construction = build_channel_construction(circuit)
    export = transverse_construction(construction)
    for i in range(len(x)):
        for curve in (export.inner_arcs[i], export.outer_arcs[i]):
            data = curve.to_data()
            assert data["kind"] == "line"
            np.testing.assert_allclose(data["start_uv_mm"], curve.point(0)*1000)
            np.testing.assert_allclose(data["end_uv_mm"], curve.point(1)*1000)
        np.testing.assert_allclose(construction.polygon(i, 11), export.polygon(i, 11))


def test_removed_section_class_and_argument_are_not_exported():
    import pyskyfire.regen as regen
    assert not hasattr(regen, "CrossSectionRoundedInternal")
    with pytest.raises(TypeError):
        CrossSectionSquared(rib_width=0.001)
