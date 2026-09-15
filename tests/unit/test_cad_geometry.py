"""Geometry-only transverse export and hot-wall assembly contracts."""

from types import SimpleNamespace

import numpy as np
import pytest

from pyskyfire.regen.cad_geometry import (
    PlanarArc,
    TangentBezier,
    build_hot_wall,
    hot_wall_regions,
    hot_wall_revolve,
    join_hot_wall,
    transverse_construction,
)
from pyskyfire.regen.channel_construction import build_channel_construction
from pyskyfire.regen.channel_layout import ChannelLayout, ResolvedChannel
from pyskyfire.regen.cross_section import CrossSectionRounded, CrossSectionSquared


def circuit(
    span=(0.0, 1.0), *, rounded=False, slope=0.0, helix=0.0, count=32, resolution=21
):
    x = np.linspace(*span, resolution)
    r = 0.1 + slope * (x - span[0])
    theta = helix * x
    return SimpleNamespace(
        name="circuit",
        reference_path=np.column_stack((x, r, theta)),
        x_domain=x,
        contour=SimpleNamespace(
            dr_dx=lambda x: slope, r=lambda x: 0.1 + slope * (x - span[0])
        ),
        channel_heights=np.full_like(x, 0.1),
        channel_local_sector=np.full_like(x, 2 * np.pi / count),
        channel_angular_width=np.full_like(x, 2 * np.pi / count),
        channel_wall_thickness=lambda x: np.full_like(np.asarray(x), 0.001),
        cross_section=CrossSectionRounded()
        if rounded
        else CrossSectionSquared(sidewall_half_thickness=(0.001) / 2),
    )


def layout_chamber(circuits, count=32, phases=None, contour_span=(0.0, 1.0)):
    channels = []
    for i, c in enumerate(circuits):
        path = c.reference_path.copy()
        if phases is not None:
            path[:, 2] += phases[i]
        channels.append(
            ResolvedChannel(f"channel_{i}", i, 0, path, count, 2 * np.pi / count)
        )
    contour = SimpleNamespace(xs=np.array(contour_span), r=circuits[0].contour.r)
    return SimpleNamespace(cooling_circuits=circuits, contour=contour), ChannelLayout(
        "test", tuple(channels)
    )


@pytest.mark.parametrize(
    "slope,helix", [(0.0, 0.0), (0.2, 0.0), (0.0, 2.0), (0.2, 2.0), (-0.03, -2.0)]
)
def test_transverse_sections_are_local_sweep_intersections_without_simulation_changes(
    slope, helix
):
    c = circuit(rounded=True, slope=slope, helix=helix)
    before = c.reference_path.copy()
    source = build_channel_construction(c)
    centroid = source.centroid_path(c.cross_section).copy()
    transverse = transverse_construction(source)
    for i in (0, 10, 20):
        original = source.polygon(i, 25)
        normal = source.flow_planes.normals[i]
        expected = (
            original
            - (original[:, :1] - source.planes.origins[i, 0]) * normal / normal[0]
        )
        np.testing.assert_allclose(transverse.polygon(i, 25), expected, atol=1e-15)
        np.testing.assert_allclose(transverse.polygon(i, 25)[:, 0], c.x_domain[i])
    np.testing.assert_array_equal(c.reference_path, before)
    np.testing.assert_array_equal(source.centroid_path(c.cross_section), centroid)
    kind = transverse.inner_arcs[0].to_data()["kind"]
    assert kind == ("circle" if slope == helix == 0 else "ellipse")


def test_principal_ellipse_parameters_preserve_arc_geometry_and_direction():
    arc = PlanarArc(
        np.array([0.1, 0.02]), np.array([0.01, 0.003]), np.array([0.0, 0.02]), 4.5, 1.7
    )
    data = arc.to_data()
    major = np.asarray(data["majorAxis"])
    minor = np.array([-major[1], major[0]])
    t0, t1 = data["startParameter"], data["endParameter"]
    angles = 2 * np.pi * np.linspace(t0, t0 + (t1 - t0) % 1, 31)
    points = (
        np.asarray(data["center_uv_mm"])
        + data["majorRadius_mm"] * np.cos(angles)[:, None] * major
        + data["minorRadius_mm"] * np.sin(angles)[:, None] * minor
    ) / 1000
    np.testing.assert_allclose(points, arc.point(np.linspace(1, 0, 31)), atol=1e-15)


@pytest.mark.parametrize("rounded", [False, True])
def test_complete_wall_profiles_and_tangent_connectors(rounded):
    c = circuit(rounded=rounded, helix=0.4, slope=0.03)
    chamber, layout = layout_chamber([c])
    wall = build_hot_wall(chamber, layout)
    if not rounded:
        assert wall["sections"] == wall["runs"] == []
        assert len(wall["revolves"]) == 1
        assert wall["fullContourAllowed"]
        return
    assert not wall["fullContourAllowed"]
    assert len(wall["sections"]) == 21
    for section in wall["sections"]:
        curves = section["curves"]
        assert sum(d["kind"] == "bezier" for d in curves) == (32 if rounded else 0)
        for i, d in enumerate(curves):
            if d["kind"] != "bezier":
                continue
            before, after = curves[i - 1], curves[(i + 1) % len(curves)]

            def from_data(d):
                return PlanarArc(
                    np.array(d["center_uv_mm"]),
                    np.array(d["axisA_uv_mm"]),
                    np.array(d["axisB_uv_mm"]),
                    np.radians(d["startAngle_deg"]),
                    np.radians(d["endAngle_deg"]),
                )

            p = np.array(d["controlPoints_uv_mm"])
            first, second = from_data(before), from_data(after)
            np.testing.assert_allclose(first.point(1), p[0], atol=1e-11)
            np.testing.assert_allclose(second.point(0), p[-1], atol=1e-11)
            np.testing.assert_allclose(
                (p[1] - p[0]) / np.linalg.norm(p[1] - p[0]),
                first.tangent(1),
                atol=1e-12,
            )
            np.testing.assert_allclose(
                (p[-1] - p[-2]) / np.linalg.norm(p[-1] - p[-2]),
                second.tangent(0),
                atol=1e-12,
            )


def test_mixed_wall_and_internal_gap_export():
    square, rounded = circuit(count=64), circuit(count=64, rounded=True)
    chamber, layout = layout_chamber([square, rounded], phases=[0, np.pi / 32])
    wall = build_hot_wall(chamber, layout, trim_fraction=0.15)
    assert len(wall["sections"][0]["curves"]) == 128
    a, b = circuit((0, 0.495)), circuit((0.505, 0.95), rounded=True)
    chamber, layout = layout_chamber([a, b])
    wall = build_hot_wall(chamber, layout, transition_fraction=0.04)
    positions = [s["plane"]["origin_mm"][2] / 1000 for s in wall["sections"]]
    assert positions[0] == 0.505 and positions[-1] == 0.95
    assert wall["revolves"][0]["span_mm"] == [0, 495]
    assert wall["excludedInterleaveIntervals_mm"] == []
    assert "bridges" not in wall


def test_bezier_trim_control_changes_shape():
    pieces = []
    for i in range(8):
        theta = i * 2 * np.pi / 8
        arc = PlanarArc(
            np.zeros(2),
            np.array([0.1, 0]),
            np.array([0, 0.1]),
            theta - 0.3,
            theta + 0.3,
        )
        pieces.append((str(i), "rounded", arc, theta))
    first, second = join_hot_wall(pieces, 0.1), join_hot_wall(pieces, 0.2)
    assert isinstance(first[1][1], TangentBezier)
    assert not np.allclose(first[1][1].control_points, second[1][1].control_points)


def test_helical_profile_seam_stays_on_the_same_channel():
    c = circuit(rounded=True, helix=6.0)
    chamber, layout = layout_chamber([c])
    wall = build_hot_wall(chamber, layout)
    names = [[curve["name"] for curve in s["curves"]] for s in wall["sections"]]
    assert all(row == names[0] for row in names)


def test_one_helical_circuit_is_a_single_run_cad_can_pair_curve_by_curve():
    """CAD pairs curve j across a run, so the run must name what it pairs."""
    c = circuit(rounded=True, helix=6.0)
    chamber, layout = layout_chamber([c])
    wall = build_hot_wall(chamber, layout)

    (run,) = wall["runs"]
    assert run["sectionIndices"] == list(range(len(wall["sections"])))
    for index in run["sectionIndices"]:
        section = wall["sections"][index]
        assert [curve["name"] for curve in section["curves"]] == run["curveNames"]


def _curve_from_data(data):
    if data["kind"] == "bezier":
        return TangentBezier(np.array(data["controlPoints_uv_mm"]))
    return PlanarArc(
        np.array(data["center_uv_mm"]),
        np.array(data["axisA_uv_mm"]),
        np.array(data["axisB_uv_mm"]),
        np.radians(data["startAngle_deg"]),
        np.radians(data["endAngle_deg"]),
    )


def test_runs_reject_a_station_that_is_not_a_construction_sample():
    """A drifting station would silently take a neighbour's arc."""
    from pyskyfire.regen.cad_geometry import _station_index

    stations = np.array([0.0, 0.5, 1.0])
    assert _station_index(stations, 0.5) == 1
    assert _station_index(stations, 1.0) == 2
    with pytest.raises(ValueError, match="not a construction sample"):
        _station_index(stations, 0.75)


@pytest.mark.parametrize("helix", [0, 6])
def test_squared_wall_revolves_only_covered_union_and_offers_full_contour(helix):
    a, b = circuit((0.1, 0.4), helix=helix), circuit((0.6, 0.9), helix=helix)
    chamber, layout = layout_chamber([a, b])
    with pytest.warns(UserWarning, match="not covered"):
        wall = build_hot_wall(chamber, layout)
    assert [p["span_mm"] for p in wall["revolves"]] == [[100, 400], [600, 900]]
    assert wall["fullContour"]["span_mm"] == [0, 1000]
    assert wall["sections"] == wall["runs"] == []
    assert wall["uncoveredLengthFraction"] == pytest.approx(0.4)


def test_squared_interleave_uses_one_unbroken_revolve():
    chamber, layout = layout_chamber([circuit((0, 0.6)), circuit((0.4, 1))])
    wall = build_hot_wall(chamber, layout, transition_fraction=0.2)
    assert [p["span_mm"] for p in wall["revolves"]] == [[0, 1000]]
    assert wall["excludedInterleaveIntervals_mm"] == []


def test_rounded_interleave_leaves_percentage_gap_and_preserves_shared_span():
    a, b = circuit((0, 0.6), count=64), circuit((0.4, 1), count=64, rounded=True)
    chamber, layout = layout_chamber([a, b], phases=[0, np.pi / 32])
    wall = build_hot_wall(chamber, layout, transition_fraction=0.04)
    np.testing.assert_allclose(
        wall["excludedInterleaveIntervals_mm"], [[380, 420], [580, 620]]
    )
    assert wall["revolves"][0]["span_mm"] == [0, 380]
    positions = np.array([s["plane"]["origin_mm"][2] for s in wall["sections"]])
    runs = wall["runs"]
    assert len(runs) == 2
    np.testing.assert_allclose(
        positions[runs[0]["sectionIndices"]][[0, -1]], [420, 580]
    )
    np.testing.assert_allclose(
        positions[runs[1]["sectionIndices"]][[0, -1]], [620, 1000]
    )
    assert len(runs[0]["curveNames"]) == 128
    assert "fullContour" not in wall and not wall["fullContourAllowed"]


def test_repeated_rounded_run_does_not_loft_across_a_removed_margin():
    # A short secondary circuit is swallowed by overlapping margins. The same
    # primary curves remain on both sides but must still be separate lofts.
    a, b = circuit(rounded=True), circuit((0.49, 0.51))
    chamber, layout = layout_chamber([a, b])
    wall = build_hot_wall(chamber, layout, transition_fraction=0.04)
    first, second = wall["runs"]
    assert first["curveNames"] == second["curveNames"]
    np.testing.assert_allclose(wall["excludedInterleaveIntervals_mm"], [[470, 530]])
    assert wall["sections"][first["sectionIndices"][-1]]["plane"]["origin_mm"][2] == 470
    assert wall["sections"][second["sectionIndices"][0]]["plane"]["origin_mm"][2] == 530


@pytest.mark.parametrize("fraction", [0, -1, 1.01, np.nan, np.inf])
def test_invalid_transition_fraction_fails(fraction):
    with pytest.raises(ValueError, match="fraction"):
        hot_wall_regions([(0, 1)], [True], (0, 1), transition_fraction=fraction)


def test_transition_fraction_uses_full_contour_length():
    regions, exclusions, _ = hot_wall_regions(
        [(0, 1.2), (0.8, 2)], [False, False], (0, 2), transition_fraction=0.04
    )
    np.testing.assert_allclose(exclusions, [[0.76, 0.84], [1.16, 1.24]])
    np.testing.assert_allclose(
        [(r["start"], r["end"]) for r in regions],
        [(0, 0.76), (0.84, 1.16), (1.24, 2)],
    )


def test_coverage_warning_uses_union_and_strict_ten_percent_threshold():
    import warnings

    with warnings.catch_warnings(record=True) as caught:
        _, _, missing = hot_wall_regions([(0, 0.5), (0.3, 0.9)], [True, True], (0, 1))
    assert missing == pytest.approx(0.1)
    assert not caught
    with pytest.warns(UserWarning, match="not covered"):
        hot_wall_regions([(0, 0.899)], [True], (0, 1))


def _revolve_curve_points(curve):
    """The (radius, axial) mm points a revolve curve has to pass through."""
    if curve["kind"] == "spline":
        return np.array(curve["points_uv_mm"])
    return np.array(curve["controlPoints_uv_mm"])[[0, -1]]


def test_revolve_fits_one_spline_per_smooth_run_and_cuts_at_the_throat():
    from pyskyfire.regen import Contour

    contour = Contour(
        xs=[-1, -0.7, -0.3, 0, 0.4, 1],
        rs=[0.3, 0.3, 0.2, 0.1, 0.2, 0.4],
        chamber_angle=0.2,
        nozzle_angle=0.3,
    )
    for span, kinds in [
        ((-1, 1), ["spline", "spline"]),
        ((-0.81, 0.23), ["spline", "bezier"]),
    ]:
        profile = hot_wall_revolve(contour, *span)
        assert profile["span_mm"] == [v * 1000 for v in span]
        # Smooth runs become one curve each; the throat's independent slopes
        # stay a corner because they end up on a curve boundary.
        assert [curve["kind"] for curve in profile["curves"]] == kinds
        pieces = [_revolve_curve_points(curve) for curve in profile["curves"]]
        for points in pieces:
            np.testing.assert_allclose(
                points[:, 0] / 1000, contour.r(points[:, 1] / 1000), atol=1e-15
            )
        np.testing.assert_allclose(pieces[0][-1], pieces[1][0], atol=1e-15)
        assert pieces[0][-1][1] == 0.0
        assert [pieces[0][0][1], pieces[-1][-1][1]] == [v * 1000 for v in span]

    # A converging run interpolates every contour knot it covers, trimmed
    # spans included, so CAD never has to infer the shape between them.
    converging = hot_wall_revolve(contour, -0.81, 0.23)["curves"][0]
    np.testing.assert_allclose(
        np.array(converging["points_uv_mm"])[:, 1], [-810, -700, -300, 0]
    )


def test_revolve_cuts_a_spline_at_every_corner_of_a_faceted_contour():
    """An interpolating contour corners at each point; splines must not round it."""
    xs = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
    rs = np.array([0.3, 0.15, 0.1, 0.25, 0.3])
    contour = SimpleNamespace(xs=xs, r=lambda x: np.interp(x, xs, rs))

    curves = hot_wall_revolve(contour, -1.0, 1.0)["curves"]
    assert [curve["kind"] for curve in curves] == ["bezier"] * 4
    for curve, (a, b) in zip(curves, zip(rs, rs[1:])):
        control = np.array(curve["controlPoints_uv_mm"])
        np.testing.assert_allclose(control[[0, -1], 0] / 1000, [a, b], atol=1e-15)


def test_multiple_seeds_of_one_circuit_are_not_interleave_boundaries():
    c = circuit(rounded=True, count=64)
    chamber, layout = layout_chamber([c])
    path = c.reference_path.copy()
    path[:, 2] += np.pi / 32
    second = ResolvedChannel("second_seed", 0, 1, path, 32, np.pi / 16)
    layout = ChannelLayout("test", (*layout.channels, second))
    wall = build_hot_wall(chamber, layout)
    assert wall["excludedInterleaveIntervals_mm"] == []
    assert len(wall["runs"]) == 1


def test_mixed_helical_wall_preserves_square_and_rounded_curve_correspondence():
    circuits = [circuit(count=64, helix=6), circuit(count=64, rounded=True, helix=6)]
    chamber, layout = layout_chamber(circuits, phases=[0, np.pi / 32])
    wall = build_hot_wall(chamber, layout)
    assert wall["revolves"] == []
    (run,) = wall["runs"]
    first = wall["sections"][0]["curves"]
    last = wall["sections"][-1]["curves"]
    rotation = np.array([[np.cos(6), -np.sin(6)], [np.sin(6), np.cos(6)]])
    for start, end in zip(first, last):
        assert start["name"] == end["name"]
        points = _curve_from_data(start).point(np.linspace(0, 1, 5))
        np.testing.assert_allclose(
            _curve_from_data(end).point(np.linspace(0, 1, 5)),
            points @ rotation.T,
            atol=1e-10,
        )
    assert len(run["curveNames"]) == 128


def test_rounded_gap_stops_lofts_at_circuit_ends_without_transition():
    chamber, layout = layout_chamber(
        [circuit((0, 0.495), rounded=True), circuit((0.505, 1), rounded=True)]
    )
    wall = build_hot_wall(chamber, layout, transition_fraction=0.04)
    positions = np.array([s["plane"]["origin_mm"][2] for s in wall["sections"]])
    assert len(wall["runs"]) == 2
    np.testing.assert_allclose(
        positions[wall["runs"][0]["sectionIndices"]][[0, -1]], [0, 495]
    )
    np.testing.assert_allclose(
        positions[wall["runs"][1]["sectionIndices"]][[0, -1]], [505, 1000]
    )
    assert wall["excludedInterleaveIntervals_mm"] == []


def test_percentage_gap_can_consume_all_hot_wall_without_generating_a_bridge():
    chamber, layout = layout_chamber([circuit(), circuit((0.4, 0.6), rounded=True)])
    wall = build_hot_wall(chamber, layout, transition_fraction=1)
    assert wall["sections"] == wall["runs"] == wall["revolves"] == []
    assert wall["excludedInterleaveIntervals_mm"] == [[0, 1000]]
