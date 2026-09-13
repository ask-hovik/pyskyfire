"""Geometry-only transverse export and hot-wall assembly contracts."""

from types import SimpleNamespace

import numpy as np
import pytest

from pyskyfire.regen.cad_geometry import (
    PlanarArc,
    TangentBezier,
    build_hot_wall,
    hot_wall_stations,
    hot_wall_bridge,
    join_hot_wall,
    topology_runs,
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
        else CrossSectionSquared(rib_width=0.001),
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
    contour = SimpleNamespace(xs=np.array(contour_span))
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


def test_transition_stations_do_not_shrink_with_resolution():
    spans = [(0.0, 0.6), (0.4, 1.0)]
    for n in (11, 10001):
        stations, intervals, uncovered = hot_wall_stations(
            spans, np.linspace(0, 1, n), (0, 1), transition_fraction=0.1
        )
        np.testing.assert_allclose(intervals, [[0.35, 0.45], [0.55, 0.65]])
        for a, b in intervals:
            assert a in stations and b in stations
            assert not np.any((stations > a) & (stations < b))
        assert uncovered == 0


def test_gap_is_one_transition_and_extent_stops_at_circuits():
    with pytest.warns(UserWarning, match="not covered"):
        stations, intervals, uncovered = hot_wall_stations(
            [(0.1, 0.45), (0.46, 0.8)],
            np.linspace(0, 1, 101),
            (0, 1),
            transition_length_m=0.1,
        )
    np.testing.assert_allclose(intervals, [[0.405, 0.505]])
    assert stations[0] == 0.1 and stations[-1] == 0.8
    assert uncovered == pytest.approx(0.31)


def test_long_gap_retains_both_circuit_endpoints():
    with pytest.warns(UserWarning):
        stations, intervals, _ = hot_wall_stations(
            [(0, 0.3), (0.6, 1)], np.linspace(0, 1, 11), (0, 1)
        )
    assert intervals == [(0.3, 0.6)]
    assert 0.3 in stations and 0.6 in stations
    assert not np.any((stations > 0.3) & (stations < 0.6))


def test_coverage_is_union_and_warning_is_strictly_more_than_ten_percent():
    import warnings

    with warnings.catch_warnings(record=True) as caught:
        _, _, missing = hot_wall_stations(
            [(0, 0.5), (0.3, 0.9)], np.linspace(0, 0.9, 30), (0, 1)
        )
    assert missing == pytest.approx(0.1)
    assert not caught
    with pytest.warns(UserWarning):
        hot_wall_stations([(0, 0.899)], [0, 0.899], (0, 1))


def test_overlapping_transitions_and_invalid_controls_fail():
    with pytest.raises(ValueError, match="Overlapping"):
        hot_wall_stations(
            [(0, 0.51), (0.5, 1)], [0, 0.5, 1], (0, 1), transition_fraction=0.1
        )
    with pytest.raises(ValueError, match="beyond"):
        hot_wall_stations(
            [(0, 0.01), (0.01, 1)], [0, 0.01, 1], (0, 1), transition_fraction=0.1
        )
    for value in (0, -1, np.nan):
        with pytest.raises(ValueError):
            hot_wall_stations([(0, 1)], [0, 1], (0, 1), transition_length_m=value)


@pytest.mark.parametrize("rounded", [False, True])
def test_complete_wall_profiles_and_tangent_connectors(rounded):
    c = circuit(rounded=rounded, helix=0.4, slope=0.03)
    chamber, layout = layout_chamber([c])
    wall = build_hot_wall(chamber, layout)
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
    wall = build_hot_wall(chamber, layout, transition_length_m=0.04)
    positions = [s["plane"]["origin_mm"][2] / 1000 for s in wall["sections"]]
    assert positions[0] == 0 and positions[-1] == 0.95
    assert 0.48 in positions and 0.52 in positions
    assert not any(0.48 < x < 0.52 for x in positions)


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


def test_hand_offs_split_runs_contiguously_where_the_channels_change():
    a, b = circuit((0, 0.495)), circuit((0.505, 0.95), rounded=True)
    chamber, layout = layout_chamber([a, b])
    wall = build_hot_wall(chamber, layout, transition_length_m=0.04)

    runs = wall["runs"]
    assert len(runs) > 1
    covered = [index for run in runs for index in run["sectionIndices"]]
    assert covered == list(range(len(wall["sections"])))
    for run in runs:
        assert run["sectionIndices"] == list(
            range(run["sectionIndices"][0], run["sectionIndices"][-1] + 1)
        )
    for before, after in zip(runs, runs[1:]):
        assert before["curveNames"] != after["curveNames"]


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


def _hand_off_wall():
    a, b = circuit((0, 0.495)), circuit((0.505, 0.95), rounded=True)
    chamber, layout = layout_chamber([a, b])
    return build_hot_wall(chamber, layout, transition_length_m=0.04)


def test_a_hand_off_splits_both_profiles_onto_one_set_of_azimuths():
    """Unlike profiles cannot pair by index until they are cut to match."""
    wall = _hand_off_wall()

    (bridge,) = wall["bridges"]
    before, after = bridge["sections"]
    assert len(before["curves"]) == len(after["curves"]) == len(bridge["curveNames"])
    # Splitting only adds curve starts, so each side gains the other's.
    assert len(bridge["curveNames"]) > max(
        len(wall["sections"][bridge["beforeSection"]]["curves"]),
        len(wall["sections"][bridge["afterSection"]]["curves"]),
    )


def test_split_hand_off_profiles_stay_closed_and_wind_once():
    for section in _hand_off_wall()["bridges"][0]["sections"]:
        curves = [_curve_from_data(curve) for curve in section["curves"]]
        for first, second in zip(curves, curves[1:] + curves[:1]):
            np.testing.assert_allclose(first.point(1), second.point(0), atol=1e-10)
        points = np.vstack(
            [curve.point(np.linspace(0, 1, 41))[:-1] for curve in curves]
        )
        points = np.vstack([points, points[0]])
        angles = np.unwrap(np.arctan2(points[:, 1], points[:, 0]))
        assert np.all(np.diff(angles) > -1e-9)
        assert angles[-1] - angles[0] == pytest.approx(2 * np.pi)


def test_paired_hand_off_spans_start_at_the_same_azimuth():
    """Pairing by azimuth is what lets one loft carry several spans."""
    before, after = _hand_off_wall()["bridges"][0]["sections"]

    for first, second in zip(before["curves"], after["curves"]):
        angles = [
            np.arctan2(*_curve_from_data(curve).point(0.0)[::-1])
            for curve in (first, second)
        ]
        assert abs((angles[0] - angles[1] + np.pi) % (2 * np.pi) - np.pi) < 1e-6


def test_splitting_a_hand_off_profile_reproduces_its_curves_exactly():
    def ring(count, phase):
        return [
            (
                f"{count}_{i}",
                PlanarArc(
                    np.zeros(2),
                    np.array([1.0, 0.0]),
                    np.array([0.0, 1.0]),
                    (i + phase) * 2 * np.pi / count,
                    (i + 1 + phase) * 2 * np.pi / count,
                ),
            )
            for i in range(count)
        ]

    # Interleaved starts, half a curve apart, so none of them are merged.
    names, first, second = hot_wall_bridge(ring(5, 0.0), ring(5, 0.5))

    assert len(names) == len(first) == len(second) == 10
    for pieces in (first, second):
        points = np.vstack([curve.point(np.linspace(0, 1, 9)) for _, curve in pieces])
        # Both profiles are unit circles, so every split point stays on one.
        np.testing.assert_allclose(np.linalg.norm(points, axis=1), 1.0, atol=1e-12)


def test_a_borrowed_cut_too_close_to_a_curve_start_is_merged_into_it():
    """Cutting there would leave one profile a sliver the other never has."""

    def ring(count, phase):
        return [
            (
                f"{count}_{i}",
                PlanarArc(
                    np.zeros(2),
                    np.array([1.0, 0.0]),
                    np.array([0.0, 1.0]),
                    (i + phase) * 2 * np.pi / count,
                    (i + 1 + phase) * 2 * np.pi / count,
                ),
            )
            for i in range(count)
        ]

    names, first, second = hot_wall_bridge(ring(5, 0.0), ring(5, 0.01))

    assert len(names) == 5
    for pieces in (first, second):
        extents = [
            np.ptp(np.unwrap(np.arctan2(*curve.point(np.linspace(0, 1, 9)).T[::-1])))
            for _, curve in pieces
        ]
        assert min(extents) > 0.9 * 2 * np.pi / 5


def test_runs_reject_a_station_that_is_not_a_construction_sample():
    """A drifting station would silently take a neighbour's arc."""
    from pyskyfire.regen.cad_geometry import _station_index

    stations = np.array([0.0, 0.5, 1.0])
    assert _station_index(stations, 0.5) == 1
    assert _station_index(stations, 1.0) == 2
    with pytest.raises(ValueError, match="not a construction sample"):
        _station_index(stations, 0.75)


def test_absolute_transition_length_overrides_fraction():
    stations, intervals, _ = hot_wall_stations(
        [(0, 0.5), (0.5, 2)],
        np.linspace(0, 2, 101),
        (0, 2),
        transition_fraction=0.1,
        transition_length_m=0.04,
    )
    np.testing.assert_allclose(intervals, [[0.48, 0.52]])
    assert 0.48 in stations and 0.52 in stations


def test_fraction_uses_full_contour_not_shorter_circuit_extent():
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        _, intervals, _ = hot_wall_stations(
            [(0.2, 0.5), (0.5, 0.8)], [0.2, 0.5, 0.8], (0, 2), transition_fraction=0.1
        )
    np.testing.assert_allclose(intervals, [[0.4, 0.6]])
