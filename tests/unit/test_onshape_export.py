"""Contracts for the advanced example's external Onshape exporter."""

import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from pyskyfire.regen import build_channel_construction
from pyskyfire.regen.cross_section import (
    CrossSectionRounded,
    CrossSectionSquared,
    SectionProfiles,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]
EXPORTER_PATH = PROJECT_ROOT / "examples" / "advanced" / "export_onshape.py"
SPEC = importlib.util.spec_from_file_location("advanced_onshape_export", EXPORTER_PATH)
exporter = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(exporter)


def _circuit(name, count, span, radius, *, cross_section=None):
    axial = np.linspace(span[0], span[1], 3)
    centerline = np.column_stack(
        (axial, np.full_like(axial, radius), np.zeros_like(axial))
    )
    if cross_section is None:
        cross_section = CrossSectionSquared(
            rib_width=lambda pitch_mid: 0.25 * pitch_mid
        )

    class Placement:
        n_channel_positions = count
        occludes = True

        @staticmethod
        def dtheta_dx(x, contour, circuit=None):
            return 0.0

    return SimpleNamespace(
        name=name,
        placement=Placement(),
        reference_path=centerline,
        x_domain=axial,
        channel_heights=np.full_like(axial, 0.004),
        channel_local_sector=np.full_like(axial, 2.0 * np.pi / (2 * count)),
        channel_angular_width=np.full_like(axial, 2.0 * np.pi / (2 * count)),
        channel_wall_thickness=lambda x: np.full_like(
            np.asarray(x, dtype=float),
            0.001,
        ),
        cross_section=cross_section,
    )


def _chamber():
    return SimpleNamespace(
        cooling_circuits=[
            _circuit("Primary", 8, (0.0, 1.0), 0.10),
            _circuit("Return Pass", 8, (0.2, 0.8), 0.12),
        ],
        contour=SimpleNamespace(),
    )


def test_export_defaults_to_analytical_and_uses_reference_origins() -> None:
    data = exporter.build_onshape_data(_chamber())

    assert data["layout"]["method"] == "analytical"
    assert len(data["channels"]) == 2
    assert data["channels"][0]["pattern"] == {
        "count": 8,
        "angularStep_deg": 45.0,
    }

    first_section = data["channels"][0]["sections"][0]
    # Sections are anchored on the channel floor, one wall thickness out.
    assert first_section["plane"]["origin_mm"] == [101.0, 0.0, 0.0]
    assert first_section["plane"]["normal"] == [0.0, 0.0, 1.0]
    assert first_section["plane"]["xAxis"] == [1.0, 0.0, 0.0]
    assert data["channels"][0]["sectionType"] == "squared"
    assert first_section["innerArc"]["center_uv_mm"] == [-101.0, 0.0]
    assert first_section["innerArc"]["radius_mm"] == pytest.approx(101.0)
    assert first_section["outerArc"]["radius_mm"] == pytest.approx(105.0)
    rib = 0.25 * (0.101 + 0.002) * np.radians(22.5)
    expected_half_angle = np.radians(22.5) / 2 - np.arcsin(rib / (2 * 0.101))
    assert first_section["innerArc"]["endAngle_deg"] == pytest.approx(
        np.degrees(expected_half_angle)
    )


def test_export_can_select_one_or_a_subset_by_index_and_name() -> None:
    chamber = _chamber()

    one = exporter.build_onshape_data(chamber, circuit_selectors=["1"])
    subset = exporter.build_onshape_data(
        chamber,
        circuit_selectors=["Return Pass", "0", "Return Pass"],
    )

    assert [channel["circuitIndex"] for channel in one["channels"]] == [1]
    # Layout order remains deterministic even when selector order differs.
    assert [channel["circuitIndex"] for channel in subset["channels"]] == [0, 1]


def test_export_can_use_the_numerical_layout() -> None:
    data = exporter.build_onshape_data(_chamber(), layout="numerical")

    assert data["layout"]["method"] == "numerical"
    assert len(data["channels"]) == 2
    assert all(channel["pattern"]["count"] == 8 for channel in data["channels"])


def test_exporter_writes_valid_json(tmp_path: Path, monkeypatch) -> None:
    input_path = tmp_path / "results.pkl"
    output_path = tmp_path / "channels.json"
    monkeypatch.setattr(
        exporter.psf.common.Results,
        "load",
        lambda path: {"thrust_chamber": _chamber()},
    )

    exporter.export_onshape_json(input_path, output_path)

    document = json.loads(output_path.read_text(encoding="utf-8"))
    assert document["schema"] == "pyskyfire.onshape.cooling-channels"
    assert document["schemaVersion"] == 6


def test_exporter_rejects_cross_sections_the_feature_script_cannot_build() -> None:
    chamber = SimpleNamespace(
        cooling_circuits=[
            _circuit(
                "Rounded",
                8,
                (0.0, 1.0),
                0.1,
                cross_section=object(),
            )
        ],
        contour=SimpleNamespace(),
    )

    with pytest.raises(TypeError, match="Unsupported channel cross section"):
        exporter.build_onshape_data(chamber)


FEATURE_SCRIPT_PATH = (
    PROJECT_ROOT / "examples" / "advanced" / "onshape_cooling_channels.fs"
)


def _sketched_area(inner_radius, height, pitch, rib_width):
    """Area the FeatureScript's sketch encloses, by Green's theorem.

    Mirrors ``createSectorSketch``: an arc at the inner radius, a straight rib
    face out to the outer radius, an arc back, and a straight face in. The
    half-angle differs between the two radii because the rib keeps its width.
    """
    outer_radius = inner_radius + height
    inner_half = pitch / 2.0 - np.arcsin(rib_width / (2.0 * inner_radius))
    outer_half = pitch / 2.0 - np.arcsin(rib_width / (2.0 * outer_radius))

    return abs(
        inner_radius**2 * inner_half
        - outer_radius**2 * outer_half
        + inner_radius * outer_radius * np.sin(outer_half - inner_half)
    )


def test_the_exported_section_encloses_the_area_the_solver_used() -> None:
    """The CAD profile and the section analytics must describe one channel."""
    circuit = _circuit(
        "Primary",
        8,
        (0.0, 1.0),
        0.10,
        cross_section=CrossSectionSquared(rib_width=1.0e-3),
    )
    data = exporter.build_onshape_data(
        SimpleNamespace(cooling_circuits=[circuit], contour=SimpleNamespace())
    )
    section = data["channels"][0]["sections"][0]

    inner_radius = section["innerArc"]["radius_mm"] / 1000.0
    height = (
        section["outerArc"]["radius_mm"] - section["innerArc"]["radius_mm"]
    ) / 1000.0
    rib_width = 1.0e-3
    pitch = circuit.channel_local_sector[0]
    assert np.radians(section["innerArc"]["endAngle_deg"]) == pytest.approx(
        pitch / 2 - np.arcsin(rib_width / (2 * inner_radius))
    )

    profiles = SectionProfiles(
        h=np.array([height]),
        theta=np.array([pitch]),
        t_wall=np.array([0.001]),
        cross_section_base_point_path=np.array([[0.0, inner_radius, 0.0]]),
    )

    assert _sketched_area(inner_radius, height, pitch, rib_width) == pytest.approx(
        circuit.cross_section.A_coolant(profiles)[0], rel=1e-6
    )


def test_the_sketch_and_the_analytics_diverge_only_at_third_order() -> None:
    """The section analytics treat the rib as a chord; the sketch does not.

    The gap between them grows as ``(rib_width / radius)**3``, so it is
    negligible for a real rib and only becomes visible for one comparable to
    the chamber radius.
    """
    inner_radius, height, pitch = 0.101, 0.004, np.radians(22.5)

    def relative_gap(rib_width):
        sketched = _sketched_area(inner_radius, height, pitch, rib_width)
        modelled = (
            0.5 * pitch * ((inner_radius + height) ** 2 - inner_radius**2)
            - rib_width * height
        )
        return abs(sketched - modelled) / modelled

    assert relative_gap(1.0e-3) < 1.0e-6
    assert relative_gap(2.0e-3) / relative_gap(1.0e-3) == pytest.approx(8.0, rel=0.05)


def test_the_sketched_channel_widens_with_radius_while_the_rib_does_not() -> None:
    inner_radius, height, pitch, rib = 0.1, 0.004, np.radians(6.0), 1.0e-3

    inner_half = pitch / 2.0 - np.arcsin(rib / (2.0 * inner_radius))
    outer_half = pitch / 2.0 - np.arcsin(rib / (2.0 * (inner_radius + height)))

    # Wider channel at the roof than at the floor ...
    assert outer_half > inner_half
    # ... yet the rib is the same thickness at both, which is the whole point.
    inner_gap = 2.0 * inner_radius * np.sin(inner_half)
    outer_gap = 2.0 * (inner_radius + height) * np.sin(outer_half)
    chord = 2.0 * inner_radius * np.sin(pitch / 2.0)
    outer_chord = 2.0 * (inner_radius + height) * np.sin(pitch / 2.0)
    assert chord - inner_gap == pytest.approx(outer_chord - outer_gap, rel=1e-3)


def test_the_feature_script_and_the_exporter_agree_on_the_contract() -> None:
    """The FeatureScript cannot run in CI, so pin the interface it reads."""
    script = FEATURE_SCRIPT_PATH.read_text(encoding="utf-8")

    assert "section.innerArc" in script
    assert "section.outerArc" in script
    assert "channel.sectionType" in script
    assert 'sectionType == "squared"' in script
    assert 'sectionType == "rounded"' in script
    assert "!section.circular" in script
    assert "data.schemaVersion != 6" in script
    assert "pitchAngle_deg" not in script
    # The hot wall is lofted curve by curve along each run of constant
    # topology, so the kernel is never left to pair profiles by proximity.
    assert "wall.runs" in script
    assert "run.curveNames" in script
    assert "wall.bridges" in script
    assert "bridge.curveNames" in script
    assert "sketchEntityQuery" in script


def test_mixed_helical_round_and_square_export_matches_native_polygons():
    from pyskyfire.viz.engine_viz import build_sections_for_centerlines
    from pyskyfire.regen.channel_layout import resolve_analytical_layout

    square = _circuit("Squared", 8, (0.0, 1.0), 0.1)
    rounded = _circuit(
        "Rounded", 8, (0.0, 1.0), 0.1, cross_section=CrossSectionRounded()
    )
    for circuit in (square, rounded):
        circuit.reference_path[:, 2] = 0.3 * circuit.x_domain
        circuit.helix_angle_vals = np.full(3, np.arctan(0.03))
        circuit.channel_local_sector *= np.cos(circuit.helix_angle_vals)
        circuit.placement.dtheta_dx = lambda x, contour, circuit=None: 0.3
    rounded.channel_heights[:] = 0.05
    chamber = SimpleNamespace(
        cooling_circuits=[square, rounded], contour=SimpleNamespace()
    )
    data = exporter.build_onshape_data(chamber)
    seeds = resolve_analytical_layout(chamber).channels
    polygons = build_sections_for_centerlines(
        [seed.reference_path for seed in seeds],
        [seed.circuit_index for seed in seeds],
        chamber,
        stride=2,
        n_arc_pts=19,
    )
    assert [channel["sectionType"] for channel in data["channels"]] == [
        "squared",
        "rounded",
    ]
    for exported, native in zip(data["channels"], polygons, strict=True):
        for index, points in zip((0, 2), native, strict=True):
            section = exported["sections"][index]
            arcs = []
            for key in ("innerArc", "outerArc"):
                arc = section[key]
                angles = np.radians(
                    np.linspace(arc["startAngle_deg"], arc["endAngle_deg"], 19)
                )
                arcs.append(
                    np.asarray(arc["center_uv_mm"])
                    + np.cos(angles)[:, None] * np.asarray(arc["axisA_uv_mm"])
                    + np.sin(angles)[:, None] * np.asarray(arc["axisB_uv_mm"])
                )
            uv = np.vstack(arcs)
            normal = np.asarray(section["plane"]["normal"])
            radial = np.asarray(section["plane"]["xAxis"])
            world = (
                np.asarray(section["plane"]["origin_mm"])
                + uv[:, :1] * radial
                + uv[:, 1:] * np.cross(normal, radial)
            )
            # Onshape XYZ corresponds to native YZX.
            np.testing.assert_allclose(world / 1000, points[:, [1, 2, 0]], atol=1e-15)
    assert data["channels"][0]["sections"][0]["plane"]["normal"] == [0, 0, 1]
    assert data["channels"][1]["sections"][0]["plane"]["normal"] == [0, 0, 1]


@pytest.mark.parametrize("base_height", [0.0, 0.004])
def test_rounded_base_and_full_height_inputs_export_the_same_arcs(base_height):
    base = _circuit(
        "Rounded",
        8,
        (0.0, 1.0),
        0.1,
        cross_section=CrossSectionRounded(height_input="base"),
    )
    base.channel_heights[:] = base_height
    full = _circuit("Rounded", 8, (0.0, 1.0), 0.1, cross_section=CrossSectionRounded())
    # Lobes are sized where they end up, so the equivalent full height follows
    # from the solved lobe-centre distance rather than from the contour radius.
    profiles = build_channel_construction(base).profiles
    half_theta = profiles.theta / 2
    full.channel_heights[:] = 2 * (
        profiles.sector_radius * np.sin(half_theta) - 0.001
    ) + base_height * (1 + np.sin(half_theta))
    documents = [
        exporter.build_onshape_data(
            SimpleNamespace(cooling_circuits=[circuit], contour=SimpleNamespace())
        )
        for circuit in (base, full)
    ]
    for base_section, full_section in zip(
        documents[0]["channels"][0]["sections"],
        documents[1]["channels"][0]["sections"],
        strict=True,
    ):
        assert base_section["circular"] == (base_height == 0.0)
        assert base_section["plane"] == full_section["plane"]
        for key in ("innerArc", "outerArc"):
            for field in (
                "center_uv_mm",
                "radius_mm",
                "startAngle_deg",
                "endAngle_deg",
            ):
                np.testing.assert_allclose(
                    base_section[key][field], full_section[key][field], atol=1e-12
                )


def test_export_optional_hot_wall_and_feature_toggle():
    chamber = _chamber()
    chamber.contour.xs = np.array([0.0, 1.0])
    for circuit in chamber.cooling_circuits:
        circuit.reference_path[:, 1] = 0.1
    assert "hotWall" not in exporter.build_onshape_data(chamber)
    data = exporter.build_onshape_data(
        chamber, include_hot_wall=True, hot_wall_transition_fraction=0.1
    )
    assert data["schemaVersion"] == 6
    np.testing.assert_allclose(
        data["hotWall"]["transitionIntervals_mm"], [[150, 250], [750, 850]]
    )
    assert all(
        curve["kind"] == "circle"
        for s in data["hotWall"]["sections"]
        for curve in s["curves"]
    )
    # The return pass starts and ends inside the chamber pass, so the wall is
    # three runs of constant topology with a hand-off bridging each seam.
    assert [run["sectionIndices"] for run in data["hotWall"]["runs"]] == [
        [0, 1],
        [2, 3, 4],
        [5, 6],
    ]
    script = FEATURE_SCRIPT_PATH.read_text()
    assert "definition.constructHotWall is boolean" in script
    assert "ToolBodyType.SURFACE" in script
    assert "skEllipticalArc" in script and "skBezier" in script
