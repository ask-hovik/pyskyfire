"""Integrity checks for the digitized RL10 reference data."""

import json
from pathlib import Path

import numpy as np
import pytest


RL10_DIR = Path(__file__).resolve().parents[2] / "validation" / "RL10"
REFERENCE_DIR = RL10_DIR / "reference_data"


@pytest.mark.parametrize(
    (
        "json_name",
        "value_key",
        "value_scale",
        "source_page",
        "point_count",
        "first_source_point",
        "last_source_point",
    ),
    [
        (
            "reference_coolant_static_pressure.json",
            "p_static",
            6894.757293168,
            "p. 56",
            20,
            (-13.193762683145334, 864.003194220636),
            (20.788938150614946, 1075.073630220304),
        ),
        (
            "reference_coolant_static_temperature.json",
            "T_static",
            5.0 / 9.0,
            "p. 56",
            20,
            (-13.112653883168392, 366.63940855459697),
            (20.734912078681376, 80.53270139016371),
        ),
        (
            "reference_heat_flux.json",
            "dQ_dA",
            1055.05585262 / 0.0254**2,
            "p. 57",
            19,
            (-11.796934449396666, 6.34112403522122),
            (16.710620719643444, 1.560560930535928),
        ),
    ],
)
def test_binder_reference_data_contains_converted_source_points(
    json_name,
    value_key,
    value_scale,
    source_page,
    point_count,
    first_source_point,
    last_source_point,
) -> None:
    reference = json.loads((REFERENCE_DIR / json_name).read_text())[
        "Binder_et_al"
    ]

    assert len(reference["x"]) == point_count
    assert len(reference[value_key]) == point_count
    assert reference["x"][0] == pytest.approx(first_source_point[0] * 0.0254)
    assert reference[value_key][0] == pytest.approx(
        first_source_point[1] * value_scale
    )
    assert reference["x"][-1] == pytest.approx(last_source_point[0] * 0.0254)
    assert reference[value_key][-1] == pytest.approx(
        last_source_point[1] * value_scale
    )
    assert reference["name"] == "Binder et. al."
    assert source_page in reference["source"]
    assert "NASA TM-107318" in reference["source"]


def test_all_rl10_reference_json_files_declare_units() -> None:
    reference_files = sorted(REFERENCE_DIR.glob("reference_*.json"))

    assert reference_files
    for path in reference_files:
        data = json.loads(path.read_text())
        assert isinstance(data.get("units"), dict), path.name
        assert data["units"], path.name


@pytest.mark.parametrize(
    "json_name",
    [
        "reference_wall_temperature.json",
        "reference_coolant_static_pressure.json",
        "reference_coolant_static_temperature.json",
    ],
)
def test_design_report_references_figure_d3_source(json_name) -> None:
    reference = json.loads((REFERENCE_DIR / json_name).read_text())[
        "RL10_design_report"
    ]

    assert "PWA FR-1769" in reference["source"]
    assert "Figure D-3" in reference["source"]
    assert "p. D-4" in reference["source"]
    assert reference["source"].endswith(
        "https://www.nasa.gov/wp-content/uploads/2025/06/"
        "design-report-for-rl-10-a-3-3-1966.pdf"
    )


@pytest.mark.parametrize(
    ("json_name", "series_name", "point_count"),
    [
        ("reference_channel_height.json", "design_points", 32),
        (
            "reference_channel_width.json",
            "pre_interleaving_180_tube_region",
            14,
        ),
        (
            "reference_channel_width.json",
            "post_interleaving_360_tube_region",
            7,
        ),
    ],
)
def test_recop_tube_dimensions_retain_source_and_si_values(
    json_name,
    series_name,
    point_count,
) -> None:
    data = json.loads((REFERENCE_DIR / json_name).read_text())
    series = data["series"][series_name]
    value_name = (
        "tube_height_od" if "height" in json_name else "tube_width_od"
    )

    assert len(series["x_in"]) == point_count
    assert len(series["x"]) == point_count
    assert len(series[f"{value_name}_in"]) == point_count
    assert len(series[value_name]) == point_count
    np.testing.assert_allclose(
        series["x"],
        np.asarray(series["x_in"]) * 0.0254,
        rtol=0.0,
        atol=1.0e-15,
    )
    np.testing.assert_allclose(
        series[value_name],
        np.asarray(series[f"{value_name}_in"]) * 0.0254,
        rtol=0.0,
        atol=1.0e-15,
    )
    assert data["source"]["url"] == (
        "https://ntrs.nasa.gov/api/citations/19950002773/"
        "downloads/19950002773.pdf"
    )
    assert data["source"]["figure"].startswith("Figure 10")


@pytest.mark.parametrize(
    ("curve_name", "point_count"),
    [("inner", 5), ("outer", 10)],
)
def test_silver_insert_reference_curve_uses_si_coordinates(
    curve_name,
    point_count,
) -> None:
    data = json.loads(
        (REFERENCE_DIR / "reference_silver_insert.json").read_text()
    )
    curve = data["curves"][curve_name]

    assert len(curve["x"]) == point_count
    assert len(curve["r"]) == point_count
    assert curve["x"][0] == pytest.approx(-0.990122656102788 * 0.0254)
    assert curve["r"][0] == pytest.approx(2.66869618184129 * 0.0254)
    assert curve["x"][-1] == pytest.approx(0.971314510580677 * 0.0254)
    assert curve["r"][-1] == pytest.approx(3.05799236631186 * 0.0254)
    assert data["units"] == {"x": "m", "r": "m"}
    assert "p. 55" in data["source"]


def test_silver_insert_curves_share_endpoints() -> None:
    data = json.loads(
        (REFERENCE_DIR / "reference_silver_insert.json").read_text()
    )
    inner = data["curves"]["inner"]
    outer = data["curves"]["outer"]

    np.testing.assert_allclose(
        [inner["x"][0], inner["r"][0]],
        [outer["x"][0], outer["r"][0]],
    )
    np.testing.assert_allclose(
        [inner["x"][-1], inner["r"][-1]],
        [outer["x"][-1], outer["r"][-1]],
    )
