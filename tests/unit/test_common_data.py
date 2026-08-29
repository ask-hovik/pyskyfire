"""Contracts for common fluid descriptors and result containers."""

import numpy as np
import pytest

from pyskyfire.common import Fluid, Results


def test_mass_fractions_are_converted_to_mole_fractions() -> None:
    fluid = Fluid(
        "coolant",
        ["Water", "Ethanol"],
        [0.5, 0.5],
        basis="mass",
        precision=6,
    )

    molar_masses = np.asarray(fluid.molar_masses())
    expected = (np.array([0.5, 0.5]) / molar_masses)
    expected /= expected.sum()

    np.testing.assert_allclose(fluid.as_mole_fractions(), expected, atol=5e-7)
    assert fluid.coolprop_string().startswith("HEOS::Water[")
    assert "&Ethanol[" in fluid.coolprop_string()


def test_mole_fraction_export_normalizes_input_and_rejects_unknown_basis() -> None:
    fluid = Fluid("coolant", ["Water", "Ethanol"], [2.0, 1.0], basis="mole")
    assert fluid.as_mole_fractions() == pytest.approx([2.0 / 3.0, 1.0 / 3.0], abs=5e-4)

    fluid.basis = "volume"
    with pytest.raises(ValueError, match="mass.*mole"):
        fluid.as_mole_fractions()


def test_results_mapping_attribute_access_and_round_trip(tmp_path) -> None:
    results = Results()
    results["pressure"] = 3.2e6
    results.temperature = np.array([100.0, 200.0])
    results.metadata["case"] = "manufactured"

    path = tmp_path / "results.pkl"
    results.save(path)
    restored = Results.load(path)

    assert restored["pressure"] == pytest.approx(3.2e6)
    np.testing.assert_allclose(restored["temperature"], [100.0, 200.0])
    assert restored.pressure == restored["pressure"]
    assert restored.metadata == {"case": "manufactured"}
    assert restored.list_items() == ["pressure", "temperature"]


def test_results_missing_attributes_follow_python_attribute_protocol() -> None:
    with pytest.raises(AttributeError, match="missing"):
        _ = Results().missing
