"""Generic solid-property model behavior independent of any material table."""

import numpy as np
import pytest

from pyskyfire.common.solids import (
    ConstantModel,
    Log10PolynomialModel,
    Material,
    PiecewiseModel,
    PolynomialModel,
    TabulatedModel,
)


def test_constant_polynomial_and_tabulated_models_vectorize() -> None:
    temperature = np.array([1.0, 2.0, 3.0])

    np.testing.assert_allclose(ConstantModel(7.0)(temperature), 7.0)
    np.testing.assert_allclose(
        PolynomialModel([1.0, 2.0, 3.0])(temperature),
        1.0 + 2.0 * temperature + 3.0 * temperature**2,
    )
    np.testing.assert_allclose(
        TabulatedModel([0.0, 10.0], [2.0, 12.0])([2.5, 7.5]),
        [4.5, 9.5],
    )


def test_property_models_validate_domains_and_tables() -> None:
    bounded = PolynomialModel([1.0], Tmin=10.0, Tmax=20.0, enforce_bounds=True)
    with pytest.raises(ValueError, match="below Tmin"):
        bounded(9.0)
    with pytest.raises(ValueError, match="above Tmax"):
        bounded(21.0)
    with pytest.raises(ValueError, match="T > 0"):
        Log10PolynomialModel([1.0])(0.0)
    with pytest.raises(ValueError, match="strictly increasing"):
        TabulatedModel([0.0, 0.0], [1.0, 2.0])


def test_piecewise_model_interpolates_gaps_averages_overlaps_and_clips() -> None:
    model = PiecewiseModel(
        [
            (0.0, 2.0, ConstantModel(2.0)),
            (4.0, 6.0, ConstantModel(6.0)),
            (5.0, 8.0, ConstantModel(10.0)),
        ],
        range_policy="clip",
    )

    assert model(3.0) == pytest.approx(4.0)
    assert model(5.5) == pytest.approx(8.0)
    assert model(-1.0) == pytest.approx(2.0)
    assert model(9.0) == pytest.approx(10.0)


def test_material_reports_missing_properties() -> None:
    material = Material("conductive only", k=ConstantModel(12.0))

    assert material.get_k(300.0) == pytest.approx(12.0)
    with pytest.raises(AttributeError, match="density not set"):
        material.get_rho(300.0)
