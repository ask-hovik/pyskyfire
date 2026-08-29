"""Spline-backed thrust-chamber contours and their plotting behavior."""

import numpy as np
import pytest

from pyskyfire.regen.contour import Contour
from pyskyfire.viz.plot_regen import PlotContour


def test_contour_leaves_both_throat_slopes_unconstrained_by_default() -> None:
    xs = np.array([-2.0, -1.0, 0.0, 0.4, 1.0, 2.0])
    rs = np.array([2.0, 1.4, 1.0, 1.1, 1.5, 2.4])

    contour = Contour(xs, rs)

    assert contour.x_t == 0.0
    assert contour.r_t == pytest.approx(1.0)
    assert contour.chamber_angle is None
    assert contour.nozzle_angle is None
    assert contour._converging_spline(0.0, 1) != pytest.approx(0.0)
    assert contour._diverging_spline(0.0, 1) != pytest.approx(0.0)
    np.testing.assert_allclose(contour.r(xs), rs)


def test_contour_clamps_splines_to_independent_throat_angles() -> None:
    xs = np.array([-2.0, -1.0, 0.0, 0.4, 1.0, 2.0])
    rs = np.array([2.0, 1.4, 1.0, 1.1, 1.5, 2.4])
    chamber_angle = np.radians(12.0)
    nozzle_angle = np.radians(18.0)

    contour = Contour(
        xs,
        rs,
        chamber_angle=chamber_angle,
        nozzle_angle=nozzle_angle,
    )

    assert contour.chamber_angle == pytest.approx(chamber_angle)
    assert contour.nozzle_angle == pytest.approx(nozzle_angle)
    assert contour._converging_spline(0.0, 1) == pytest.approx(
        -np.tan(chamber_angle)
    )
    assert contour._diverging_spline(0.0, 1) == pytest.approx(
        np.tan(nozzle_angle)
    )
    assert contour.dr_dx(0.0) == pytest.approx(-np.tan(chamber_angle))
    np.testing.assert_allclose(contour.r(xs), rs)


def test_contour_preserves_leading_cylindrical_chamber() -> None:
    xs = np.array([-3.0, -2.0, -1.0, 0.0, 0.5, 1.0])
    rs = np.array([2.0, 2.0, 1.6, 1.0, 1.2, 1.8])
    contour = Contour(xs, rs)
    chamber_x = np.linspace(-3.0, -2.0, 101)

    np.testing.assert_allclose(contour.r(chamber_x), 2.0)
    np.testing.assert_allclose(contour.dr_dx(chamber_x), 0.0)


def test_contour_interpolates_missing_throat_from_three_closest_points() -> None:
    xs = np.array([-2.0, -0.5, 0.25, 1.0, 2.0])
    rs = 1.0 + xs**2

    with pytest.warns(UserWarning, match="throat radius is not defined"):
        contour = Contour(xs, rs)

    np.testing.assert_array_equal(contour.input_xs, xs)
    np.testing.assert_array_equal(contour.input_rs, rs)
    assert np.count_nonzero(contour.xs == 0.0) == 1
    assert contour.r_t == pytest.approx(1.0)
    assert np.isfinite(contour.dr_dx(0.0))


def test_contour_applies_angles_to_an_inserted_throat() -> None:
    xs = np.array([-2.0, -0.5, 0.25, 1.0, 2.0])
    rs = 1.0 + xs**2
    chamber_angle = np.radians(9.0)
    nozzle_angle = np.radians(14.0)

    with pytest.warns(UserWarning, match="throat radius is not defined"):
        contour = Contour(
            xs,
            rs,
            chamber_angle=chamber_angle,
            nozzle_angle=nozzle_angle,
        )

    assert np.count_nonzero(contour.xs == 0.0) == 1
    assert contour.r_t == pytest.approx(1.0)
    assert contour._converging_spline(0.0, 1) == pytest.approx(
        -np.tan(chamber_angle)
    )
    assert contour._diverging_spline(0.0, 1) == pytest.approx(
        np.tan(nozzle_angle)
    )


@pytest.mark.parametrize(
    ("parameter", "value"),
    [
        ("chamber_angle", -0.1),
        ("nozzle_angle", np.pi / 2.0),
        ("chamber_angle", np.nan),
    ],
)
def test_contour_rejects_invalid_throat_angles(parameter, value) -> None:
    xs = np.array([-1.0, 0.0, 1.0])
    rs = np.array([2.0, 1.0, 2.0])

    with pytest.raises(ValueError, match=parameter):
        Contour(xs, rs, **{parameter: value})


def test_plot_contour_samples_spline_and_can_show_input_points() -> None:
    xs = np.array([-2.0, -1.0, 0.0, 0.5, 1.0])
    rs = np.array([2.0, 1.4, 1.0, 1.2, 1.8])
    contour = Contour(xs, rs, name="test contour")

    plot = PlotContour(
        contour,
        show_input_points=True,
        interpolation_points=101,
    )

    assert len(plot.fig.data) == 4
    assert len(plot.fig.data[0].x) > len(xs)
    np.testing.assert_allclose(plot.fig.data[0].y, contour.r(plot.fig.data[0].x))
    np.testing.assert_array_equal(plot.fig.data[2].x, xs)
    np.testing.assert_array_equal(plot.fig.data[2].y, rs)
    assert plot.fig.data[2].marker.symbol == "x"
