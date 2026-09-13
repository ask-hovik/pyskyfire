"""The representative documentation chart uses the fluids backend."""

import fluids
import numpy as np
import pytest

from pyskyfire.viz.lookup_tables import PlotMoodyDiagram


def _straight_factor(reynolds, relative_roughness):
    critical = 2_040.0
    end = critical * 3_500.0 / 2_300.0
    fraction = np.clip((reynolds - critical) / (end - critical), 0.0, 1.0)
    weight = fraction**2 * (3.0 - 2.0 * fraction)
    laminar = fluids.friction.friction_laminar(reynolds)
    if weight == 0.0:
        return laminar
    turbulent = fluids.friction.Clamond(reynolds, relative_roughness)
    return (1.0 - weight) * laminar + weight * turbulent


@pytest.mark.parametrize("limits", [(700.0, 1e8), (700.0, 1500.0), (3000.0, 1e8)])
def test_moody_curves_match_independent_smooth_blend_without_gaps(limits):
    roughnesses = [0.0, 1e-4, 0.01]
    plot = PlotMoodyDiagram(
        rel_rough_list=roughnesses, Re_min=limits[0], Re_max=limits[1], n_pts=40,
    )
    assert len(plot.fig.data) == len(roughnesses)
    for trace, roughness, annotation in zip(
        plot.fig.data, roughnesses, plot.fig.layout.annotations,
    ):
        reynolds = np.asarray(trace.x, dtype=float)
        factors = np.asarray(trace.y, dtype=float)
        assert np.isfinite(reynolds).all()
        assert np.isfinite(factors).all()
        expected = [_straight_factor(Re, roughness) for Re in reynolds]
        np.testing.assert_allclose(factors, expected, rtol=1e-12)
        assert annotation.x == pytest.approx(np.log10(reynolds[-1]))
        assert annotation.y == pytest.approx(np.log10(factors[-1]))
