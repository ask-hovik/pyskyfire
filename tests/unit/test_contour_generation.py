"""High-level thrust-chamber contour contracts."""

import numpy as np
import pytest

from pyskyfire.regen.contour import (
    compute_chamber_volume,
    compute_cutoff_length,
    get_contour,
    integrate_area,
)


@pytest.mark.parametrize("nozzle", ["conical", "rao"])
def test_generated_contour_hits_throat_and_exit_area_ratio(nozzle) -> None:
    throat_radius = 0.05
    area_ratio = 10.0
    x, radius = get_contour(
        r_t=throat_radius,
        area_ratio=area_ratio,
        r_c=0.1,
        L_c=0.2,
        nozzle=nozzle,
    )

    assert np.all(np.diff(x) > 0.0)
    assert np.all(radius > 0.0)
    assert x[np.argmin(radius)] == pytest.approx(0.0, abs=1e-12)
    assert np.min(radius) == pytest.approx(throat_radius)
    assert (radius[-1] / throat_radius) ** 2 == pytest.approx(area_ratio)
    assert compute_chamber_volume(x, radius) > 0.0


def test_contour_volume_mode_reproduces_requested_chamber_volume() -> None:
    target_volume = 0.005
    x, radius = get_contour(
        r_t=0.05,
        area_ratio=8.0,
        r_c=0.1,
        V_c=target_volume,
        nozzle="conical",
    )

    assert compute_chamber_volume(x, radius) == pytest.approx(
        target_volume, rel=2e-4
    )


def test_cutoff_length_and_integrated_area_have_closed_forms_on_a_cylinder() -> None:
    x = np.linspace(-2.0, 0.0, 201)
    radius = np.full_like(x, 0.5)
    target_length = 0.8
    target_volume = np.pi * 0.5**2 * target_length

    assert compute_cutoff_length(target_volume, x, radius) == pytest.approx(target_length)
    assert integrate_area(target_length, x, radius) == pytest.approx(0.5 * target_length)


def test_contour_rejects_ambiguous_input_modes() -> None:
    with pytest.raises(ValueError, match="Invalid input combination"):
        get_contour(r_t=0.05, area_ratio=8.0, r_c=0.1)
