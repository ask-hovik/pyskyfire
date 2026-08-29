"""Thrust-chamber assembly behavior not tied to a numerical solve."""

from types import SimpleNamespace

import numpy as np
import pytest

from pyskyfire.regen import FilmCooling, ThrustChamber


class _Contour:
    xs = np.array([-2.0, 0.0, 4.0])
    x_t = 0.0


def _empty_chamber(*, film_cooling=None):
    return ThrustChamber(
        contour=_Contour(),
        cooling_circuits=[],
        combustion_transport=SimpleNamespace(),
        compute_gas=False,
        film_cooling=film_cooling,
    )


@pytest.mark.parametrize(
    ("fraction", "expected"),
    [(-1.0, -2.0), (-0.5, -1.0), (0.0, 0.0), (0.5, 2.0), (1.0, 4.0)],
)
def test_signed_film_injection_fraction_maps_about_the_throat(fraction, expected) -> None:
    film = FilmCooling(
        coolant_transport=SimpleNamespace(),
        x_fraction=fraction,
        coolant_mass_flow_rate=0.1,
        film_injection_perimeter=0.2,
        liquid_absorptivity=0.5,
        mole_fraction_H2O=0.2,
        mole_fraction_CO2=0.1,
    )

    _empty_chamber(film_cooling=film)

    assert film.x == pytest.approx(expected)


def test_signed_film_injection_fraction_rejects_locations_outside_contour() -> None:
    film = FilmCooling(
        coolant_transport=SimpleNamespace(),
        x_fraction=1.1,
        coolant_mass_flow_rate=0.1,
        film_injection_perimeter=0.2,
        liquid_absorptivity=0.5,
        mole_fraction_H2O=0.2,
        mole_fraction_CO2=0.1,
    )

    with pytest.raises(ValueError, match="between -1 and 1"):
        _empty_chamber(film_cooling=film)
