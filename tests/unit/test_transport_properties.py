"""Transport-property dispatch and CoolProp wrapper contracts."""

import CoolProp.CoolProp as CP
import pytest

from pyskyfire.common import Fluid
from pyskyfire.skycea.coolant_transport import CoolantTransport, TransportProperties


def test_transport_properties_support_constants_and_callables() -> None:
    properties = TransportProperties(
        Pr=0.8,
        mu=lambda T, p: T / p,
        k=0.2,
        cp=lambda T, p: T + p * 1e-6,
        rho=10.0,
        gamma_coolant=1.4,
    )

    assert properties.Pr(300.0, 1e5) == pytest.approx(0.8)
    assert properties.mu(300.0, 1e5) == pytest.approx(0.003)
    assert properties.cp(300.0, 1e5) == pytest.approx(300.1)
    assert properties.rho(300.0, 1e5) == pytest.approx(10.0)
    assert properties.gamma_coolant(300.0, 1e5) == pytest.approx(1.4)
    assert properties.compressible_coolant


def test_transport_properties_require_coolant_specific_fields() -> None:
    properties = TransportProperties(Pr=1.0, mu=1.0, k=1.0)

    with pytest.raises(ValueError, match="density"):
        properties.rho(300.0, 1e5)
    with pytest.raises(ValueError, match="specific heat"):
        properties.cp(300.0, 1e5)
    with pytest.raises(ValueError, match="gamma"):
        properties.gamma_coolant(300.0, 1e5)


def test_coolant_transport_matches_direct_coolprop_queries() -> None:
    transport = CoolantTransport(Fluid("coolant", ["Water"], [1.0]))
    temperature, pressure = 300.0, 2.0e5

    assert transport.get_rho(temperature, pressure) == pytest.approx(
        CP.PropsSI("DMASS", "T", temperature, "P", pressure, transport.fluid)
    )
    assert transport.get_cp(temperature, pressure) > transport.get_cv(
        temperature, pressure
    )
    assert transport.get_mu(temperature, pressure) > 0.0
    assert transport.get_k(temperature, pressure) > 0.0
