"""Pump losses charged component by component: disk friction and the pump block."""
import CoolProp.CoolProp as CP
import numpy as np
import pytest

from pyskyfire.network import PumpBlock, Station
from pyskyfire.pump import DiskFriction, disk_friction_power


def test_disk_friction_follows_gulich_table_3_6_by_hand():
    rho, mu, n, R, s, Rn = 1000., 1e-3, 3000., .15, .004, .04
    omega = 2*np.pi*n/60
    Re = omega*R**2*rho/mu
    k = np.pi*R/(2*Re*s)+.02/Re**.2*(1+s/R)/(1+s/(2*R))
    expected = k*rho*omega**3*R**5*(1-(Rn/R)**5)
    assert disk_friction_power(rho=rho, mu=mu, n=n, radius=R, axial_gap=s,
                               inner_radius=Rn) == pytest.approx(expected)
    # A rough face costs more than a smooth one.
    assert disk_friction_power(rho=rho, mu=mu, n=n, radius=R, axial_gap=s, inner_radius=Rn,
                               roughness=20e-6) > expected
    # Two faces are the sum of each.
    pair = DiskFriction(radius=R, axial_gap=s, inner_radii=(Rn, 0.))
    single = [disk_friction_power(rho=rho, mu=mu, n=n, radius=R, axial_gap=s, inner_radius=r)
              for r in (Rn, 0.)]
    assert pair.power(rho=rho, mu=mu, n=n) == pytest.approx(sum(single))


@pytest.mark.parametrize("changes, message", [
    (dict(axial_gap=0), "axial_gap"), (dict(inner_radius=.2), "inner_radius"),
    (dict(roughness=-1), "roughness"), (dict(cone_angle_deg=50), "45"),
    (dict(mu=1e3), "Re > 10"),
])
def test_disk_friction_inputs_are_checked(changes, message):
    args = dict(rho=1000., mu=1e-3, n=3000., radius=.15, axial_gap=.004)
    args.update(changes)
    with pytest.raises(ValueError, match=message):
        disk_friction_power(**args)


def test_pump_block_charges_each_loss_once():
    disk = DiskFriction(radius=.1, axial_gap=.003, inner_radii=(.03, .015))
    pump = PumpBlock("p", "in", "out", ["pipe"], 2e5, 2e5, 1., 6000., .9, "Water",
                     mechanical_efficiency=.97, disk_friction=disk)
    stations = {"in": Station(2e5, 300., 10.)}
    signals = {"dp_pipe": 5e6}
    out, power = pump.compute(stations, signals)
    rho = CP.PropsSI("Dmass", "T", 300.-1e-3, "P", 2e5, "Water")
    mu = CP.PropsSI("V", "T", 300.-1e-3, "P", 2e5, "Water")
    blade = 10.*5e6/rho/.9
    friction = disk.power(rho=rho, mu=mu, n=6000.)
    assert out["out"].p == pytest.approx(52e5)
    assert power["P_p"] == pytest.approx((blade+friction)/.97)
    # Blade work and disk friction heat the fluid; the bearings do not.
    h = CP.PropsSI("Hmass", "T", 300.-1e-3, "P", 2e5, "Water")+(blade+friction)/10.
    assert out["out"].T == pytest.approx(CP.PropsSI("T", "P", 52e5, "Hmass", h, "Water"))
    split = pump.post_process(stations, signals)
    assert split["blade_power"] == pytest.approx(blade)
    assert split["disk_friction_power"] == pytest.approx(friction)
    assert split["shaft_power"] == pytest.approx(power["P_p"])
    assert split["mechanical_loss"] == pytest.approx(power["P_p"]*.03)
    # Without disk friction or bearings the blade work is the whole shaft power.
    bare = PumpBlock("p", "in", "out", ["pipe"], 2e5, 2e5, 1., 6000., .9, "Water")
    assert bare.compute(stations, signals)[1]["P_p"] == pytest.approx(blade)
    with pytest.raises(ValueError, match="hydraulic_efficiency"):
        PumpBlock("p", "in", "out", [], 0., 0., 1., 6000., 1.2, "Water")
