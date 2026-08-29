"""Wetted coolant volume of a cooling circuit.

The volume integral has to agree with the flow area the solver divides mass
flow by, follow the channel's own path rather than the axis, and stay correct
when the channels are helical -- where a narrower normal section and a longer
path have to cancel.
"""

import numpy as np
import pytest

from pyskyfire.regen.channel_placement import SurfacePlacement
from pyskyfire.regen.cross_section import ChannelSection
from pyskyfire.regen.thrust_chamber import CoolingCircuit, ThrustChamber


class _StripSection(ChannelSection):
    """Thin rectangular strip: ``A = r * theta * h``, exactly linear in theta.

    Linearity matters for the helical case. It makes ``n * A`` depend on the
    channel count only through the azimuthal share, so the closed-form volume
    below is exact and any error the test sees comes from the integration
    rather than from the section model.
    """

    def A_coolant(self, prof):
        return prof.centerline[:, 1] * prof.theta * prof.h

    def Dh_coolant(self, prof):
        return 2.0 * prof.h

    def P_thermal(self, prof):
        return prof.centerline[:, 1] * prof.theta

    def P_coolant(self, prof):
        return prof.centerline[:, 1] * prof.theta


class _Cylinder:
    """Constant-radius contour, so ``ds/dx`` isolates the helix."""

    def __init__(self, radius=0.1, length=1.0):
        self.radius = float(radius)
        self.xs = np.array([0.0, length])
        self.x_t = 0.0

    def r(self, x):
        return np.full_like(np.asarray(x, dtype=float), self.radius)

    def dr_dx(self, x):
        return np.zeros_like(np.asarray(x, dtype=float))


class _Wall:
    def __init__(self, thickness=0.0):
        self._thickness = float(thickness)

    def thickness(self, x):
        return np.full_like(np.asarray(x, dtype=float), self._thickness)


def _chamber(
    *,
    helix_angle=0.0,
    n_channel_positions=60,
    height=0.002,
    radius=0.1,
    length=1.0,
    n_nodes=200,
):
    contour = _Cylinder(radius=radius, length=length)
    circuit = CoolingCircuit(
        name="test",
        contour=contour,
        cross_section=_StripSection(),
        span=[0.0, 1.0],
        placement=SurfacePlacement(
            n_channel_positions=n_channel_positions,
            helix_angle=helix_angle,
        ),
        channel_height=lambda x: height,
        walls=[_Wall()],
        coolant_transport=None,
        roughness=0.0,
    )
    chamber = ThrustChamber(
        contour=contour,
        cooling_circuits=[circuit],
        combustion_transport=None,
        n_nodes=n_nodes,
        compute_gas=False,
    )
    return chamber, circuit


def _analytic_volume(*, radius, height, length):
    # n channels each of arc width (2 pi cos(gamma) / n) * r and height h, over
    # a path of length L / cos(gamma): the cosines cancel and the circuit holds
    # a full annulus of thickness h. Cross-sections are referenced to the
    # section centreline, which for a surface placement is the hot wall, so
    # ``radius`` here is the wall radius.
    return 2.0 * np.pi * radius * height * length


def test_axial_circuit_matches_the_annulus_it_sweeps() -> None:
    radius, height, length = 0.1, 0.002, 1.0
    chamber, circuit = _chamber(radius=radius, height=height, length=length)

    expected = _analytic_volume(radius=radius, height=height, length=length)
    assert circuit.volume == pytest.approx(expected, rel=1e-12)
    assert chamber.coolant_volume == pytest.approx(expected, rel=1e-12)


def test_volume_is_independent_of_channel_count() -> None:
    # Splitting the same annulus into more channels cannot change how much
    # coolant it holds.
    _, few = _chamber(n_channel_positions=12)
    _, many = _chamber(n_channel_positions=240)

    assert few.volume == pytest.approx(many.volume, rel=1e-12)


@pytest.mark.parametrize("helix_angle", [0.0, 0.3, np.radians(45.0)])
def test_helical_circuit_holds_the_same_volume_but_runs_further(
    helix_angle,
) -> None:
    radius, height, length = 0.1, 0.002, 1.0
    _, axial = _chamber(helix_angle=0.0, radius=radius, height=height, length=length)
    _, helical = _chamber(
        helix_angle=helix_angle, radius=radius, height=height, length=length
    )

    assert helical.volume == pytest.approx(axial.volume, rel=1e-9)

    # The path really is longer -- the volume is unchanged because the normal
    # section narrows by exactly the same factor, not because the helix was
    # ignored.
    path_length = np.trapezoid(helical.ds_dx_vals, helical.x_domain)
    assert path_length == pytest.approx(length / np.cos(helix_angle), rel=1e-6)
    assert helical.is_helical == (helix_angle != 0.0)


def test_volume_converges_with_axial_resolution() -> None:
    # A tapering height makes the integrand non-linear, so the trapezoidal rule
    # is no longer exact and the error has to shrink with refinement.
    radius, length = 0.1, 1.0

    def taper(x):
        return 0.002 * (1.0 + 0.5 * np.sin(3.0 * x))

    volumes = []
    for n_nodes in (25, 50, 100, 200, 400):
        contour = _Cylinder(radius=radius, length=length)
        circuit = CoolingCircuit(
            name="taper",
            contour=contour,
            cross_section=_StripSection(),
            span=[0.0, 1.0],
            placement=SurfacePlacement(n_channel_positions=60),
            channel_height=taper,
            walls=[_Wall()],
            coolant_transport=None,
            roughness=0.0,
        )
        ThrustChamber(
            contour=contour,
            cooling_circuits=[circuit],
            combustion_transport=None,
            n_nodes=n_nodes,
            compute_gas=False,
        )
        volumes.append(circuit.volume)

    x = np.linspace(0.0, length, 200_001)
    reference = 2.0 * np.pi * radius * np.trapezoid(taper(x), x)

    errors = [abs(volume - reference) / reference for volume in volumes]
    assert errors[-1] < 1e-6
    # Trapezoidal error is O(n^-2), so each doubling should cut it by ~4x.
    for coarse, fine in zip(errors, errors[1:]):
        assert fine < coarse / 3.0


def test_stacked_channels_per_leaf_are_counted() -> None:
    _, circuit = _chamber(n_channel_positions=60)
    single_leaf_volume = circuit.volume

    circuit.placement.n_channels_per_leaf = 3
    circuit.compute_volume()

    assert circuit.n_channels == 180
    assert circuit.volume == pytest.approx(3.0 * single_leaf_volume, rel=1e-12)
