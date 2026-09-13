"""Design-input functions as they are actually reached through a chamber.

The unit tests cover dispatch in isolation. These cover the wiring: that each
slot is handed a context containing the quantities its tier promises, and that
the promise holds for quantities the circuit cannot compute on its own, such as
the cross-circuit channel count.
"""

import numpy as np
import pytest

import pyskyfire as psf
from pyskyfire.regen.design_inputs import DesignInputError


class _Cylinder:
    """Constant-radius contour, so every tier A quantity is easy to predict."""

    def __init__(self, radius=0.1, length=1.0):
        self.radius = float(radius)
        self.xs = np.array([0.0, length])
        self.x_t = 0.0
        self.r_t = float(radius)
        self.r_c = float(radius)
        self.r_e = float(radius)

    def r(self, x):
        return np.full_like(np.asarray(x, dtype=float), self.radius)

    def dr_dx(self, x):
        return np.zeros_like(np.asarray(x, dtype=float))

    def normal_angle(self, x):
        return np.zeros_like(np.asarray(x, dtype=float))


class _Material:
    @staticmethod
    def get_k(T):
        return 350.0


def _chamber(
    *,
    rib_width=5.0e-4,
    channel_height=0.004,
    helix_angle=0.0,
    wall_thickness=0.001,
    roughness=0.0,
    n_channel_positions=60,
    radius=0.1,
):
    contour = _Cylinder(radius=radius)
    circuit = psf.regen.CoolingCircuit(
        name="test",
        contour=contour,
        cross_section=psf.regen.CrossSectionSquared(rib_width=rib_width),
        span=[0.0, 1.0],
        placement=psf.regen.SurfacePlacement(
            n_channel_positions=n_channel_positions,
            helix_angle=helix_angle,
        ),
        channel_height=channel_height,
        walls=[psf.regen.Wall(_Material(), wall_thickness)],
        coolant_transport=None,
        roughness=roughness,
    )
    chamber = psf.regen.ThrustChamber(
        contour=contour,
        cooling_circuits=[circuit],
        combustion_transport=None,
        n_nodes=20,
        compute_gas=False,
    )
    return chamber, circuit


# --- tier A -----------------------------------------------------------------


def test_a_wall_thickness_can_follow_the_contour() -> None:
    _, circuit = _chamber(wall_thickness=lambda r: 0.01 * r)

    assert circuit.channel_wall_thickness(0.5) == pytest.approx(0.001)


def test_a_channel_height_that_ignores_position_needs_no_parameters() -> None:
    _, constant = _chamber(channel_height=lambda: 0.004)
    _, plain = _chamber(channel_height=0.004)

    np.testing.assert_allclose(constant.channel_heights, plain.channel_heights)


def test_a_helix_angle_that_ignores_position_needs_no_parameters() -> None:
    _, circuit = _chamber(helix_angle=lambda: np.radians(30.0))

    np.testing.assert_allclose(
        circuit.helix_angle_vals, np.radians(30.0), rtol=1e-12
    )
    assert circuit.is_helical


# --- tier B -----------------------------------------------------------------


def test_a_channel_height_can_depend_on_the_wall_stack_and_channel_count() -> None:
    _, circuit = _chamber(
        channel_height=lambda t_wall, n_channels: 2.0 * t_wall + 1.0e-6 * n_channels,
        wall_thickness=0.001,
        n_channel_positions=60,
    )

    np.testing.assert_allclose(circuit.channel_heights, 0.002 + 60.0e-6)


def test_a_rib_width_can_depend_on_the_cross_circuit_channel_count() -> None:
    _, circuit = _chamber(
        rib_width=lambda circumference, n_channels_all_circuits: (
            0.2 * circumference / n_channels_all_circuits
        ),
        n_channel_positions=60,
    )

    expected = 0.2 * 2.0 * np.pi * 0.1 / 60.0
    widths = circuit.cross_section.rib_width(
        circuit._prof(circuit.cross_section_base_point_path)
    )
    np.testing.assert_allclose(widths, expected, rtol=1e-12)


def test_a_cross_circuit_count_needs_a_chamber() -> None:
    circuit = psf.regen.CoolingCircuit(
        name="orphan",
        contour=_Cylinder(),
        cross_section=psf.regen.CrossSectionSquared(rib_width=1e-4),
        span=[0.0, 1.0],
        placement=psf.regen.SurfacePlacement(n_channel_positions=60),
        channel_height=lambda n_channels_all_circuits: 1e-3,
        walls=[psf.regen.Wall(_Material(), 0.001)],
        coolant_transport=None,
        roughness=0.0,
    )

    with pytest.raises(DesignInputError, match="belongs to a ThrustChamber"):
        circuit.channel_height(0.5)


# --- tier C -----------------------------------------------------------------


def test_a_rib_width_proportional_to_pitch_reproduces_constant_blockage() -> None:
    """The migration path off the old ``blockage_ratio`` is exact."""
    fraction = 0.1
    _, circuit = _chamber(rib_width=lambda pitch_mid: fraction * pitch_mid)
    _, open_channel = _chamber(rib_width=0.0)

    # Not bit-exact: the section metrics integrate the real straight rib faces
    # rather than treating the rib as a wedge of the same angular width.
    np.testing.assert_allclose(
        circuit.A_coolant_vals,
        (1.0 - fraction) * open_channel.A_coolant_vals,
        rtol=1e-6,
    )


def test_a_constant_rib_leaves_a_channel_that_widens_with_radius() -> None:
    rib = 5.0e-4
    _, circuit = _chamber(rib_width=rib, channel_height=0.004)
    _, open_channel = _chamber(rib_width=0.0, channel_height=0.004)

    np.testing.assert_allclose(
        open_channel.A_coolant_vals - circuit.A_coolant_vals,
        rib * 0.004,
        rtol=1e-5,
    )


def test_a_roughness_can_depend_on_the_resolved_section() -> None:
    _, circuit = _chamber(roughness=lambda h: 0.01 * h)

    assert circuit.roughness(0.5) == pytest.approx(0.01 * 0.004)


def test_the_helix_angle_narrows_the_pitch_a_rib_width_sees() -> None:
    _, axial = _chamber(rib_width=lambda pitch: pitch * 0.0 + 1.0e-4)
    _, helical = _chamber(
        rib_width=lambda pitch: pitch * 0.0 + 1.0e-4,
        helix_angle=np.radians(45.0),
    )

    # theta already carries cos(gamma), so a helical channel is narrower.
    np.testing.assert_allclose(
        helical.channel_local_sector,
        axial.channel_local_sector * np.cos(np.radians(45.0)),
        rtol=1e-12,
    )


# --- rejection --------------------------------------------------------------


def test_a_channel_height_cannot_depend_on_the_section_it_defines() -> None:
    with pytest.raises(DesignInputError, match="tier C"):
        _chamber(channel_height=lambda h: h)


def test_a_helix_angle_cannot_depend_on_the_channel_pitch() -> None:
    with pytest.raises(DesignInputError, match="tier C"):
        _chamber(helix_angle=lambda pitch: 0.0)


# --- shared placements ------------------------------------------------------


def test_a_shared_placement_evaluates_against_each_circuit_in_turn() -> None:
    """A placement may serve several circuits, so it must not bind to one.

    The helix law below reads the channel height, which differs between the two
    circuits. If the placement remembered a single circuit, whichever was
    constructed last would set the helix for both.
    """
    contour = _Cylinder()
    placement = psf.regen.SurfacePlacement(
        n_channel_positions=60,
        helix_angle=lambda h: np.arctan(h / 0.01),
    )

    def circuit(name, height):
        return psf.regen.CoolingCircuit(
            name=name,
            contour=contour,
            cross_section=psf.regen.CrossSectionSquared(rib_width=2.0e-4),
            span=[0.0, 1.0],
            placement=placement,
            channel_height=height,
            walls=[psf.regen.Wall(_Material(), 0.001)],
            coolant_transport=None,
            roughness=0.0,
        )

    shallow = circuit("shallow", 0.002)
    deep = circuit("deep", 0.006)
    psf.regen.ThrustChamber(
        contour=contour,
        cooling_circuits=[shallow, deep],
        combustion_transport=None,
        n_nodes=20,
        compute_gas=False,
    )

    np.testing.assert_allclose(shallow.helix_angle_vals, np.arctan(0.2), rtol=1e-12)
    np.testing.assert_allclose(deep.helix_angle_vals, np.arctan(0.6), rtol=1e-12)


def test_an_area_preserving_helix_holds_the_section_constant() -> None:
    """cos(gamma) * (2*r_base + h) constant preserves area at constant height.

    This is the law the advanced example uses. A cone is needed here so that
    the radius, and therefore the pitch, genuinely varies along the circuit.
    """

    class _Cone(_Cylinder):
        def r(self, x):
            return 0.05 + 0.15 * np.asarray(x, dtype=float)

        def dr_dx(self, x):
            return np.full_like(np.asarray(x, dtype=float), 0.15)

    contour = _Cone()
    height, wall = 0.004, 0.001
    # The section stands off the contour along its normal, so the radius the
    # flow area is set at is the contour radius plus the projected wall stack
    # and half the projected height.
    scale = np.sqrt(1.0 + 0.15**2)
    reference = np.cos(np.radians(20.0)) * 2.0 * (0.05 + (wall + height / 2) * scale)

    def _cone_circuit(rib_width):
        circuit = psf.regen.CoolingCircuit(
            name="cone",
            contour=contour,
            cross_section=psf.regen.CrossSectionSquared(rib_width=rib_width),
            span=[0.0, 1.0],
            placement=psf.regen.SurfacePlacement(
                n_channel_positions=60,
                helix_angle=lambda r, t_wall, dr_dx: np.arccos(
                    np.clip(
                        reference
                        / (2.0 * (r + (t_wall + height / 2) * np.sqrt(1.0 + dr_dx**2))),
                        -1.0,
                        1.0,
                    )
                ),
            ),
            channel_height=height,
            walls=[psf.regen.Wall(_Material(), wall)],
            coolant_transport=None,
            roughness=0.0,
        )
        psf.regen.ThrustChamber(
            contour=contour,
            cooling_circuits=[circuit],
            combustion_transport=None,
            n_nodes=40,
            compute_gas=False,
        )
        return circuit

    circuit = _cone_circuit(2.0e-4)
    open_channel = _cone_circuit(0.0)

    assert np.degrees(circuit.helix_angle_vals[0]) == pytest.approx(20.0)
    # The radius quadruples along the cone, so an unhelixed channel would open
    # out with it. The helix holds the section still instead.
    assert circuit.helix_angle_vals[-1] > circuit.helix_angle_vals[0]
    np.testing.assert_allclose(
        open_channel.A_coolant_vals, open_channel.A_coolant_vals[0], rtol=1e-12
    )
    # A rib of constant material thickness does not remove exactly rib*h from
    # a curved, twisted section, and how much it does remove varies with the
    # helix angle, so the ribbed section drifts by about a part in a thousand.
    np.testing.assert_allclose(
        circuit.A_coolant_vals, circuit.A_coolant_vals[0], rtol=2e-3
    )


def test_construction_preserves_full_rib_context_on_a_resampled_grid() -> None:
    def rib(h, gamma, n_channels_all_circuits, t_wall_total, t_hot_layer):
        assert h == pytest.approx(0.004)
        assert gamma == pytest.approx(0.5)
        assert n_channels_all_circuits == pytest.approx(60)
        assert t_wall_total == pytest.approx(0.001)
        assert t_hot_layer == pytest.approx(0.0)
        return 0.0005

    _, circuit = _chamber(rib_width=rib, helix_angle=0.5)
    before_area = circuit.cross_section.A_coolant(circuit._prof(circuit.reference_path)).copy()
    before_height = circuit.channel_heights.copy()
    x = np.linspace(0.1, 0.9, 7)
    nominal = circuit.reference_path
    centerline = np.column_stack((
        x,
        np.interp(x, nominal[:, 0], nominal[:, 1]),
        np.interp(x, nominal[:, 0], nominal[:, 2]) + 0.3,
    ))
    construction = psf.regen.build_channel_construction(circuit, centerline)
    assert len(construction.planes) == 7
    np.testing.assert_allclose(construction.rib_width, 0.0005)
    np.testing.assert_array_equal(circuit.channel_heights, before_height)
    np.testing.assert_array_equal(
        circuit.cross_section.A_coolant(circuit._prof(circuit.reference_path)), before_area,
    )
