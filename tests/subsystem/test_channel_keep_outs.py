"""Keep-outs: room in the channel layer that surface channels leave clear.

A keep-out narrows every channel in its span, which the thermal model sees,
and the numerical layout bends the channels around it without letting any of
them overlap it, each other, or change order.
"""

import numpy as np
import pytest

import pyskyfire as psf
from pyskyfire.regen.channel_layout import (
    AnalyticalLayoutError,
    resolve_analytical_layout,
    resolve_numerical_layout,
)

TWO_PI = 2.0 * np.pi


class _Cylinder:
    def __init__(self, radius=0.1, length=1.0):
        self.radius = float(radius)
        self.xs = np.array([0.0, length])
        self.x_t = 0.0
        self.r_t = self.r_c = self.r_e = float(radius)

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


def _circuit(contour, name, span, placement):
    return psf.regen.CoolingCircuit(
        name=name,
        contour=contour,
        cross_section=psf.regen.CrossSectionSquared(sidewall_half_thickness=2.0e-4),
        span=span,
        placement=placement,
        channel_height=0.004,
        walls=[psf.regen.Wall(_Material(), 0.001)],
        coolant_transport=None,
        roughness=0.0,
    )


def _chamber(keep_outs=(), *, count=24):
    contour = _Cylinder()
    circuit = _circuit(contour, "only", [0.0, 1.0], psf.regen.SurfacePlacement(count))
    chamber = psf.regen.ThrustChamber(
        contour=contour,
        cooling_circuits=[circuit],
        combustion_transport=None,
        n_nodes=21,
        compute_gas=False,
        keep_outs=list(keep_outs),
    )
    return chamber, circuit


def _keep_out(**overrides):
    options = dict(x_start=0.3, x_end=0.7, width=0.03, count=4, angle=0.0, taper=0.1)
    options.update(overrides)
    return psf.regen.ChannelKeepOut(
        options.pop("x_start"), options.pop("x_end"), options.pop("width"), **options
    )


def test_a_keep_out_takes_its_room_from_every_channel_in_its_span() -> None:
    keep_out = _keep_out()
    chamber, circuit = _chamber([keep_out])
    x = circuit.x_domain
    blocked = np.array([keep_out.blocked_angle(v, chamber.contour) for v in x])

    np.testing.assert_allclose(24 * circuit.channel_angular_width + blocked, TWO_PI)
    inside = (x >= 0.3) & (x <= 0.7)
    np.testing.assert_allclose(
        circuit.channel_angular_width[inside], (TWO_PI - 4 * 0.03 / 0.1) / 24
    )
    outside = (x <= 0.2) | (x >= 0.8)
    np.testing.assert_allclose(circuit.channel_angular_width[outside], TWO_PI / 24)


def test_the_circuit_grid_has_a_station_on_every_keep_out_corner() -> None:
    chamber, circuit = _chamber([_keep_out(x_start=0.33, x_end=0.61, taper=0.07)])

    for corner in (0.26, 0.33, 0.61, 0.68):
        assert np.min(np.abs(circuit.x_domain - corner)) < 1e-12
    assert np.all(np.diff(circuit.x_domain) > 0)


def test_a_keep_out_that_leaves_no_room_is_rejected() -> None:
    with pytest.raises(ValueError, match="no circumference"):
        _chamber([_keep_out(width=0.2)])


def _layout_items(chamber, layout, x):
    """Every channel and keep-out sector at ``x``: (angle, width, label)."""
    items = []
    for channel in layout.expanded_channels():
        circuit = chamber.cooling_circuits[channel.circuit_index]
        path = channel.reference_path
        angle = np.interp(x, path[:, 0], path[:, 2]) % TWO_PI
        width = np.interp(x, circuit.x_domain, circuit.channel_angular_width)
        items.append((angle, width, channel.name))
    for keep_out in chamber.keep_outs:
        width = keep_out.sector_angle(x, chamber.contour)
        if width > 0.0:
            for index in range(keep_out.count):
                angle = (keep_out.angle + index * TWO_PI / keep_out.count) % TWO_PI
                items.append((angle, width, f"keep-out {index}"))
    return sorted(items)


def test_the_numerical_layout_bends_the_channels_around_a_keep_out() -> None:
    chamber, circuit = _chamber([_keep_out()])

    layout = resolve_numerical_layout(chamber)

    assert layout.metadata["g"] == 4
    assert len(layout.channels) == 6
    previous = None
    for x in circuit.x_domain:
        items = _layout_items(chamber, layout, x)
        angles = np.array([item[0] for item in items])
        widths = np.array([item[1] for item in items])
        clearance = (np.roll(angles, -1) - angles) % TWO_PI - 0.5 * (
            widths + np.roll(widths, -1)
        )
        assert clearance.min() > -1e-4 * TWO_PI / 24
        order = [item[2] for item in items if not item[2].startswith("keep-out")]
        if previous is not None:
            start = previous.index(order[0])
            assert previous[start:] + previous[:start] == order
        previous = order


def test_the_analytical_layout_cannot_honour_a_keep_out() -> None:
    chamber, _ = _chamber([_keep_out()])

    with pytest.raises(AnalyticalLayoutError, match="numerical layout"):
        resolve_analytical_layout(chamber)


def test_a_stated_angular_width_holds_where_an_interleaved_pass_ends() -> None:
    contour = _Cylinder()
    placement = psf.regen.SurfacePlacement(12, channel_angular_width=TWO_PI / 24)
    long = _circuit(contour, "long", [0.0, 1.0], placement)
    short = _circuit(contour, "short", [0.2, 1.0], placement)
    psf.regen.ThrustChamber(
        contour=contour,
        cooling_circuits=[long, short],
        combustion_transport=None,
        n_nodes=21,
        compute_gas=False,
    )

    np.testing.assert_allclose(long.channel_angular_width, TWO_PI / 24)


def test_by_default_a_lone_pass_fills_the_room_the_other_leaves() -> None:
    contour = _Cylinder()
    placement = psf.regen.SurfacePlacement(12)
    long = _circuit(contour, "long", [0.0, 1.0], placement)
    short = _circuit(contour, "short", [0.2, 1.0], placement)
    psf.regen.ThrustChamber(
        contour=contour,
        cooling_circuits=[long, short],
        combustion_transport=None,
        n_nodes=21,
        compute_gas=False,
    )

    assert long.channel_angular_width[0] == pytest.approx(TWO_PI / 12)
    assert long.channel_angular_width[-1] == pytest.approx(TWO_PI / 24)


def test_a_keep_out_width_runs_out_over_its_taper() -> None:
    keep_out = _keep_out()

    assert keep_out.width_at(0.5) == pytest.approx(0.03)
    assert keep_out.width_at(0.25) == pytest.approx(0.015)
    assert keep_out.width_at(0.75) == pytest.approx(0.015)
    assert keep_out.width_at(0.1) == 0.0
    np.testing.assert_allclose(keep_out.breakpoints, [0.2, 0.3, 0.7, 0.8])


@pytest.mark.parametrize(
    "overrides",
    [dict(width=0.0), dict(count=0), dict(count=1.5), dict(taper=-0.1),
     dict(x_start=0.5, x_end=0.5)],
)
def test_a_keep_out_rejects_meaningless_dimensions(overrides) -> None:
    with pytest.raises(ValueError):
        _keep_out(**overrides)


def test_a_chamber_pickled_before_keep_outs_blocks_nothing() -> None:
    chamber, _ = _chamber()
    del chamber.keep_outs  # as a chamber unpickled from before they existed

    assert chamber.blocked_angle(0.5) == 0.0
    assert resolve_numerical_layout(chamber).metadata["g"] == 24
