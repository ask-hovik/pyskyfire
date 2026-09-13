"""Analytical and numerical cooling-channel layout contracts."""

from types import SimpleNamespace

import numpy as np
import pytest

from pyskyfire.regen.channel_layout import (
    AnalyticalLayoutError,
    resolve_analytical_layout,
    resolve_numerical_layout,
)


def _circuit(
    name,
    count,
    span,
    *,
    angular_rate=0.0,
    initial_angle=0.0,
    occludes=True,
):
    axial = np.linspace(span[0], span[1], 7)
    angle = initial_angle + angular_rate * axial
    centerline = np.column_stack((axial, np.full_like(axial, 0.1), angle))

    class Placement:
        n_channel_positions = count

        def __init__(self):
            self.occludes = occludes

        def dtheta_dx(self, x, contour, circuit=None):
            return angular_rate

    return SimpleNamespace(
        name=name,
        placement=Placement(),
        reference_path=centerline,
        channel_angular_width=np.full_like(axial, 2.0 * np.pi / count),
    )


def _chamber(*circuits):
    return SimpleNamespace(
        cooling_circuits=list(circuits),
        contour=SimpleNamespace(),
    )


def test_analytical_layout_interleaves_equal_counts_without_deforming_paths() -> None:
    first = _circuit("first", 8, (0.0, 1.0), angular_rate=0.2)
    second = _circuit("second", 8, (0.25, 0.8), angular_rate=0.2)

    layout = resolve_analytical_layout(_chamber(first, second))

    assert layout.method == "analytical"
    assert len(layout.channels) == 2
    assert layout.channels[0].angular_offset_rad == 0.0
    assert layout.channels[1].angular_offset_rad == pytest.approx(np.pi / 8.0)
    np.testing.assert_allclose(
        layout.channels[0].reference_path,
        first.reference_path,
    )
    np.testing.assert_allclose(
        layout.channels[1].reference_path[:, :2],
        second.reference_path[:, :2],
    )
    np.testing.assert_allclose(
        layout.channels[1].reference_path[:, 2],
        second.reference_path[:, 2] + np.pi / 8.0,
    )
    assert len(layout.expanded_channels()) == 16


def test_analytical_layout_compensates_a_nominal_constant_phase() -> None:
    first = _circuit("first", 4, (0.0, 1.0), angular_rate=0.1)
    second = _circuit(
        "second",
        4,
        (0.0, 1.0),
        angular_rate=0.1,
        initial_angle=0.2,
    )

    layout = resolve_analytical_layout(_chamber(first, second))

    np.testing.assert_allclose(
        layout.channels[1].reference_path[:, 2] - layout.channels[0].reference_path[:, 2],
        np.pi / 4.0,
    )


def test_analytical_layout_rejects_unequal_interleaving_counts() -> None:
    chamber = _chamber(
        _circuit("first", 8, (0.0, 1.0)),
        _circuit("second", 10, (0.2, 0.8)),
    )

    with pytest.raises(AnalyticalLayoutError, match="equal channel counts"):
        resolve_analytical_layout(chamber)


def test_analytical_layout_rejects_three_different_spans() -> None:
    chamber = _chamber(
        _circuit("first", 8, (0.0, 1.0)),
        _circuit("second", 8, (0.0, 0.8)),
        _circuit("third", 8, (0.1, 1.0)),
    )

    with pytest.raises(AnalyticalLayoutError, match="identical axial spans"):
        resolve_analytical_layout(chamber)


def test_analytical_layout_rejects_mismatched_helix_evolution() -> None:
    chamber = _chamber(
        _circuit("first", 8, (0.0, 1.0), angular_rate=0.1),
        _circuit("second", 8, (0.0, 1.0), angular_rate=0.2),
    )

    with pytest.raises(AnalyticalLayoutError, match="angular/helix evolution"):
        resolve_analytical_layout(chamber)


def test_disjoint_circuits_do_not_need_equal_counts() -> None:
    layout = resolve_analytical_layout(
        _chamber(
            _circuit("first", 8, (0.0, 0.4)),
            _circuit("second", 10, (0.5, 1.0)),
        )
    )

    assert [channel.pattern_count for channel in layout.channels] == [8, 10]
    assert [channel.angular_offset_rad for channel in layout.channels] == [0.0, 0.0]


def test_numerical_layout_returns_patternable_unique_wedge() -> None:
    chamber = _chamber(
        _circuit("first", 6, (0.0, 1.0)),
        _circuit("second", 4, (0.2, 0.8)),
    )

    layout = resolve_numerical_layout(chamber, max_iters=2)

    assert layout.method == "numerical"
    assert len(layout.channels) == 5
    assert all(channel.pattern_count == 2 for channel in layout.channels)
    assert all(
        channel.pattern_angle_rad == pytest.approx(np.pi) for channel in layout.channels
    )
    assert len(layout.expanded_channels()) == 10


def test_resolved_channel_planes_are_lazy_and_cached() -> None:
    channel = resolve_analytical_layout(
        _chamber(_circuit("only", 6, (0.0, 1.0)))
    ).channels[0]

    assert "section_planes" not in channel.__dict__
    first = channel.section_planes

    assert channel.section_planes is first
    assert "section_planes" in channel.__dict__
