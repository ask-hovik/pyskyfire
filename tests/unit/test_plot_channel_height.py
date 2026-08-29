"""Cooling-channel height visualization."""

from types import SimpleNamespace

import numpy as np

from pyskyfire.viz.plot_regen import PlotChannelHeight, PlotChannelWidth


def test_channel_height_plot_uses_raw_precomputed_profile() -> None:
    first = SimpleNamespace(
        name="Half Pass",
        x_domain=np.array([0.0, 0.5, 1.0]),
        channel_heights=np.array([0.001, 0.002, 0.003]),
    )
    second = SimpleNamespace(
        name="Full Pass",
        x_domain=np.array([-0.5, 0.0, 0.5]),
        channel_heights=np.array([0.004, 0.005, 0.006]),
    )
    chamber = SimpleNamespace(cooling_circuits=[first, second])

    plot = PlotChannelHeight(chamber)

    assert len(plot.fig.data) == 2
    assert [trace.name for trace in plot.fig.data] == ["Half Pass", "Full Pass"]
    np.testing.assert_allclose(plot.fig.data[0].x, first.x_domain)
    np.testing.assert_allclose(plot.fig.data[0].y, first.channel_heights)
    np.testing.assert_allclose(plot.fig.data[1].x, second.x_domain)
    np.testing.assert_allclose(plot.fig.data[1].y, second.channel_heights)
    assert plot.fig.layout.yaxis.title.text == "Channel Height, h (m)"


def test_channel_height_plot_can_select_one_circuit() -> None:
    circuits = [
        SimpleNamespace(
            name=f"Circuit {index}",
            x_domain=np.array([0.0, 1.0]),
            channel_heights=np.array([index + 1.0, index + 2.0]),
        )
        for index in range(2)
    ]

    plot = PlotChannelHeight(
        SimpleNamespace(cooling_circuits=circuits),
        circuit_index=1,
    )

    assert len(plot.fig.data) == 1
    assert plot.fig.data[0].name == "Circuit 1"
    np.testing.assert_allclose(plot.fig.data[0].y, [2.0, 3.0])


def test_channel_width_plot_compares_hot_side_chord_with_reference() -> None:
    theta = np.array([0.2, 0.1])
    radius = np.array([0.05, 0.08])
    circuit = SimpleNamespace(
        name="Full Pass",
        x_domain=np.array([0.0, 1.0]),
        centerlines=[
            np.column_stack(
                [np.array([0.0, 1.0]), radius, np.zeros(2)]
            )
        ],
        channel_local_sector=theta,
    )
    reference = SimpleNamespace(
        name="RECOP design points",
        x=np.array([0.0, 0.5]),
        tube_width_od=np.array([0.01, 0.02]),
    )

    plot = PlotChannelWidth(
        SimpleNamespace(cooling_circuits=[circuit]),
        reference_profiles=[reference],
    )

    assert len(plot.fig.data) == 2
    assert plot.fig.data[0].name == "Full Pass (calculated)"
    np.testing.assert_allclose(
        plot.fig.data[0].y,
        2.0 * radius * np.sin(theta / 2.0),
    )
    assert plot.fig.data[1].name == "RECOP design points"
    assert plot.fig.data[1].mode == "markers"
    np.testing.assert_allclose(plot.fig.data[1].y, reference.tube_width_od)
    assert plot.fig.layout.yaxis.title.text == "Tube Width OD (m)"


def test_channel_width_plot_rejects_mismatched_reference_arrays() -> None:
    circuit = SimpleNamespace(
        name="Full Pass",
        x_domain=np.array([0.0]),
        centerlines=[np.array([[0.0, 0.05, 0.0]])],
        channel_local_sector=np.array([0.1]),
    )
    reference = SimpleNamespace(
        x=np.array([0.0, 1.0]),
        tube_width_od=np.array([0.01]),
    )

    with np.testing.assert_raises_regex(ValueError, "arrays must match"):
        PlotChannelWidth(
            SimpleNamespace(cooling_circuits=[circuit]),
            reference_profiles=[reference],
        )
