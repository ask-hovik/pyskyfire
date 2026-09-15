from types import SimpleNamespace
from unittest.mock import Mock
import inspect

import numpy as np
import pytest

from pyskyfire.regen.cross_section import CrossSectionSquared
from pyskyfire.viz.engine_viz import (
    Engine3DViewer,
    build_sections_for_centerlines,
    make_engine_3d,
)


def test_engine_3d_viewer_close_releases_plotter():
    plotter = Mock()
    viewer = Engine3DViewer(plotter=plotter, data_url="data:text/html;base64,")

    viewer.close()
    viewer.close()

    plotter.close.assert_called_once_with()
    plotter.deep_clean.assert_called_once_with()
    assert viewer.plotter is None

    with pytest.raises(RuntimeError, match="closed Engine3DViewer"):
        viewer.show()


def test_engine_3d_viewer_requires_a_color_for_each_circuit():
    chamber = SimpleNamespace(cooling_circuits=[object(), object()])

    with pytest.raises(ValueError, match="color for every cooling circuit"):
        make_engine_3d(chamber, show=False, circuit_colors=["white"])


def test_engine_visualizer_defaults_to_numerical_layout():
    assert inspect.signature(make_engine_3d).parameters["layout"].default == "numerical"


def test_squared_visualization_draws_the_wall_stack_and_the_rib():
    x = np.linspace(0.0, 0.02, 5)
    circuit = SimpleNamespace(
        cross_section=CrossSectionSquared(sidewall_half_thickness=lambda pitch_mid: (0.25 * pitch_mid) / 2),
        reference_path=np.column_stack((x, np.full_like(x, 0.1), np.zeros_like(x))),
        x_domain=x,
        channel_heights=np.full_like(x, 0.004),
        channel_local_sector=np.full_like(x, 0.2),
        helix_angle_vals=np.zeros_like(x),
        channel_wall_thickness=lambda v: np.full_like(v, 0.001),
        hot_side_layer_thickness=lambda v: np.full_like(v, 0.002),
        total_thickness=lambda v: np.full_like(v, 0.003),
    )
    chamber = SimpleNamespace(cooling_circuits=[circuit])

    section = build_sections_for_centerlines([circuit.reference_path], [0], chamber)[0][0]

    # The passage is drawn on top of the whole wall stack, liner included, not
    # sunk into it and not offset by the channel wall alone.
    floor, roof = section[:8], section[8:]
    np.testing.assert_allclose(np.linalg.norm(floor[:, 1:], axis=1), 0.103)
    np.testing.assert_allclose(np.linalg.norm(roof[:, 1:], axis=1), 0.107)

    # A rib taking a quarter of the mid-height pitch leaves the rest open.
    rib = 0.25 * (0.103 + 0.002) * 0.2
    opening = np.ptp(np.arctan2(floor[:, 2], floor[:, 1])) * 0.103
    assert opening == pytest.approx(0.103 * 0.2 - rib, rel=1e-4)
