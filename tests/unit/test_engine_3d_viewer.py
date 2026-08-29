from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from pyskyfire.viz.engine_viz import Engine3DViewer, make_engine_3d


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
