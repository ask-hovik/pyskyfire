"""Exercise the user-facing impeller example with its optional dependency."""
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest


def test_impeller_example_creates_scene_report_and_splitter_geometry(tmp_path):
    pytest.importorskip("pyturbo")
    pytest.importorskip("plotly")
    script = Path(__file__).resolve().parents[2]/"examples/impeller/impeller_report.py"
    result = subprocess.run(
        [sys.executable, str(script), "--output-dir", str(tmp_path)],
        cwd=tmp_path, capture_output=True, text=True, timeout=120,
        env={**os.environ, "MPLBACKEND":"Agg", "MPLCONFIGDIR":str(tmp_path/"mpl")},
    )
    assert result.returncode == 0, result.stdout+result.stderr
    for name in ("impeller-scene.html", "impeller-report.html", "impeller-geometry.npz"):
        assert (tmp_path/name).stat().st_size > 1000
    report = (tmp_path/"impeller-report.html").read_text()
    assert "<iframe" in report and "Volumetric efficiency (supplied)" in report
    with np.load(tmp_path/"impeller-geometry.npz", allow_pickle=False) as geometry:
        assert geometry["suction"].shape == (33,129,3)
        assert geometry["splitter_suction"].shape == (33,129,3)
        assert np.isfinite(geometry["suction"]).all()
