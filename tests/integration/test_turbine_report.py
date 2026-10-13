"""Exercise the user-facing turbine example with its optional dependencies."""
import os
from pathlib import Path
import subprocess
import sys

import pytest


def test_turbine_example_creates_scene_and_report(tmp_path):
    pytest.importorskip("pyturbo")
    pytest.importorskip("plotly")
    script = Path(__file__).resolve().parents[2]/"examples/turbine/turbine_report.py"
    result = subprocess.run(
        [sys.executable, str(script), "--output-dir", str(tmp_path)],
        cwd=tmp_path, capture_output=True, text=True, timeout=180,
        env={**os.environ, "MPLBACKEND": "Agg", "MPLCONFIGDIR": str(tmp_path/"mpl")},
    )
    assert result.returncode == 0, result.stdout + result.stderr
    for name in ("turbine-scene.html", "turbine-report.html"):
        assert (tmp_path/name).stat().st_size > 1000
    report = (tmp_path/"turbine-report.html").read_text()
    assert "<iframe" in report and "Efficiency (supplied)" in report
    assert "Velocity triangles" in report and "Gauging − metal" in report
