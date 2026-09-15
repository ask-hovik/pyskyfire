"""Render the impeller through the same iframe that the report reader sees."""
from pathlib import Path
import runpy

import pytest

import test_scene_viewer as viewer


def test_impeller_report_iframe_renders(tmp_path, monkeypatch):
    pytest.importorskip("pyturbo")
    pytest.importorskip("plotly")
    script = Path(__file__).resolve().parents[2]/"examples/impeller/impeller_report.py"
    example = runpy.run_path(str(script))
    _, _, report = example["create_report"](tmp_path)
    # The saved report decodes its scene into a blob iframe. Wait for that
    # viewer's own ready signal, so loading the outer report alone cannot pass.
    monkeypatch.setattr(viewer, "READY", """new Promise(resolve => {
        let attempts = 0;
        const timer = setInterval(() => {
            const body = document.querySelector('iframe')?.contentDocument?.body;
            if (body?.dataset.sceneReady || body?.dataset.sceneError || ++attempts > 200) {
                clearInterval(timer);
                resolve({ready: body?.dataset.sceneReady, error: body?.dataset.sceneError});
            }
        }, 50);
    })""")
    result = viewer.run_viewer(tmp_path, report, """(() => {
        const frame = document.querySelector('iframe');
        const canvas = frame.contentDocument.querySelector('canvas');
        return {width: canvas.width, height: canvas.height,
                title: frame.contentDocument.title};
    })()""", tmp_path/"impeller-browser.png")
    assert result["width"] > 0 and result["height"] > 0
    assert "impeller" in result["title"].lower()
    assert (tmp_path/"impeller-browser.png").stat().st_size > 10000
