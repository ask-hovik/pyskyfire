"""Execute the standalone viewer in a locally available Chromium browser."""
import os
import asyncio
import base64
from pathlib import Path
import shutil

import numpy as np
import pytest

from pyskyfire import ducts
from pyskyfire.viz import FieldData, PipeComponent, Scene

READY = """new Promise(resolve => {
    let attempts = 0;
    const timer = setInterval(() => {
        if (document.body?.dataset.sceneReady || document.body?.dataset.sceneError || ++attempts > 200) {
            clearInterval(timer);
            resolve({ready: document.body?.dataset.sceneReady, error: document.body?.dataset.sceneError});
        }
    }, 50);
})"""


def run_viewer(tmp_path, output, expression, screenshot=None):
    """Load a saved viewer, wait for it to display, then evaluate ``expression``."""
    choreographer = pytest.importorskip("choreographer")
    from choreographer.browsers.chromium import Chromium
    browser = os.environ.get("PYSKYFIRE_TEST_BROWSER")
    if browser is None:
        browser = next((path for name in ("chromium", "chromium-browser", "google-chrome", "brave-browser")
                        if (path := shutil.which(name))), None)
    if browser is None:
        pytest.skip("Install a Chromium browser or set PYSKYFIRE_TEST_BROWSER")
    # Brave's shell wrapper masks exit codes and leaves pipe-owning children
    # behind on timeout. Use its adjacent binary when present.
    executable = Path(browser).resolve()
    if executable.name == "brave-browser" and executable.with_name("brave").is_file():
        browser = str(executable.with_name("brave"))

    class TestChromium(Chromium):
        def get_cli(self):
            return [*super().get_cli(), "--use-gl=angle", "--use-angle=swiftshader",
                    "--enable-unsafe-swiftshader", "--no-first-run", "--no-default-browser-check"]

    async def check():
        async with choreographer.Browser(path=browser, browser_cls=TestChromium,
                                          enable_gpu=True, tmp_dir=str(tmp_path)) as driver:
            tab = await driver.create_tab(output.as_uri(), width=1000, height=750, window=True)

            async def evaluate(source):
                result = await tab.send_command("Runtime.evaluate", params={
                    "expression": source, "awaitPromise": True, "returnByValue": True})
                assert "exceptionDetails" not in result["result"], result["result"]["exceptionDetails"]
                return result["result"]["result"]["value"]

            state = await evaluate(READY)
            assert state.get("ready") == "true", state
            assert not state.get("error"), state
            value = await evaluate(expression)
            if screenshot is not None:
                image = await tab.send_command("Page.captureScreenshot")
                screenshot.write_bytes(base64.b64decode(image["result"]["data"]))
            return value

    return asyncio.run(asyncio.wait_for(check(), timeout=60))


def solved_scene(**kwargs):
    # A solver pipe, so the refined rounded bend reaches the renderer too.
    stations, connections = ducts.generate_path_stations([
        ducts.PathNode((0., 0., 0.)),
        ducts.PathNode((1., 0., 0.), bend=ducts.RoundedBend(radius=0.3)),
        ducts.PathNode((1., 1., 0.)),
    ])
    pipe = ducts.Pipe(stations, connections, diameter=0.2)
    scene = Scene(color_by="p_static", **kwargs)
    scene.add(PipeComponent(pipe, {
        "p_static": FieldData(np.linspace(3e6, 1e6, len(stations)), "Pa", "Static pressure"),
        "T_static": FieldData(np.linspace(120, 300, len(stations)), "K", "Static temperature"),
    }), name="points")
    other = scene.add(PipeComponent(pipe), name="cells", position=(0, 0, 0.4), rotation=(0, 0.5, 0))
    other.fields["p_static"] = FieldData(np.linspace(1e6, 3e6, len(other.mesh.triangles)), "Pa", "Static pressure", "cell")
    return scene


def test_offline_scene_renders_point_and_cell_fields(tmp_path):
    output = solved_scene(show_axes=True).save_html(tmp_path/"scene.html")
    screenshot = tmp_path/"scene.png"
    legend = run_viewer(tmp_path, output, "document.getElementById('legend').textContent", screenshot)
    assert "Static pressure [Pa]" in legend
    assert Path(screenshot).stat().st_size > 1000


def test_click_reports_placement_and_field_selector_recolors(tmp_path):
    output = solved_scene().save_html(tmp_path/"scene.html")
    # Sweep the canvas until a click lands on geometry, so the assertions do not
    # depend on where the framed view happens to put the pipes.
    expression = """(() => {
        const canvas = document.getElementById('c');
        const click = (x, y) => ['pointerdown', 'pointerup'].forEach(type =>
            canvas.dispatchEvent(new PointerEvent(type, {clientX: x, clientY: y, button: 0})));
        let hit = null;
        for (let y = 60; y < 700 && !hit; y += 20)
            for (let x = 60; x < 950 && !hit; x += 20) {
                click(x, y);
                if (document.body.dataset.selected) hit = {x: x, y: y};
            }
        const info = document.getElementById('info');
        const box = info.getBoundingClientRect();
        const before = document.body.dataset.selected;
        const text = info.textContent;
        const selector = document.getElementById('field');
        const options = [...selector.options].map(o => o.textContent);
        selector.value = 'T_static';
        selector.dispatchEvent(new Event('change'));
        const legend = document.getElementById('legend').textContent;
        click(2, 2);
        return {hit: hit, selected: before, text: text, options: options, legend: legend,
                near: hit && Math.abs(box.left - hit.x) < 200 && Math.abs(box.top - hit.y) < 200,
                cleared: !document.body.dataset.selected, hidden: info.hidden};
    })()"""
    value = run_viewer(tmp_path, output, expression)
    assert value["hit"], "no click landed on the pipes"
    # A click anywhere on a component selects the node Scene.add() created.
    assert value["selected"] in ("points", "cells"), value["selected"]
    assert "position=" in value["text"] and "rotation=" in value["text"]
    if value["selected"] == "cells":
        assert "(0, 0, 0.4)" in value["text"] and "(0, 0.5, 0)" in value["text"]
    else:
        assert "(0, 0, 0)" in value["text"]
    assert "inlet" in value["text"] and "outlet" in value["text"]
    assert value["near"], "the info box did not appear beside the click"
    assert value["options"] == ["None (component colors)", "Static pressure", "Static temperature"]
    assert "Static temperature [K]" in value["legend"]
    assert value["cleared"] and value["hidden"], "clicking empty space did not deselect"
