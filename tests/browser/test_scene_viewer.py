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


def test_offline_scene_renders_point_and_cell_fields(tmp_path):
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
    # A solver pipe, so the refined rounded bend reaches the renderer too.
    stations, connections = ducts.generate_pipe_stations([
        ducts.PathNode((0., 0., 0.), 0.2),
        ducts.PathNode((1., 0., 0.), 0.2, bend=ducts.RoundedBend(radius=0.3)),
        ducts.PathNode((1., 1., 0.), 0.2),
    ])
    pipe = ducts.Pipe(stations, connections)
    scene = Scene(color_by="p_static")
    scene.add(PipeComponent(pipe, {"p_static": FieldData(np.linspace(3e6, 1e6, len(stations)), "Pa", "Static pressure")}), name="points")
    other = scene.add(PipeComponent(pipe), name="cells", position=(0, 0, 0.4))
    other.fields["p_static"] = FieldData(np.linspace(1e6, 3e6, len(other.mesh.triangles)), "Pa", "Static pressure", "cell")
    output = scene.save_html(tmp_path/"scene.html")
    screenshot = tmp_path/"scene.png"
    class TestChromium(Chromium):
        def get_cli(self):
            return [*super().get_cli(), "--use-gl=angle", "--use-angle=swiftshader",
                    "--enable-unsafe-swiftshader", "--no-first-run", "--no-default-browser-check"]

    async def check():
        async with choreographer.Browser(path=browser, browser_cls=TestChromium,
                                          enable_gpu=True, tmp_dir=str(tmp_path)) as driver:
            tab = await driver.create_tab(output.as_uri(), width=1000, height=750, window=True)
            result = await tab.send_command("Runtime.evaluate", params={
                "expression": """new Promise(resolve => {
                    let attempts = 0;
                    const timer = setInterval(() => {
                        if (document.body?.dataset.sceneReady || document.body?.dataset.sceneError || ++attempts > 200) {
                            clearInterval(timer);
                            resolve({ready: document.body?.dataset.sceneReady,
                                     error: document.body?.dataset.sceneError,
                                     legend: document.getElementById('legend')?.textContent});
                        }
                    }, 50);
                })""", "awaitPromise": True, "returnByValue": True,
            })
            value = result["result"]["result"]["value"]
            assert value.get("ready") == "true", value
            assert not value.get("error"), value
            assert "Static pressure [Pa]" in value["legend"]
            image = await tab.send_command("Page.captureScreenshot")
            screenshot.write_bytes(base64.b64decode(image["result"]["data"]))

    asyncio.run(asyncio.wait_for(check(), timeout=40))
    assert Path(screenshot).stat().st_size > 1000
