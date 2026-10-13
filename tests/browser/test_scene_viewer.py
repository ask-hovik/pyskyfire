"""Execute the standalone viewer in a locally available Chromium browser."""
import os
import asyncio
import base64
from pathlib import Path
import shutil
from types import SimpleNamespace

import numpy as np
import pytest

from pyskyfire import ducts
from pyskyfire.viz import FieldData, PipeComponent, Scene, SceneNode, Transform

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


def test_heat_exchanger_and_cooled_chamber_share_one_offline_scene(tmp_path):
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        'baffle_scene_example', Path(__file__).resolve().parents[2] / 'examples/heat_exchanger/baffles.py'
    )
    example = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(example)
    chamber, installation, inlets = example.build_baffles(nodes=4)
    result = installation.solve(inlets)
    scene = example.build_scene(chamber, installation, result)
    output = scene.save_html(tmp_path / 'baffle-scene.html')
    observed = run_viewer(tmp_path, output, """(() => {
        const data = JSON.parse(document.getElementById('psf-scene').textContent);
        const selector = document.getElementById('field');
        return {kinds: data.nodes.map(n => n.metadata.kind),
                field: selector.value,
                options: Array.from(selector.options).map(o => o.value),
                legend: document.getElementById('legend').textContent};
    })()""", tmp_path / 'baffle-scene.png')
    assert observed['kinds'] == ['thrust_chamber', 'heat_exchanger']
    assert observed['field'] == 'T_coolant_wall'
    assert {'T_static', 'p_static'} <= set(observed['options'])
    assert 'T_hot_wall' not in observed['options']
    assert 'Coolant wall temperature [K]' in observed['legend']


class NamedTransport:
    """Constant properties, carrying a fluid name the way a real transport does."""

    def __init__(self, fluid):
        self.fluid = fluid

    def get_rho(self, T, p):
        return 70.0

    def get_mu(self, T, p):
        return 1e-5

    def get_cp(self, T, p):
        return 1e4


class Assembly:
    """Named circuit of repeated channels inside a translucent liner, as a chamber is built."""

    def build(self):
        tube = PipeComponent(SimpleNamespace(centerline=np.array([[0., 0, 0], [1., 0, 0]]), radius=0.01),
                             n_theta=8).build().mesh
        liner = PipeComponent(SimpleNamespace(centerline=np.array([[0., 0.3, 0], [1., 0.3, 0]]), radius=0.5),
                              n_theta=24).build().mesh
        circuit = SceneNode("circuit_0", metadata={"kind": "circuit", "label": "Fuel", "fluid": "Hydrogen"})
        # Channels 0.2 apart in y: a ray from the framed view meets one at most.
        circuit.children = [SceneNode(f"channel_{i}", mesh=tube, transform=Transform((0, 0.2*i, 0)),
                                      fields={"p_static": FieldData(np.full(len(tube.points), 5e6), "Pa",
                                                                    "Static pressure")},
                                      metadata={"instance_index": i}) for i in range(4)]
        return SceneNode("assembly", children=[circuit, SceneNode("liner", mesh=liner, opacity=0.3)])


def test_groups_selection_and_parts_under_the_cursor(tmp_path):
    from pyskyfire.thermo import PHASE_UNKNOWN, FluidState
    scene = Scene(color_by="p_static")
    scene.add(Assembly(), name="assembly")
    results = {}
    for name, fluid, p, z in (("feed", "HEOS::Hydrogen[1.000]", 3e6, 0.9), ("ox", "Oxygen", 1e6, -0.9)):
        stations, connections = ducts.generate_path_stations([
            ducts.PathNode((0., 0., z)), ducts.PathNode((1., 0., z), bend=ducts.RoundedBend(radius=0.2)),
            ducts.PathNode((1., 0.6, z))])
        duct = ducts.Pipe(stations, connections, diameter=0.1)
        state = FluidState(T=100., p=p, h=np.nan, quality=np.nan, rho=70., phase=PHASE_UNKNOWN, T_sat=np.nan)
        results[name] = duct.solve(inlet_state=state, mass_flow=5.0, transport=NamedTransport(fluid))
        scene.add(PipeComponent(duct, results[name]), name=name, tags="oxidiser side" if name == "ox" else ())
    output = scene.save_html(tmp_path/"scene.html")
    expression = """(() => {
        const canvas = document.getElementById('c');
        const send = (type, x, y, extra) => canvas.dispatchEvent(new PointerEvent(type, {clientX: x, clientY: y, ...extra}));
        const click = (x, y, extra = {}) => ['pointerdown', 'pointerup'].forEach(t => send(t, x, y, {button: 0, ...extra}));
        const right = (x, y) => ['pointerdown', 'pointerup'].forEach(t => send(t, x, y, {button: 2}));
        const selected = () => (document.body.dataset.selected || '').split('\\n').filter(Boolean);
        const pick = document.getElementById('pick');
        const rows = () => [...pick.querySelectorAll('.item')].map(e => e.textContent);
        // Find where each top-level component can be clicked.
        const at = {};
        for (let y = 30; y < 740; y += 15)
            for (let x = 30; x < 990; x += 15) {
                click(x, y);
                const [id] = selected();
                if (id && !at[id]) at[id] = [x, y];
            }
        click(2, 2);
        const out = {found: Object.keys(at).sort(),
                     groups: document.getElementById('groups').innerText.split('\\n').join(' ')};
        click(...at.feed);
        click(...at.assembly, {shiftKey: true});
        out.added = selected();
        out.listing = document.getElementById('info').innerText;
        click(...at.assembly, {shiftKey: true});
        out.toggled = selected();
        // The first spot a channel is under the cursor lists the whole assembly.
        let spot = null;
        for (let y = 30; y < 740 && !spot; y += 10)
            for (let x = 30; x < 990 && !spot; x += 10) {
                right(x, y);
                if (!pick.hidden && rows().some(r => r.startsWith('channel_'))) spot = [x, y];
            }
        out.rows = rows();
        const circuit = [...pick.querySelectorAll('.item')].find(e => e.textContent.startsWith('circuit_0'));
        circuit.dispatchEvent(new MouseEvent('click', {bubbles: true, shiftKey: true}));
        out.picked = selected();
        out.pickHidden = pick.hidden;
        const form = document.getElementById('create');
        form.querySelector('input').value = 'Rig';
        form.requestSubmit();
        click(2, 2);
        out.cleared = selected();
        const rig = [...document.querySelectorAll('#groups .item')].find(e => e.textContent.includes('Rig'));
        out.rig = rig.textContent;
        rig.dispatchEvent(new MouseEvent('click', {bubbles: true}));
        out.reselected = selected();
        right(...spot);
        const liner = [...pick.querySelectorAll('.item')].find(e => e.textContent === 'liner');
        liner.dispatchEvent(new MouseEvent('click', {bubbles: true}));
        out.sub = selected();
        out.subInfo = document.getElementById('info').innerText;
        const range = document.getElementById('range');
        out.ranges = [...range.options].map(o => o.textContent);
        const limits = () => [...document.querySelectorAll('#limits span')].map(s => parseFloat(s.textContent));
        out.scene = limits();
        range.value = 'Fluid:Hydrogen';
        range.dispatchEvent(new Event('change'));
        out.hydrogen = limits();
        out.caption = document.getElementById('legend').innerText;
        out.stored = JSON.parse(localStorage.getItem('pyskyfire-groups:' + document.title));
        return out;
    })()"""
    value = run_viewer(tmp_path, output, expression)
    assert value["found"] == ["assembly", "feed", "ox"], value["found"]
    # The heat exchanger style HEOS spelling and the plain name are one fluid.
    assert "Fluid Hydrogen (2)" in value["groups"] and "Fluid Oxygen (1)" in value["groups"]
    assert "Tag oxidiser side (1)" in value["groups"]
    assert value["added"] == ["feed", "assembly"] and "2 parts selected" in value["listing"]
    assert value["toggled"] == ["feed"]
    # Named parts are all listed in tree order; of the repeated channels, only
    # the one under the cursor.
    rows = value["rows"]
    assert rows[:2] == ["assembly", "circuit_0 (Fuel) Hydrogen"] and rows[3] == "liner", rows
    assert rows[2].startswith("channel_") and sum(r.startswith("channel_") for r in rows) == 1
    assert value["picked"] == ["feed", "assembly/circuit_0"] and value["pickHidden"]
    assert value["cleared"] == []
    assert value["rig"].startswith("Saved") and "Rig" in value["rig"] and "(2)" in value["rig"]
    assert value["reselected"] == ["feed", "assembly/circuit_0"]
    assert value["stored"] == [{"name": "Rig", "members": ["feed", "assembly/circuit_0"]}]
    assert value["sub"] == ["assembly/liner"] and "Part of assembly" in value["subInfo"]
    assert "Range: Fluid · Hydrogen" in value["ranges"] and "Range: Saved · Rig" in value["ranges"]
    # The scene-wide range spans the oxygen pipe; hydrogen's only its own parts.
    low = float(results["ox"].p_static.min())
    assert value["scene"][0] == pytest.approx(low, rel=1e-4)
    assert value["hydrogen"] == pytest.approx([float(results["feed"].p_static.min()), 5e6], rel=1e-4)
    assert "Static pressure [Pa] · Hydrogen" in value["caption"]
