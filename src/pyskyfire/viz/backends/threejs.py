"""Offline HTML rendering with mesh references and original scalar values.

The scene is serialized directly instead of exporting a VTK rendering window
through glTF, so fields and node identities survive the browser boundary.
"""
import html
import json
from importlib.resources import files

import numpy as np
import pyvista as pv

from .common import field_style


def scene_payload(scene):
    meshes, mesh_ids = [], {}

    def node_data(node):
        mesh_id = None
        if node.mesh is not None:
            key = id(node.mesh)
            if key not in mesh_ids:
                mesh_ids[key] = len(meshes)
                meshes.append({"points": node.mesh.points.ravel().tolist(),
                               "triangles": node.mesh.triangles.ravel().tolist()})
            mesh_id = mesh_ids[key]
        fields = {}
        for name, f in node.fields.items():
            fields[name] = {"values": [float(x) if np.isfinite(x) else None for x in f.values],
                            "association": f.association, "units": f.units, "label": f.label}
        return {"id": node.id, "name": node.name, "mesh": mesh_id, "fields": fields,
                "transform": node.transform.matrix.T.ravel().tolist(),
                "anchors": {k: a.transform.matrix.T.ravel().tolist() for k, a in node.anchors.items()},
                "color": pv.Color(node.color).float_rgb, "opacity": node.opacity,
                "visible": node.visible, "metadata": node.metadata,
                "children": [node_data(child) for child in node.children]}

    style = field_style(scene)
    nodes = [node_data(node) for node in scene.nodes]
    lut = pv.LookupTable(cmap="viridis", n_values=256)
    return {"version": 1, "title": scene.title, "meshes": meshes, "nodes": nodes,
            "style": style, "color_table": np.asarray(lut.values)[:, :3].tolist()}


def to_html(scene):
    assets = files("pyskyfire.viz").joinpath("assets")
    payload = json.dumps(scene_payload(scene), allow_nan=False, separators=(",", ":")).replace("<", "\\u003c")
    def script(name):
        return assets.joinpath(name).read_text(encoding="utf-8").replace("</script", "<\\/script")
    title = html.escape(scene.title)
    return f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title><style>
html,body{{height:100%;margin:0;background:white;overflow:hidden}}
canvas{{display:block;width:100%;height:100%}}
#hud,#legend{{position:absolute;padding:8px;background:#ffffffe0;border-radius:5px;font:13px system-ui;z-index:1}}
#hud{{top:10px;left:10px}} #legend{{bottom:12px;left:12px;min-width:240px}}
#gradient{{height:14px;margin-top:5px}} #limits{{display:flex;justify-content:space-between}}
</style></head><body>
<div id="hud">Loading scene…</div><div id="legend" hidden></div><canvas id="c"></canvas>
<script id="psf-scene" type="application/json">{payload}</script>
<script>{script("three.min.js")}</script>
<script>{script("OrbitControls.js")}</script>
<script>{script("scene_viewer.js")}</script>
</body></html>'''
