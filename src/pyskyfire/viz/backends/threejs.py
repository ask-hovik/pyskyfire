"""Offline HTML rendering with mesh references and original scalar values.

The scene is serialized directly instead of exporting a VTK rendering window
through glTF, so fields and node identities survive the browser boundary.
"""
import html
import json
from importlib.resources import files

import numpy as np
import pyvista as pv
from scipy.spatial.transform import Rotation

from .common import field_style, field_styles


def scene_payload(scene):
    # Repeated parts share mesh and field arrays by identity: 280 cooling
    # channels built from three seeds carry three sets of values, not 280.
    meshes, mesh_ids = [], {}
    values, value_ids = [], {}

    def values_ref(f):
        key = id(f.values)
        if key not in value_ids:
            value_ids[key] = len(values)
            values.append([float(x) if np.isfinite(x) else None for x in f.values])
        return value_ids[key]

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
            fields[name] = {"values": values_ref(f), "association": f.association,
                            "units": f.units, "label": f.label}
        return {"id": node.id, "name": node.name, "mesh": mesh_id, "fields": fields,
                "transform": node.transform.matrix.T.ravel().tolist(),
                "anchors": {k: a.transform.matrix.T.ravel().tolist() for k, a in node.anchors.items()},
                "color": pv.Color(node.color).float_rgb, "opacity": node.opacity,
                "visible": node.visible, "metadata": node.metadata, "tags": sorted(node.tags),
                "children": [node_data(child) for child in node.children]}

    style = field_style(scene)
    nodes = []
    for node in scene.nodes:
        data = node_data(node)
        # Only top-level nodes are selectable, and their local transform is the
        # one Scene.add() placed, so these are literally the arguments that
        # reproduce the placement. Extrinsic xyz matches Transform's convention.
        matrix = node.transform.matrix
        data["placement"] = {
            "position": matrix[:3, 3].tolist(),
            "rotation": Rotation.from_matrix(matrix[:3, :3]).as_euler("xyz").tolist(),
        }
        nodes.append(data)
    lut = pv.LookupTable(cmap="viridis", n_values=256)
    return {"version": 2, "title": scene.title, "meshes": meshes, "values": values, "nodes": nodes,
            "style": style, "styles": field_styles(scene),
            "show_axes": bool(getattr(scene, "show_axes", False)),
            "color_table": np.asarray(lut.values)[:, :3].tolist()}


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
#hud,#legend,#info,#groups,#pick{{position:absolute;padding:8px;background:#ffffffe0;border-radius:5px;font:13px system-ui;z-index:1}}
[hidden]{{display:none!important}}
#hud{{top:10px;left:10px;max-width:calc(100% - 280px)}} #legend{{bottom:12px;left:12px;min-width:240px}}
#gradient{{height:14px;margin-top:5px}} #limits{{display:flex;justify-content:space-between}}
#field,#range{{display:block;width:100%;margin-bottom:5px;font:13px system-ui}}
#info{{z-index:2;max-width:360px;box-shadow:0 2px 8px #00000030;padding:0}}
#head{{padding:6px 28px 6px 8px;cursor:move;user-select:none;border-bottom:1px solid #0000001a}}
#title{{font-weight:600}} #body{{padding:6px 8px 8px;max-height:60vh;overflow:auto}}
#info .mono{{font-family:ui-monospace,monospace;font-size:12px;white-space:nowrap}}
#close{{position:absolute;top:3px;right:6px;border:0;background:none;font:16px system-ui;cursor:pointer;color:#555}}
#create{{display:flex;gap:4px;margin-top:8px}} #create input{{flex:1;min-width:0;font:13px system-ui}}
#create button{{font:13px system-ui}}
#groups{{top:10px;right:10px;max-width:260px;max-height:45vh;overflow:auto}}
#pick{{z-index:3;padding:4px 0;box-shadow:0 2px 10px #00000040;max-height:70vh;overflow:auto;background:white}}
.heading{{font-weight:600;padding:2px 6px 4px}} #groups .heading{{padding:0 0 4px}}
.item{{padding:2px 6px;border-radius:3px;cursor:pointer;white-space:nowrap;display:flex;align-items:center;gap:2px}}
.item:hover{{background:#e8eef8}} .note{{color:#777;font-size:12px}}
.drop{{margin-left:auto;border:0;background:none;cursor:pointer;color:#777;font:14px system-ui;padding:0 4px}}
</style></head><body>
<div id="hud">Loading scene…</div><div id="legend" hidden></div><div id="groups" hidden></div>
<div id="info" hidden></div><div id="pick" hidden></div><canvas id="c"></canvas>
<script id="psf-scene" type="application/json">{payload}</script>
<script>{script("three.min.js")}</script>
<script>{script("TrackballControls.js")}</script>
<script>{script("scene_viewer.js")}</script>
</body></html>'''
