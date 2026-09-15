"""Convenience plotting and surface export through the shared impeller scene."""
from pathlib import Path

import numpy as np
import pyvista as pv

from .components import ImpellerComponent
from .scene import Scene

__all__ = ["make_impeller_meshes", "make_blade_meshes", "make_shroud_meshes",
           "plot_impeller", "export_impeller_surfaces"]


def make_impeller_meshes(impeller, *, n_theta=160):
    root = ImpellerComponent(impeller, n_theta=n_theta).build()
    meshes = {"blades": [], "splitters": []}
    for node in root.children:
        faces = np.column_stack((np.full(len(node.mesh.triangles), 3), node.mesh.triangles))
        mesh = pv.PolyData(node.transform.apply(node.mesh.points), faces)
        if node.name.startswith("blade_"):
            meshes["blades"].append(mesh)
        elif node.name.startswith("splitter_"):
            meshes["splitters"].append(mesh)
        else:
            meshes[node.name] = mesh
    return meshes


def make_blade_meshes(impeller):
    return make_impeller_meshes(impeller)["blades"]


def make_shroud_meshes(impeller, *, n_theta=160):
    meshes = make_impeller_meshes(impeller, n_theta=n_theta)
    return meshes["hub"], meshes["shroud"]


def plot_impeller(impeller, *, show_edges=False, show_shrouds=True, show_blades=True,
                  blade_opacity=1., shroud_opacity=0.22, window_size=(1300,900),
                  screenshot=None, return_plotter=False):
    scene = Scene(title="Impeller")
    root = scene.add(ImpellerComponent(impeller))
    for node in root.children:
        endwall = node.name in ("hub", "shroud")
        node.visible = show_shrouds if endwall else show_blades
        node.opacity = shroud_opacity if endwall else blade_opacity
    pl = scene.to_plotter(off_screen=screenshot is not None)
    pl.window_size = window_size
    for actor in pl.actors.values():
        if hasattr(actor, "prop"):
            actor.prop.show_edges = show_edges
    if screenshot is not None:
        pl.show(screenshot=str(screenshot), auto_close=not return_plotter)
    elif not return_plotter:
        pl.show()
    return pl if return_plotter else None


def export_impeller_surfaces(impeller, path, *, binary=True):
    """Export the flow-path endwalls and closed blade meshes as a reference mesh."""
    meshes = make_impeller_meshes(impeller)
    merged = pv.PolyData()
    for mesh in [meshes["hub"], meshes["shroud"], *meshes["blades"], *meshes["splitters"]]:
        merged = merged.merge(mesh)
    out = Path(path)
    merged.save(out, binary=binary)
    return out
