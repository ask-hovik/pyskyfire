"""Desktop rendering of a scene, retaining named scalar arrays and node IDs."""
import numpy as np
import pyvista as pv

from .common import field_style, is_visible


def make_plotter(scene, *, off_screen=False):
    style = field_style(scene)
    plotter = pv.Plotter(off_screen=off_screen)
    try:
        plotter.set_background("white")
        for node in scene.walk():
            if node.mesh is None or not is_visible(node):
                continue
            source = node.mesh
            faces = np.column_stack((np.full(len(source.triangles), 3), source.triangles))
            mesh = pv.PolyData(node.world_transform.apply(source.points), faces.ravel())
            for name, data in node.fields.items():
                target = mesh.point_data if data.association == "point" else mesh.cell_data
                target[name] = data.values.copy()
            mesh.field_data["scene_node_id"] = [node.id]
            kwargs = {"color": node.color}
            if style and style["name"] in node.fields:
                data = node.fields[style["name"]]
                kwargs = dict(scalars=style["name"], preference=data.association,
                              cmap="viridis", clim=style["limits"], nan_color="#808080",
                              scalar_bar_args={"title": f'{style["label"]} [{style["units"]}]'})
            plotter.add_mesh(mesh, name=node.id, opacity=node.opacity,
                             smooth_shading=True, show_edges=False, **kwargs)
        plotter.hide_axes()
        plotter.view_isometric()
        plotter.reset_camera()
        return plotter
    except Exception:
        plotter.close()
        plotter.deep_clean()
        raise
