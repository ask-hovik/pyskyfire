# pyskyfire/viz/mesh3d.py
from __future__ import annotations
from typing import Optional, Tuple, Iterable
import numpy as np
import plotly.graph_objects as go
#import trimesh

from .core import PlotBase

def _require_trimesh():
    try:
        import trimesh  # type: ignore
        return trimesh
    except Exception as e:
        raise RuntimeError(
            "Optional dependency 'trimesh' is required for EmbedSTL. "
            "Install with: pip install trimesh"
        ) from e

class EmbedSTL(PlotBase):
    """
    Thin wrapper that turns an STL on disk into a Plotly Mesh3d figure.
    Keeps the figure open for further Plotly commands via PlotBase.
    """

    def __init__(self,
                 stl_path: str,
                 *,
                 color: str = "#9cc4ff",
                 opacity: float = 1.0,
                 show_wireframe: bool = False,
                 template: str = "plotly_white",
                 process: bool = False):
        """
        Parameters
        ----------
        stl_path : str
            Path to the STL file on disk.
        color : str
            Mesh color (any Plotly color string).
        opacity : float
            Mesh opacity in [0, 1].
        show_wireframe : bool
            If True, also overlays an approximate wireframe using triangle edges.
        template : str
            Plotly layout template.
        process : bool
            Let trimesh 'process' the mesh (repairs, merges). Defaults off to preserve input.
        """
        super().__init__(go.Figure())
        tm = _require_trimesh()

        m = tm.load(stl_path, force="mesh", process=process)
        if isinstance(m, tm.Scene):
            # Merge scene geometry into a single mesh
            m = tm.util.concatenate(tuple(g for g in m.geometry.values()))
        if not isinstance(m, tm.Trimesh):
            raise RuntimeError("Loaded STL is not a single triangle mesh.")

        self._vertices = np.asarray(m.vertices)
        self._faces = np.asarray(m.faces, dtype=int)

        self.template(template)
        self._add_mesh_trace(color=color, opacity=opacity)

        if show_wireframe:
            self.add_wireframe()

        # sensible 3D defaults
        self.layout(
            scene=dict(
                aspectmode="data",
                xaxis_title="x",
                yaxis_title="y",
                zaxis_title="z",
            ),
            title="STL view"
        )

    # -----------------------------
    # public API / fluent helpers
    # -----------------------------
    @property
    def vertices(self) -> np.ndarray:
        return self._vertices

    @property
    def faces(self) -> np.ndarray:
        return self._faces

    def add_wireframe(self,
                      *,
                      name: str = "Wireframe",
                      width: float = 1.0) -> "EmbedSTL":
        """Overlay a wireframe by drawing all triangle edges as line segments."""
        V, F = self._vertices, self._faces
        edges = np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]], axis=0)
        seg = V[edges.reshape(-1)]
        self.fig.add_trace(go.Scatter3d(
            x=seg[:, 0], y=seg[:, 1], z=seg[:, 2],
            mode="lines",
            line=dict(width=width),
            name=name,
            showlegend=False
        ))
        return self

    def add_section_plane(self,
                          plane_point: Tuple[float, float, float],
                          plane_normal: Tuple[float, float, float],
                          *,
                          name: str = "Section",
                          width: float = 3.0) -> "EmbedSTL":
        """
        Slice the mesh with a plane and overlay polyline(s) of the intersection.
        """
        tm = _require_trimesh()
        m = tm.Trimesh(vertices=self._vertices, faces=self._faces, process=False)
        section = m.section(plane_origin=np.array(plane_point, float),
                            plane_normal=np.array(plane_normal, float))
        if section is None:
            return self
        # discretize into polylines
        try:
            paths = section.discrete
        except Exception:
            paths = [np.asarray(section.vertices)]
        for pts in paths:
            if pts is None or len(pts) == 0:
                continue
            self.fig.add_trace(go.Scatter3d(
                x=pts[:, 0], y=pts[:, 1], z=pts[:, 2],
                mode="lines",
                name=name,
                line=dict(width=width)
            ))
        return self

    def recenter(self) -> "EmbedSTL":
        """Translate so the model centroid is near the origin."""
        c = self._vertices.mean(axis=0)
        self.traces(x=lambda x: x - c[0], y=lambda y: y - c[1], z=lambda z: z - c[2])
        return self

    # -----------------------------
    # internal
    # -----------------------------
    def _add_mesh_trace(self, *, color: str, opacity: float):
        V, F = self._vertices, self._faces
        self.fig.add_trace(go.Mesh3d(
            x=V[:, 0], y=V[:, 1], z=V[:, 2],
            i=F[:, 0], j=F[:, 1], k=F[:, 2],
            color=color,
            opacity=opacity,
            lighting=dict(ambient=0.45, diffuse=0.9, specular=0.2),
            name="STL"
        ))
