"""Deprecated engine viewer compatibility API and legacy mesh helpers."""
import warnings
from collections.abc import Sequence
from pathlib import Path
import base64
from dataclasses import dataclass

import numpy as np
import pyvista as pv

from pyskyfire.regen.channel_construction import build_channel_construction
from pyskyfire.regen.cad_geometry import transverse_construction
from pyskyfire.regen.channel_layout import (
    ChannelLayout,
    circular_pattern as circular_pattern,
    collect_channels as collect_channels,
    gcd_list as gcd_list,
    redistribute_channels as redistribute_channels,
    resolve_analytical_layout,
    resolve_numerical_layout,
    tile_centerlines as tile_centerlines,
    wedge_info as wedge_info,
)

# ---------------------------------------------------------------------------
# Returned dataclass
# ---------------------------------------------------------------------------


@dataclass
class Engine3DViewer:
    """Deprecated compatibility handle; new code should use :class:`Scene`."""

    plotter: pv.Plotter | None
    data_url: str

    def __post_init__(self):
        warnings.warn("Engine3DViewer is deprecated; use pyskyfire.viz.Scene", DeprecationWarning, stacklevel=2)

    def save_html(self, path: str | Path) -> Path:
        """Write the interactive viewer to a standalone HTML file."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)

        prefix = "data:text/html;base64,"
        if not self.data_url.startswith(prefix):
            raise ValueError("Expected a base64-encoded HTML data URL")

        encoded_html = self.data_url.removeprefix(prefix)
        html = base64.b64decode(encoded_html).decode("utf-8")

        path.write_text(html, encoding="utf-8")
        return path

    def show(self):
        """Open the PyVista desktop viewer."""
        if self.plotter is None:
            raise RuntimeError("Cannot show a closed Engine3DViewer")
        return self.plotter.show()

    def close(self) -> None:
        """Release PyVista rendering resources."""
        plotter = self.plotter
        if plotter is None:
            return

        # Release our reference before cleanup so a repeated call is harmless
        # and the plotter is destroyed while PyVista's module globals still
        # exist, rather than during interpreter shutdown.
        self.plotter = None
        try:
            plotter.close()
        finally:
            plotter.deep_clean()


# ---------------------------------------------------------------------------
# Geometry helpers for cross-sections
# ---------------------------------------------------------------------------


def build_sections_for_centerlines(
    centerlines,
    cid_full,
    thrust_chamber,
    *,
    stride: int = 1,
    n_arc_pts: int = 8,
):
    """
    Build 3D cross-section polygons for each centerline, using the actual
    CoolingCircuit + cross_section definitions.

    Parameters
    ----------
    centerlines :
        List of (Ni,3) centerlines in cylindrical [x, r, theta].
    cid_full :
        List/array of circuit indices for each centerline.
    thrust_chamber :
        Thrust chamber providing cooling circuit geometry.
    stride : int
        Subsampling stride along centerlines (increase to reduce mesh size).
    n_arc_pts : int
        Number of angular samples on inner and outer arcs.

    Returns
    -------
    sections_per_channel : list[list[np.ndarray]]
        For each channel, a list of (Nk,3) section polygons in XYZ.
    """
    if stride < 1:
        raise ValueError("stride must be at least 1")
    if n_arc_pts < 2:
        raise ValueError("n_arc_pts must be at least 2")
    circuits = thrust_chamber.cooling_circuits
    sections_per_channel: list[list[np.ndarray]] = []
    for centerline, circuit_index in zip(centerlines, cid_full, strict=True):
        construction = transverse_construction(
            build_channel_construction(circuits[int(circuit_index)], centerline)
        )
        count = len(construction.planes)
        indices = list(range(0, count, stride))
        if indices[-1] != count - 1:
            indices.append(count - 1)
        sections_per_channel.append([
            construction.polygon(index, n_arc_pts) for index in indices
        ])
    return sections_per_channel



# ---------------------------------------------------------------------------
# PyVista lofting
# ---------------------------------------------------------------------------


def loft_sections_to_polydata(sections, *, cap_ends: bool = False) -> pv.PolyData:
    """
    Loft a list of cross-section polygons into a surface.

    Parameters
    ----------
    sections : list of (Ni, 3) arrays
        All sections must have the same number of points N,
        ordered consistently around the perimeter.
    cap_ends : bool
        If True, add a polygon at each end to close the channel.

    Returns
    -------
    pv.PolyData
        Lofted surface mesh.
    """
    sections = [np.asarray(sec, float) for sec in sections]
    if len(sections) < 2:
        raise ValueError("Need at least 2 sections to loft")

    N = sections[0].shape[0]
    for sec in sections:
        if sec.shape != (N, 3):
            raise ValueError("All sections must have shape (N,3) with the same N")

    M = len(sections)
    points = np.vstack(sections)

    def vid(i_sec, j_pt):
        return i_sec * N + j_pt

    faces = []
    for i_sec in range(M - 1):
        for j in range(N):
            j_next = (j + 1) % N
            v0 = vid(i_sec, j)
            v1 = vid(i_sec, j_next)
            v2 = vid(i_sec + 1, j_next)
            v3 = vid(i_sec + 1, j)
            faces.extend([4, v0, v1, v2, v3])

    if cap_ends:
        start_face = [N]
        start_face.extend(vid(0, j) for j in range(N))
        faces.extend(start_face)

        end_face = [N]
        end_face.extend(vid(M - 1, j) for j in reversed(range(N)))
        faces.extend(end_face)

    faces = np.asarray(faces, dtype=np.int64)
    return pv.PolyData(points, faces)


# ---------------------------------------------------------------------------
# PyVista entry point
# ---------------------------------------------------------------------------


def make_engine_3d(
    thrust_chamber,
    *,
    stride: int = 2,
    n_arc_pts: int = 8,
    show: bool = True,
    circuit_colors: Sequence[str] | None = None,
    layout: str | ChannelLayout = "numerical",
):
    """
    Deprecated: use ``Scene.add(ThrustChamberComponent(...))``.

    This compatibility wrapper delegates to the shared scene implementation,
    returns the historical viewer handle, and honours ``show``. New scene
    callers should explicitly call ``scene.show()`` to open a window.

    Parameters
    ----------
    thrust_chamber
        Thrust chamber object with cooling circuit definitions.
    stride : int
        Subsampling stride along each centerline (increase to reduce mesh size).
    n_arc_pts : int
        Number of angular samples used for inner/outer arcs.
    show : bool
        If True, open the interactive PyVista window.
    circuit_colors : sequence of str, optional
        Colors indexed by cooling-circuit order. When omitted, the default
        categorical palette is used.
    layout : {"numerical", "analytical"} or ChannelLayout, optional
        Channel layout to draw. Numerical redistribution remains the default.
        A previously resolved layout may also be supplied directly.

    Returns
    -------
    Engine3DViewer
        Interactive viewer and its serialized HTML data URL.
    """

    warnings.warn(
        "make_engine_3d is deprecated; use Scene.add(ThrustChamberComponent(...))",
        DeprecationWarning, stacklevel=2,
    )
    from .scene import Scene
    from .components import ThrustChamberComponent
    scene = Scene()
    scene.add(ThrustChamberComponent(thrust_chamber, stride=stride, n_arc_pts=n_arc_pts,
                                     circuit_colors=circuit_colors, layout=layout))
    plotter = scene.to_plotter(off_screen=not show)
    try:
        viewer = Engine3DViewer(plotter=plotter, data_url=scene.data_url)
        if show:
            plotter.show(auto_close=False)
        return viewer
    except Exception:
        plotter.close()
        plotter.deep_clean()
        raise
