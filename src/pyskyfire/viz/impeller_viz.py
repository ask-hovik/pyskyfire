"""PyVista visualisation helpers for pyskyfire pump impellers.

The functions here deliberately consume an ``Impeller`` instance rather than
recomputing geometry.  This mirrors the thrust-chamber visualisation style: the
engineering object prepares geometry; the viz module only turns it into meshes
and plots.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pyvista as pv

try:  # normal package import
    from pyskyfire.pump.impeller_new import Impeller, cylindrical_to_cartesian
except Exception:  # useful when developing this file next to impeller.py
    from impeller_new import Impeller, cylindrical_to_cartesian  # type: ignore


def _structured_surface_from_xyz(xyz: np.ndarray) -> pv.StructuredGrid:
    """Create a PyVista structured surface from ``(ni, nj, 3)`` points."""
    pts = np.asarray(xyz, dtype=float)
    if pts.ndim != 3 or pts.shape[2] != 3:
        raise ValueError("xyz must have shape (ni, nj, 3)")
    ni, nj, _ = pts.shape
    grid = pv.StructuredGrid()
    grid.dimensions = (nj, ni, 1)
    grid.points = pts.reshape(-1, 3)
    return grid


def make_blade_meshes(impeller: Impeller) -> list[pv.StructuredGrid]:
    """Return one camber-surface mesh per blade."""
    return [
        _structured_surface_from_xyz(cylindrical_to_cartesian(blade))
        for blade in impeller.geometry.blade_surfaces
    ]


def make_shroud_meshes(
    impeller: Impeller,
    *,
    n_theta: int = 160,
) -> tuple[pv.StructuredGrid, pv.StructuredGrid]:
    """Create hub and shroud surfaces by revolving their meridional curves."""
    theta = np.linspace(0.0, 2.0 * np.pi, int(n_theta), endpoint=True)

    def revolve(curve_2d: np.ndarray) -> pv.StructuredGrid:
        x = curve_2d[:, 0]
        r = curve_2d[:, 1]
        X, TH = np.meshgrid(x, theta, indexing="ij")
        R, _ = np.meshgrid(r, theta, indexing="ij")
        xyz = np.stack([X, R * np.cos(TH), R * np.sin(TH)], axis=-1)
        return _structured_surface_from_xyz(xyz)

    return revolve(impeller.geometry.hub_curve), revolve(impeller.geometry.shroud_curve)


def make_impeller_meshes(
    impeller: Impeller,
    *,
    n_theta: int = 160,
) -> dict[str, object]:
    """Build all currently available impeller visualisation meshes."""
    hub, shroud = make_shroud_meshes(impeller, n_theta=n_theta)
    blades = make_blade_meshes(impeller)
    return {"hub": hub, "shroud": shroud, "blades": blades}


def plot_impeller(
    impeller: Impeller,
    *,
    show_edges: bool = True,
    show_shrouds: bool = True,
    show_blades: bool = True,
    blade_opacity: float = 0.75,
    shroud_opacity: float = 0.22,
    window_size: tuple[int, int] = (1300, 900),
    screenshot: str | Path | None = None,
    return_plotter: bool = False,
) -> pv.Plotter | None:
    """Plot a preliminary 3D impeller using PyVista.

    Parameters
    ----------
    impeller:
        ``pyskyfire.pump.impeller.Impeller`` instance.
    screenshot:
        Optional path.  If supplied, an off-screen screenshot is written there.
    return_plotter:
        If true, return the ``pv.Plotter`` after adding meshes.  Useful for GUI
        embedding or adding custom annotations.
    """
    off_screen = screenshot is not None
    pl = pv.Plotter(window_size=window_size, off_screen=off_screen)
    pl.set_background("white")

    meshes = make_impeller_meshes(impeller)

    if show_shrouds:
        pl.add_mesh(
            meshes["hub"],
            color="lightsteelblue",
            opacity=shroud_opacity,
            show_edges=show_edges,
            edge_color="gray",
        )
        pl.add_mesh(
            meshes["shroud"],
            color="lightsteelblue",
            opacity=shroud_opacity,
            show_edges=show_edges,
            edge_color="gray",
        )

    if show_blades:
        for blade in meshes["blades"]:  # type: ignore[index]
            pl.add_mesh(
                blade,
                color="tomato",
                opacity=blade_opacity,
                show_edges=show_edges,
                edge_color="black",
            )

    d = impeller.dimensions
    title = (
        f"Impeller: Q={impeller.inputs.Q:.4g} m³/s, H={impeller.inputs.H:.4g} m, "
        f"n={impeller.inputs.n:.0f} rpm, d2={d.d2:.4g} m, z={d.blade_count}"
    )
    pl.add_text(title, position="upper_left", font_size=10, color="black")
    pl.add_axes()
    pl.view_isometric()
    pl.camera.zoom(1.2)

    if screenshot is not None:
        pl.show(screenshot=str(screenshot), auto_close=not return_plotter)
    elif not return_plotter:
        pl.show()

    return pl if return_plotter else None


def export_impeller_surfaces(
    impeller: Impeller,
    path: str | Path,
    *,
    binary: bool = True,
) -> Path:
    """Export current visualisation surfaces as a single PolyData mesh.

    This is a visual/CAD-reference mesh, not a watertight manufacturing model.
    """
    meshes = make_impeller_meshes(impeller)
    merged = pv.PolyData()
    for mesh in [meshes["hub"], meshes["shroud"], *meshes["blades"]]:  # type: ignore[index]
        merged = merged.merge(mesh.extract_surface().triangulate())
    out = Path(path)
    merged.save(out, binary=binary)
    return out
