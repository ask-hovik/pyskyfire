# src/pyskyfire/viz/transport_property.py
import numpy as np
import plotly.graph_objects as go
from .core import PlotBase

_PROP_INFO = {
    "M":   ("M", ""),
    "gamma": ("γ", ""),
    "T":   ("T", "K"),
    "p":   ("p", "Pa"),
    "h":   ("h", "J/kg"),
    "cp":  ("cₚ (mass)", "J/(kg·K)"),
    "cv":  ("cᵥ (mass)", "J/(kg·K)"),
    "k":   ("k", "W/(m·K)"),
    "mu":  ("μ", "Pa·s"),
    "Pr":  ("Pr", "–"),
    "rho": ("ρ", "kg/m³"),
    "a":   ("a", "m/s"),
    "v":   ("v", "m/s"),
}


def _iter_results(results):
    if results is None:
        return []
    if isinstance(results, dict):
        if any(key in results for key in ("x", "x_wall", "x_heat_flux")):
            return [results]
        return list(results.values())
    if isinstance(results, (list, tuple)):
        return list(results)
    return [results]


def result_gas_x(results) -> np.ndarray:
    """Return the unique gas-side coordinates used by solved results."""
    grids = []
    for result in _iter_results(results):
        def field(name):
            if isinstance(result, dict):
                return result.get(name)
            return getattr(result, name, None)

        result_grids = []
        for attribute in ("x_wall", "x_heat_flux"):
            grid = field(attribute)
            if grid is not None:
                result_grids.append(np.asarray(grid, dtype=float))
        if not result_grids or all(grid.size == 0 for grid in result_grids):
            grid = field("x")
            if grid is not None:
                result_grids.append(np.asarray(grid, dtype=float))
        grids.extend(result_grids)
    if not grids:
        raise ValueError("results do not contain an axial calculation grid")
    return np.unique(np.concatenate(grids))


def _visualization_x(at, *, x, results, nodes):
    supplied = sum(value is not None for value in (x, results, nodes))
    if supplied > 1:
        raise ValueError("provide only one of x, results, or nodes")
    if results is not None:
        return result_gas_x(results)
    if nodes is not None:
        if not isinstance(nodes, (int, np.integer)) or nodes < 2:
            raise ValueError("nodes must be an integer of at least 2")
        contour = getattr(at, "contour", None)
        if contour is None:
            raise ValueError("attach a contour before using nodes")
        return np.linspace(contour.xs[0], contour.xs[-1], int(nodes))
    if x is not None:
        values = np.asarray(x, dtype=float)
        if values.ndim != 1 or values.size == 0:
            raise ValueError("x must be a non-empty one-dimensional array")
        return values

    # Compatibility for precomputed/legacy transport objects. New live CEA
    # objects do not create x_nodes when a contour is attached.
    legacy = getattr(at, "x_nodes", None)
    if legacy is not None:
        return np.asarray(legacy, dtype=float)
    raise ValueError(
        "choose the visualization grid with results=..., nodes=..., or x=..."
    )

class PlotTransportProperty(PlotBase):
    """
    Plot one equilibrium transport property against axial position.

    The current live-solve interface is sampled through ``get_<prop>``. Pass
    ``results`` to reuse the gas-side coordinates of an actual regenerative
    run, ``nodes`` for a uniform visualization grid, or ``x`` for explicit
    coordinates. Repeated property plots reuse the transport object's live
    state cache. Old precomputed-map objects remain supported for compatibility.
    """

    def __init__(
        self,
        *ats,
        prop: str,
        results=None,
        nodes: int | None = None,
        x=None,
        template: str = "plotly_white",
    ):
        if prop not in _PROP_INFO:
            raise ValueError(f"Unknown property '{prop}'. Valid keys: {list(_PROP_INFO)}")

        super().__init__(go.Figure())
        self.template(template)

        map_attr = f"{prop}_map"
        y_label, unit = _PROP_INFO[prop]

        for i, at in enumerate(ats):
            x_values = _visualization_x(
                at, x=x, results=results, nodes=nodes
            )
            getter = getattr(at, f"get_{prop}", None)
            if getter is not None:
                y = np.asarray([getter(x_i) for x_i in x_values], dtype=float)
            elif prop == "v":
                # Legacy maps carry no velocity map; rebuild it from M and a.
                M = np.asarray(getattr(at, "M_map"), dtype=float)[:, 0]
                a = np.asarray(getattr(at, "a_map"), dtype=float)[:, 0]
                y = M * a
            else:
                Z = np.asarray(getattr(at, map_attr), dtype=float)
                y = Z[:, 0]
                if prop == "p":       # legacy maps store bar
                    y = y * 1.0e5
                elif prop in {"h", "cp", "cv"}:  # legacy maps store kJ units
                    y = y * 1.0e3

            name = getattr(at, "name", f"Set {i+1}")
            self.fig.add_trace(go.Scatter(
                x=x_values, y=y, mode="lines", name=name, showlegend=True
            ))

        self.fig.update_layout(
            title=f"{y_label} profile",
            xaxis=dict(title="Axial position, x (m)"),
            yaxis=dict(title=f"{y_label}" + (f" ({unit})" if unit else "")),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=50, b=55),
        )


class PlotTransportPropertyField(PlotBase):
    """
    Plot a full precomputed transport-property field as z = property(x, T).

    Only legacy precomputed-map objects carry these attributes; live CEA
    objects solve states on demand and have no property field to plot.

    Uses:
      - at.x_nodes
      - at.T_grid
      - at.<prop>_map
    """

    def __init__(
        self,
        at,
        *,
        prop: str,
        mode: str = "surface",  # "surface", "heatmap", or "contour"
        template: str = "plotly_white",
    ):
        if prop not in _PROP_INFO:
            raise ValueError(f"Unknown property '{prop}'. Valid keys: {list(_PROP_INFO)}")

        super().__init__(go.Figure())
        self.template(template)

        map_attr = f"{prop}_map"
        label, unit = _PROP_INFO[prop]

        x = np.asarray(at.x_nodes, dtype=float)
        T = np.asarray(at.T_grid, dtype=float)
        Z = np.asarray(getattr(at, map_attr), dtype=float)

        # Make X have the same shape as T and Z: (Nx, Nt)
        X = np.repeat(x[:, None], T.shape[1], axis=1)

        title = f"{label} field"
        z_title = f"{label}" + (f" ({unit})" if unit else "")

        if mode == "surface":
            self.fig.add_trace(go.Surface(
                x=X,
                y=T,
                z=Z,
                colorbar=dict(title=z_title),
            ))

            self.fig.update_layout(
                title=title,
                scene=dict(
                    xaxis_title="Axial position, x (m)",
                    yaxis_title="Temperature, T (K)",
                    zaxis_title=z_title,
                ),
                margin=dict(l=0, r=0, t=50, b=0),
            )

        elif mode == "heatmap":
            self.fig.add_trace(go.Heatmap(
                x=x,
                y=T[0, :],
                z=Z.T,
                colorbar=dict(title=z_title),
            ))

            self.fig.update_layout(
                title=title,
                xaxis=dict(title="Axial position, x (m)"),
                yaxis=dict(title="Temperature, T (K)"),
                margin=dict(l=70, r=20, t=50, b=60),
            )

        elif mode == "contour":
            self.fig.add_trace(go.Contour(
                x=x,
                y=T[0, :],
                z=Z.T,
                colorbar=dict(title=z_title),
                contours=dict(showlabels=True),
            ))

            self.fig.update_layout(
                title=title,
                xaxis=dict(title="Axial position, x (m)"),
                yaxis=dict(title="Temperature, T (K)"),
                margin=dict(l=70, r=20, t=50, b=60),
            )

        else:
            raise ValueError("mode must be 'surface', 'heatmap', or 'contour'")
