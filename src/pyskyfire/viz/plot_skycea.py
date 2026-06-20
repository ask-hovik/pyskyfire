# src/pyskyfire/viz/transport_property.py
import numpy as np
import plotly.graph_objects as go
from .core import PlotBase

_PROP_INFO = {
    "M":   ("M", ""),
    "gamma": ("γ", ""),
    "T":   ("T", "K"),
    "p":   ("p", "Pa"),          # note: p_map in Aerothermodynamics stores bar in maps;
                                 # here we still plot the equilibrium column raw unless you prefer Pa.
    "h":   ("h", "kJ/kg"),       # maps are kJ/kg; change label if you convert
    "cp":  ("cₚ (mass)", "kJ/(kg·K)"),
    "cv":  ("cᵥ (mass)", "kJ/(kg·K)"),
    "k":   ("k", "W/(m·K)"),
    "mu":  ("μ", "Pa·s"),
    "Pr":  ("Pr", "–"),
    "rho": ("ρ", "kg/m³"),
    "a":   ("a", "m/s"),
}

class PlotTransportProperty(PlotBase):
    """
    Plot a single transport-property map (equilibrium column vs x) for one or more
    Aerothermodynamics objects.

    Each object must have:
      - .x_nodes (built by compute_aerothermodynamics)
      - .<prop>_map with shape (Nx, Nt); we use column 0 (equilibrium).
    """

    def __init__(self, *ats, prop: str, template: str = "plotly_white"):
        if prop not in _PROP_INFO:
            raise ValueError(f"Unknown property '{prop}'. Valid keys: {list(_PROP_INFO)}")

        super().__init__(go.Figure())
        self.template(template)

        map_attr = f"{prop}_map"
        y_label, unit = _PROP_INFO[prop]

        for i, at in enumerate(ats):
            x = np.asarray(getattr(at, "x_nodes"), dtype=float)
            Z = np.asarray(getattr(at, map_attr), dtype=float)   # (Nx, Nt)
            y = Z[:, 0]  # equilibrium column

            name = getattr(at, "name", f"Set {i+1}")
            self.fig.add_trace(go.Scatter(x=x, y=y, mode="lines", name=name, showlegend=True))

        self.fig.update_layout(
            title=f"{y_label} map",
            xaxis=dict(title="Axial position, x (m)"),
            yaxis=dict(title=f"{y_label}" + (f" ({unit})" if unit else "")),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=50, b=55),
        )


class PlotTransportPropertyField(PlotBase):
    """
    Plot a full precomputed transport-property field as z = property(x, T).

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