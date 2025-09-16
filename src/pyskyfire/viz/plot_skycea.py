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
