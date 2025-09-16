from __future__ import annotations
from typing import Sequence, Optional, Any, Iterable, Union
import numpy as np
import plotly.graph_objects as go
from .core import PlotBase


class PlotContour(PlotBase):
    """
    Plot bell-nozzle or toroidal-aerospike contours, mirrored about the x-axis.

    Each contour can be either:
      • classic: has .xs and .rs
      • aerospike: has .xs_outer, .rs_outer, .xs_inner, .rs_inner

    If a contour has a .name attribute, it's used for the legend label.
    """

    def __init__(
        self,
        *contours: Any,
        show_labels: bool = True,
        title: str = "Contour Profiles",
        template: str = "plotly_white",
    ):
        super().__init__(go.Figure())
        self.template(template)

        # ------------ colorway helper (same color for +r and -r; rotate per contour) ------------
        try:
            tpl = pio.templates[template] if template in pio.templates else pio.templates["plotly"]
            colorway = list(tpl.layout.colorway) if tpl.layout.colorway else []
        except Exception:
            colorway = []

        if not colorway:
            # plotly's default colorway fallback
            colorway = ["#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A",
                        "#19D3F3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52"]

        color_idx = 0
        def next_color():
            nonlocal color_idx
            c = colorway[color_idx % len(colorway)]
            color_idx += 1
            return c

        for idx, contour in enumerate(contours):
            label = getattr(contour, "name", f"Contour {idx + 1}")
            group = f"contour-{idx+1}"
            color = next_color()

            # ---------- Classic bell nozzle: xs / rs ----------
            if hasattr(contour, "xs") and hasattr(contour, "rs"):
                xs = np.asarray(contour.xs)
                rs = np.asarray(contour.rs)

                # +r (legend)
                self.fig.add_trace(go.Scatter(
                    x=xs, y=rs,
                    mode="lines",
                    name=label,
                    legendgroup=group,
                    showlegend=show_labels,
                    line=dict(color=color),
                ))
                # -r (no legend, same color)
                self.fig.add_trace(go.Scatter(
                    x=xs, y=-rs,
                    mode="lines",
                    name=f"{label} (mirror)",
                    legendgroup=group,
                    showlegend=False,
                    line=dict(color=color),
                ))

            # ---------- Aerospike: outer + inner walls ----------
            elif (hasattr(contour, "xs_outer") and hasattr(contour, "rs_outer")
                  and hasattr(contour, "xs_inner") and hasattr(contour, "rs_inner")):

                xs_o = np.asarray(contour.xs_outer)
                rs_o = np.asarray(contour.rs_outer)
                xs_i = np.asarray(contour.xs_inner)
                rs_i = np.asarray(contour.rs_inner)

                # Outer wall — in legend (+r)
                self.fig.add_trace(go.Scatter(
                    x=xs_o, y=rs_o,
                    mode="lines",
                    name=label,
                    legendgroup=group,
                    showlegend=show_labels,
                    line=dict(color=color),
                ))
                # Outer wall mirror — same color, no legend
                self.fig.add_trace(go.Scatter(
                    x=xs_o, y=-rs_o,
                    mode="lines",
                    name=f"{label} (mirror)",
                    legendgroup=group,
                    showlegend=False,
                    line=dict(color=color),
                ))
                # Inner wall — same color but dashed, no legend (keep legend clean)
                self.fig.add_trace(go.Scatter(
                    x=xs_i, y=rs_i,
                    mode="lines",
                    name=f"{label} (inner)",
                    legendgroup=group,
                    showlegend=False,
                    line=dict(color=color, dash="dot"),
                ))
                # Inner wall mirror — same color/dash, no legend
                self.fig.add_trace(go.Scatter(
                    x=xs_i, y=-rs_i,
                    mode="lines",
                    name=f"{label} (inner mirror)",
                    legendgroup=group,
                    showlegend=False,
                    line=dict(color=color, dash="dot"),
                ))

            else:
                raise TypeError(f"Object at index {idx} lacks the required contour attributes.")

        # Axes, aspect, labels
        self.layout(
            title=title or None,
            xaxis=dict(title="Axial position, x (m)"),
            yaxis=dict(title="Radius, r (m)"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )
        # Equal aspect: lock y to x
        self.fig.update_yaxes(scaleanchor="x", scaleratio=1)

class PlotWallTemperature(PlotBase):
    """
    Build a wall-temperature plot from one or more cooling-data dicts.

    Each input dict must contain:
      - 'x': 1D array of axial positions
      - 'T': 2D array (len(x) × n_layers)

    Optional per-dataset label:
      - 'name': str  (used as legendgroup; defaults to 'Set {i+1}')
    """
    def __init__(self,
                 *cooling_data_dicts,
                 plot_hot: bool = True,
                 plot_interfaces: bool = False,
                 plot_coolant_wall: bool = False,
                 template: str = "plotly_white"):
        super().__init__(go.Figure())
        self.template(template)

        # --- domain building ---
        for i, data in enumerate(cooling_data_dicts):
            x = np.asarray(data["x"])
            T = np.asarray(data["T"])
            dataset_name = data.get("name", f"Set {i+1}")

            cols = []
            if plot_coolant_wall and T.shape[1] > 1:
                cols.append(1)
            if plot_interfaces and T.shape[1] > 2:
                cols += list(range(2, T.shape[1] - 1))
            if plot_hot:
                cols.append(T.shape[1] - 1)

            for col in cols:
                dash = "dash" if (plot_interfaces and 2 <= col <= T.shape[1] - 2) else "solid"
                self.fig.add_trace(go.Scatter(
                    x=x,
                    y=T[:, col],
                    name=(
                        "Cold Wall" if (col == 1 and plot_coolant_wall)
                        else "Hot Wall" if (col == T.shape[1] - 1 and plot_hot)
                        else f"Interface {col - 1}"
                    ),
                    legendgroup=dataset_name,
                    showlegend=True,
                    line=dict(dash=dash),
                ))

        # sensible defaults
        self.layout(
            title="Wall Temperature",
            xaxis=dict(title="Axial Position (m)"),
            yaxis=dict(title="Temperature (K)")
        )


class PlotCoolantTemperature(PlotBase):
    """
    Expects dicts with keys: 'x', 'T_static', optional 'name'.
    """
    def __init__(self, *cooling_data_dicts):
        super().__init__(go.Figure())
        self.template("plotly_white")

        for i, data in enumerate(cooling_data_dicts):
            x = np.asarray(data["x"])
            y = np.asarray(data["T_static"])
            name = data.get("name", f"Set {i+1}")
            self.fig.add_trace(go.Scatter(x=x, y=y, mode="lines", name=name, showlegend=True))

        self.fig.update_layout(
            title="Coolant Temperature",
            xaxis=dict(title="Axial Position (m)"),
            yaxis=dict(title="Coolant Temperature (K)"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )


class PlotCoolantPressure(PlotBase):
    """
    Expects dicts with keys: 'x', 'p_static', 'p_stagnation' (either/both may exist), optional 'name'.
    Use flags to include static and/or stagnation traces.
    """
    def __init__(self, *cooling_data_dicts, static: bool = True, stagnation: bool = True):
        super().__init__(go.Figure())
        self.template("plotly_white")

        for i, data in enumerate(cooling_data_dicts):
            x = np.asarray(data["x"])
            name = data.get("name", f"Set {i+1}")

            if static and "p_static" in data:
                y_static = np.asarray(data["p_static"])
                self.fig.add_trace(go.Scatter(
                    x=x, y=y_static, mode="lines",
                    name=f"{name} — Static", showlegend=True
                ))

            if stagnation and "p_stagnation" in data:
                y_stag = np.asarray(data["p_stagnation"])
                self.fig.add_trace(go.Scatter(
                    x=x, y=y_stag, mode="lines",
                    name=f"{name} — Stagnation", showlegend=True
                ))

        self.fig.update_layout(
            title="Coolant Pressure",
            xaxis=dict(title="Axial Position (m)"),
            yaxis=dict(title="Pressure (Pa)"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )


class PlotHeatFlux(PlotBase):
    """
    Expects dicts with keys: 'x', 'dQ_dA', optional 'name'.
    """
    def __init__(self, *cooling_data_dicts):
        super().__init__(go.Figure())
        self.template("plotly_white")

        for i, data in enumerate(cooling_data_dicts):
            x = np.asarray(data["x"])
            y = np.asarray(data["dQ_dA"])
            name = data.get("name", f"Set {i+1}")
            self.fig.add_trace(go.Scatter(x=x, y=y, mode="lines", name=name, showlegend=True))

        self.fig.update_layout(
            title="Heat Flux",
            xaxis=dict(title="Axial Position (m)"),
            yaxis=dict(title="Heat Flux (W/m²)"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )


class PlotVelocity(PlotBase):
    """
    Expects dicts with keys: 'x', 'velocity', optional 'name'.
    """
    def __init__(self, *cooling_data_dicts):
        super().__init__(go.Figure())
        self.template("plotly_white")

        for i, data in enumerate(cooling_data_dicts):
            x = np.asarray(data["x"])
            y = np.asarray(data["velocity"])
            name = data.get("name", f"Set {i+1}")
            self.fig.add_trace(go.Scatter(x=x, y=y, mode="lines", name=name, showlegend=True))

        self.fig.update_layout(
            title="Coolant Velocity",
            xaxis=dict(title="Axial Position (m)"),
            yaxis=dict(title="Velocity (m/s)"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )


IndexLike = Optional[Union[int, Iterable[int]]]

def _normalize_indices(thrust_chamber, circuit_index: IndexLike):
    circuits = thrust_chamber.cooling_circuit_group.circuits
    n = len(circuits)
    if circuit_index is None:
        return list(range(n))
    if isinstance(circuit_index, int):
        return [circuit_index]
    # assume iterable
    idxs = list(circuit_index)
    return idxs

class PlotdAdxThermalHotGas(PlotBase):
    def __init__(self, thrust_chamber, circuit_index: IndexLike = None):
        super().__init__(go.Figure())
        self.template("plotly_white")

        circuits = thrust_chamber.cooling_circuit_group.circuits
        indices = _normalize_indices(thrust_chamber, circuit_index)

        for idx in indices:
            circuit = circuits[idx]
            x = np.asarray(circuit.x_domain)
            y = np.asarray(circuit.dA_dx_thermal_exhaust_vals)
            name = getattr(circuit, "name", f"Circuit {idx}")
            self.fig.add_trace(go.Scatter(x=x, y=y, mode="lines", name=name, showlegend=True))

        self.fig.update_layout(
            title="dA/dx Thermal Hot Gas",
            xaxis=dict(title="Axial Position, x (m)"),
            yaxis=dict(title="dA/dx"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )

class PlotdAdxThermalCoolant(PlotBase):
    def __init__(self, thrust_chamber, circuit_index: IndexLike = None):
        super().__init__(go.Figure())
        self.template("plotly_white")

        circuits = thrust_chamber.cooling_circuit_group.circuits
        indices = _normalize_indices(thrust_chamber, circuit_index)

        for idx in indices:
            circuit = circuits[idx]
            x = np.asarray(circuit.x_domain)
            y = np.asarray(circuit.dA_dx_thermal_coolant_vals)
            name = getattr(circuit, "name", f"Circuit {idx}")
            self.fig.add_trace(go.Scatter(x=x, y=y, mode="lines", name=name, showlegend=True))

        self.fig.update_layout(
            title="dA/dx Thermal Coolant",
            xaxis=dict(title="Axial Position, x (m)"),
            yaxis=dict(title="dA/dx"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )

class PlotCoolantArea(PlotBase):
    def __init__(self, thrust_chamber, circuit_index: IndexLike = None):
        super().__init__(go.Figure())
        self.template("plotly_white")

        circuits = thrust_chamber.cooling_circuit_group.circuits
        indices = _normalize_indices(thrust_chamber, circuit_index)

        for idx in indices:
            circuit = circuits[idx]
            x = np.asarray(circuit.x_domain)
            y = np.asarray(circuit.A_coolant_vals)
            name = getattr(circuit, "name", f"Circuit {idx}")
            self.fig.add_trace(go.Scatter(x=x, y=y, mode="lines", name=name, showlegend=True))

        self.fig.update_layout(
            title="Coolant Cross-Sectional Area",
            xaxis=dict(title="Axial Position, x (m)"),
            yaxis=dict(title="A (m²)"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )

class PlotdAdxCoolantArea(PlotBase):
    def __init__(self, thrust_chamber, circuit_index: IndexLike = None):
        super().__init__(go.Figure())
        self.template("plotly_white")

        circuits = thrust_chamber.cooling_circuit_group.circuits
        indices = _normalize_indices(thrust_chamber, circuit_index)

        for idx in indices:
            circuit = circuits[idx]
            x = np.asarray(circuit.x_domain)
            y = np.asarray(circuit.dA_dx_coolant_vals)
            name = getattr(circuit, "name", f"Circuit {idx}")
            self.fig.add_trace(go.Scatter(x=x, y=y, mode="lines", name=name, showlegend=True))

        self.fig.update_layout(
            title="dA/dx Coolant Area",
            xaxis=dict(title="Axial Position, x (m)"),
            yaxis=dict(title="dA/dx"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )

class PlotHydraulicDiameter(PlotBase):
    def __init__(self, thrust_chamber, circuit_index: IndexLike = None):
        super().__init__(go.Figure())
        self.template("plotly_white")

        circuits = thrust_chamber.cooling_circuit_group.circuits
        indices = _normalize_indices(thrust_chamber, circuit_index)

        for idx in indices:
            circuit = circuits[idx]
            x = np.asarray(circuit.x_domain)
            y = np.asarray(circuit.Dh_coolant_vals)
            name = getattr(circuit, "name", f"Circuit {idx}")
            self.fig.add_trace(go.Scatter(x=x, y=y, mode="lines", name=name, showlegend=True))

        self.fig.update_layout(
            title="Coolant Hydraulic Diameter",
            xaxis=dict(title="Axial Position, x (m)"),
            yaxis=dict(title="Dₕ (m)"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )

class PlotRadiusOfCurvature(PlotBase):
    def __init__(self, thrust_chamber, circuit_index: IndexLike = None):
        super().__init__(go.Figure())
        self.template("plotly_white")

        circuits = thrust_chamber.cooling_circuit_group.circuits
        indices = _normalize_indices(thrust_chamber, circuit_index)

        for idx in indices:
            circuit = circuits[idx]
            x = np.asarray(circuit.x_domain)
            y = np.asarray(circuit.radius_of_curvature_vals)
            name = getattr(circuit, "name", f"Circuit {idx}")
            self.fig.add_trace(go.Scatter(x=x, y=y, mode="lines", name=name, showlegend=True))

        self.fig.update_layout(
            title="Radius of Curvature",
            xaxis=dict(title="Axial Position, x (m)"),
            yaxis=dict(title="R (m)"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )

# src/pyskyfire/viz/temperature_profile.py
import numpy as np
import plotly.graph_objects as go
from .core import PlotBase

class PlotTemperatureProfile(PlotBase):
    """
    Temperature profile T(y) across gas boundary layer, wall, and coolant
    boundary layer at a given axial location x_query.

    Inputs:
      - results: dict with keys "x", "T" (N×n_layers), "T_static", "dQ_dA", "p_static"
      - thrust_chamber: provides combustion_transport and cooling circuits
      - circuit_index: which cooling circuit to use for coolant props
      - x_query: axial location (m) to sample
    """

    def __init__(self, results, thrust_chamber, circuit_index: int, x_query: float, n_bl: int = 1000):
        super().__init__(go.Figure())
        self.template("plotly_white")

        # ---- 1) nearest axial node ----
        x_arr = np.asarray(results["x"], dtype=float)
        i = int(np.argmin(np.abs(x_arr - x_query)))
        x0 = float(x_arr[i])

        # ---- 2) extract temps & heat flux ----
        T_inf = float(thrust_chamber.combustion_transport.get_T(x0))
        T_hw  = float(results["T"][i, -1])    # hot-wall
        T_cw  = float(results["T"][i,  1])    # coolant-side wall
        T_c   = float(results["T_static"][i]) # bulk coolant
        qpp   = float(results["dQ_dA"][i])

        # ---- 3) gas-side BL (1/7th power law) ----
        k_g   = float(thrust_chamber.combustion_transport.get_k(x0))
        h_g   = qpp / max(T_inf - T_hw, 1e-12)
        delta_g = 7.0 * k_g / max(h_g, 1e-30)
        y_g = np.linspace(-delta_g, 0.0, n_bl)
        T_g = np.where(
            np.abs(y_g) <= delta_g,
            T_inf + (T_hw - T_inf) * (1.0 - (np.abs(y_g) / max(delta_g, 1e-30)) ** (1.0 / 7.0)),
            T_inf,
        )

        # ---- 4) wall layers (hot→coolant) ----
        Ts_rev   = results["T"][i, 1:]       # coolant-side → hot-side
        Ts_wall  = Ts_rev[::-1]              # hot-side → coolant-side (interfaces included)
        walls    = thrust_chamber.wall_group.walls
        thicknesses = [float(w.thickness(x0)) for w in walls]
        y_w = np.insert(np.cumsum(thicknesses), 0, 0.0) if thicknesses else np.array([0.0])
        wall_thickness = float(y_w[-1])

        # ---- 5) coolant BL (1/7th power law) ----
        p_static = float(results["p_static"][i])
        T_film   = 0.5 * (T_c + T_cw)
        coolant  = thrust_chamber.cooling_circuit_group.circuits[circuit_index].coolant_transport
        k_c      = float(coolant.get_k(T_film, p_static))
        h_c      = qpp / max(T_cw - T_c, 1e-12)
        delta_c  = 7.0 * k_c / max(h_c, 1e-30)
        y_c = np.linspace(0.0, delta_c, n_bl)
        T_cBL = np.where(
            y_c <= delta_c,
            T_c + (T_cw - T_c) * (1.0 - (y_c / max(delta_c, 1e-30)) ** (1.0 / 7.0)),
            T_c,
        )

        # ---- combine domains ----
        y_all = np.concatenate([y_g, y_w, y_c + wall_thickness])
        T_all = np.concatenate([T_g, Ts_wall, T_cBL])

        # Extents for freestreams
        x_min = -delta_g - 2.0 * wall_thickness
        x_max = wall_thickness + delta_c + 2.0 * wall_thickness

        # ---- shaded regions (under the curves) ----
        self.fig.add_shape(
            type="rect", xref="x", yref="paper",
            x0=-delta_g, x1=0.0, y0=0.0, y1=1.0,
            fillcolor="lightgray", opacity=0.7, line_width=0, layer="below"
        )
        self.fig.add_shape(
            type="rect", xref="x", yref="paper",
            x0=0.0, x1=wall_thickness, y0=0.0, y1=1.0,
            fillcolor="gray", opacity=0.8, line_width=0, layer="below"
        )
        self.fig.add_shape(
            type="rect", xref="x", yref="paper",
            x0=wall_thickness, x1=wall_thickness + delta_c, y0=0.0, y1=1.0,
            fillcolor="lightgray", opacity=0.7, line_width=0, layer="below"
        )

        # ---- freestream lines ----
        self.fig.add_trace(go.Scatter(
            x=[x_min, -delta_g], y=[T_inf, T_inf],
            mode="lines",
            line=dict(color="crimson", width=2, dash="dash"),
            name="Gas freestream", showlegend=False
        ))
        self.fig.add_trace(go.Scatter(
            x=[wall_thickness + delta_c, x_max], y=[T_c, T_c],
            mode="lines",
            line=dict(color="crimson", width=2, dash="dash"),
            name="Coolant freestream", showlegend=False
        ))

        # ---- temperature profile ----
        self.fig.add_trace(go.Scatter(
            x=y_all, y=T_all,
            mode="lines",
            line=dict(color="crimson", width=2),
            name="Temperature profile"
        ))

        # ---- region labels ----
        y_min = float(np.nanmin([np.nanmin(T_all), T_inf, T_c]))
        y_max = float(np.nanmax([np.nanmax(T_all), T_inf, T_c]))
        y_mid = 0.5 * (y_min + y_max)

        self.fig.add_annotation(x=x_min, y=y_mid, xref="x", yref="y",
                                text="Freestream gas", showarrow=False,
                                xanchor="left", textangle=45)
        self.fig.add_annotation(x=-0.5 * delta_g, y=y_mid, xref="x", yref="y",
                                text="Gas BL", showarrow=False,
                                xanchor="center", textangle=45)
        self.fig.add_annotation(x=0.5 * wall_thickness, y=y_mid, xref="x", yref="y",
                                text="Wall", showarrow=False,
                                xanchor="center", textangle=45)
        self.fig.add_annotation(x=wall_thickness + 0.75 * delta_c, y=y_mid, xref="x", yref="y",
                                text="Coolant BL", showarrow=False,
                                xanchor="center", textangle=45)
        self.fig.add_annotation(x=x_max, y=y_mid, xref="x", yref="y",
                                text="Freestream coolant", showarrow=False,
                                xanchor="right", textangle=45)

        # ---- axes & layout ----
        self.fig.update_layout(
            title=f"Temperature profile at x = {x0:.3f} m, circuit {circuit_index}",
            xaxis=dict(title="Distance from hot-wall interface, y (m)"),
            yaxis=dict(title="Temperature, T (K)"),
            margin=dict(l=70, r=20, t=60, b=60),
        )
        self.fig.update_xaxes(range=[x_min, x_max])
