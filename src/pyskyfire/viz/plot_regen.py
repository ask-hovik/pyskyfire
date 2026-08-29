from __future__ import annotations
from typing import Sequence, Optional, Any, Iterable, Union
import warnings
import numpy as np
import plotly.graph_objects as go
from pyskyfire.regen import physics as regen_physics
from .core import PlotBase
from .plot_skycea import _visualization_x


def _result_grid(result, attribute: str) -> np.ndarray:
    """Return a result's native grid, falling back to legacy ``x``."""
    grid = getattr(result, attribute, None)
    return np.asarray(result.x if grid is None else grid, dtype=float)


class PlotContour(PlotBase):
    """
    Plot bell-nozzle or toroidal-aerospike contours, mirrored about the x-axis.

    Each contour can be either:
      • classic: has .xs and .rs
      • aerospike: has .xs_outer, .rs_outer, .xs_inner, .rs_inner

    Classic contours are sampled through their ``r(x)`` method when present.
    Set ``show_input_points=True`` to overlay the original input stations as
    ``x`` markers; ``interpolation_points`` controls the line resolution.
    If a contour has a .name attribute, it's used for the legend label.
    """

    def __init__(
        self,
        *contours: Any,
        show_labels: bool = True,
        show_input_points: bool = False,
        interpolation_points: int = 500,
        title: str = "Contour Profiles",
        template: str = "plotly_white",
    ):
        super().__init__(go.Figure())
        self.template(template)
        if interpolation_points < 2:
            raise ValueError("interpolation_points must be at least 2")

        # ------------ colorway helper (same color for +r and -r; rotate per contour) ------------
        #try:
        #    tpl = pio.templates[template] if template in pio.templates else pio.templates["plotly"]
        #    colorway = list(tpl.layout.colorway) if tpl.layout.colorway else []
        #except Exception:
        #    colorway = []

        #if not colorway:
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
                knots = np.asarray(contour.xs, dtype=float)
                if callable(getattr(contour, "r", None)):
                    xs = np.unique(
                        np.concatenate(
                            [
                                np.linspace(
                                    knots[0], knots[-1], interpolation_points
                                ),
                                knots,
                            ]
                        )
                    )
                    rs = np.asarray(contour.r(xs), dtype=float)
                else:
                    xs = knots
                    rs = np.asarray(contour.rs, dtype=float)

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

                if show_input_points:
                    input_xs = np.asarray(
                        getattr(contour, "input_xs", contour.xs),
                        dtype=float,
                    )
                    input_rs = np.asarray(
                        getattr(contour, "input_rs", contour.rs),
                        dtype=float,
                    )
                    marker = dict(color=color, symbol="x", size=8)
                    self.fig.add_trace(go.Scatter(
                        x=input_xs,
                        y=input_rs,
                        mode="markers",
                        name=f"{label} input points",
                        legendgroup=group,
                        showlegend=show_labels,
                        marker=marker,
                    ))
                    self.fig.add_trace(go.Scatter(
                        x=input_xs,
                        y=-input_rs,
                        mode="markers",
                        name=f"{label} input points (mirror)",
                        legendgroup=group,
                        showlegend=False,
                        marker=marker,
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


class PlotCoolantTemperature(PlotBase):
    """Plot bulk static coolant temperature from one or more ``RegenResult`` objects."""
    def __init__(self, *regen_results):
        super().__init__(go.Figure())
        self.template("plotly_white")

        for result in regen_results:
            x = _result_grid(result, "x_coolant")
            y = np.asarray(result.T_static)
            self.fig.add_trace(
                go.Scatter(x=x, y=y, mode="lines", name=result.circuit_name, showlegend=True)
            )

        self.fig.update_layout(
            title="Coolant Temperature",
            xaxis=dict(title="Axial Position (m)"),
            yaxis=dict(title="Coolant Temperature (K)"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )

class PlotCoolantPressure(PlotBase):
    """
    Plot static and/or stagnation coolant pressure for one or more circuits.

    Each input is a ``RegenResult``. The circuit name is used as the legend-group title.
    """

    def __init__(
        self,
        *regen_results,
        static: bool = True,
        stagnation: bool = True,
        template: str = "plotly_white",
    ):
        super().__init__(go.Figure())
        self.template(template)

        for result in regen_results:
            x = _result_grid(result, "x_coolant")
            circuit_name = result.circuit_name

            if static:
                self.fig.add_trace(
                    go.Scatter(
                        x=x,
                        y=np.asarray(result.p_static),
                        mode="lines",
                        name="Static Pressure",
                        legendgroup=circuit_name,
                        legendgrouptitle=dict(text=circuit_name),
                        showlegend=True,
                        line=dict(dash="solid"),
                    )
                )

            if stagnation:
                self.fig.add_trace(
                    go.Scatter(
                        x=x,
                        y=np.asarray(result.p_stagnation),
                        mode="lines",
                        name="Stagnation Pressure",
                        legendgroup=circuit_name,
                        legendgrouptitle=dict(text=circuit_name),
                        showlegend=True,
                        line=dict(dash="dash"),
                    )
                )

        self.fig.update_layout(
            title="Coolant Pressure",
            xaxis=dict(title="Axial Position (m)"),
            yaxis=dict(title="Pressure (Pa)"),
            legend=dict(
                title=None,
                tracegroupgap=10,
                groupclick="toggleitem",
            ),
            margin=dict(l=60, r=20, t=60, b=60),
        )

class PlotWallTemperature(PlotBase):
    """
    Build a wall-temperature plot from one or more ``RegenResult`` objects.

    ``T`` must have shape ``(len(x), n_layers)``. The circuit name is used as the
    legend-group title. Nodes whose final scaled wall-balance residual exceeds
    the solver tolerance are shown with red ``x`` markers by default.

    The optional layers extend the plot outwards from the wall stack, cold side
    first: ``plot_coolant_bulk`` (the bulk coolant the circuit marches on),
    ``plot_recovery`` (the hot-side recovery temperature ``T_aw_hot``, which is
    the driving potential behind the reported ``h_hot``), and
    ``plot_combustion`` (gas static and total temperature, which needs
    ``thrust_chamber``). Turning them all on gives the full cold-to-hot
    overview at the cost of stretching the y-axis over the whole range, which
    flattens the detail inside the wall stack.
    """

    def __init__(
        self,
        *regen_results,
        plot_hot: bool = True,
        plot_interfaces: bool = False,
        plot_coolant_wall: bool = False,
        plot_coolant_bulk: bool = False,
        plot_recovery: bool = False,
        plot_combustion: bool = False,
        thrust_chamber=None,
        mark_nonconverged: bool = True,
        template: str = "plotly_white",
    ):
        if plot_combustion and thrust_chamber is None:
            raise ValueError(
                "plot_combustion needs thrust_chamber to reach the combustion "
                "transport model"
            )

        super().__init__(go.Figure())
        self.template(template)

        for result in regen_results:
            x = _result_grid(result, "x_wall")
            T = np.asarray(result.T)
            circuit_name = result.circuit_name

            cols = []

            if plot_coolant_bulk:
                cols.append(0)

            if plot_coolant_wall and T.shape[1] > 1:
                cols.append(1)

            if plot_interfaces and T.shape[1] > 2:
                cols.extend(range(2, T.shape[1] - 1))

            if plot_hot:
                cols.append(T.shape[1] - 1)

            for col in cols:
                is_coolant_bulk = col == 0 and plot_coolant_bulk
                is_cold_wall = col == 1 and plot_coolant_wall
                is_hot_wall = col == T.shape[1] - 1 and plot_hot
                is_interface = plot_interfaces and 2 <= col <= T.shape[1] - 2

                if is_coolant_bulk:
                    trace_name = "Coolant Bulk"
                elif is_cold_wall:
                    trace_name = "Cold Wall"
                elif is_hot_wall:
                    trace_name = "Hot Wall"
                else:
                    trace_name = f"Interface {col - 1}"

                self.fig.add_trace(
                    go.Scatter(
                        x=x,
                        y=T[:, col],
                        mode="lines",
                        name=trace_name,
                        legendgroup=circuit_name,
                        legendgrouptitle=dict(text=circuit_name),
                        showlegend=True,
                        line=dict(dash="dash" if is_interface else "solid"),
                    )
                )

            T_aw = getattr(result, "T_aw_hot", None)
            if plot_recovery and T_aw is not None:
                self.fig.add_trace(
                    go.Scatter(
                        x=_result_grid(result, "x_heat_flux"),
                        y=np.asarray(T_aw, dtype=float),
                        mode="lines",
                        name="Recovery (T_aw)",
                        legendgroup=circuit_name,
                        legendgrouptitle=dict(text=circuit_name),
                        showlegend=True,
                        line=dict(dash="dot"),
                    )
                )

            if plot_combustion:
                transport = thrust_chamber.combustion_transport
                x_gas = _visualization_x(transport, x=None, results=result, nodes=None)
                series = (
                    ("Combustion Static", transport.get_T),
                    ("Combustion Total", transport.get_T0),
                )
                for trace_name, getter in series:
                    self.fig.add_trace(
                        go.Scatter(
                            x=x_gas,
                            y=np.asarray(
                                [getter(float(x_i)) for x_i in x_gas],
                                dtype=float,
                            ),
                            mode="lines",
                            name=trace_name,
                            legendgroup=circuit_name,
                            legendgrouptitle=dict(text=circuit_name),
                            showlegend=True,
                            line=dict(dash="dot"),
                        )
                    )

            converged = getattr(result, "wall_converged", None)
            if mark_nonconverged and cols and converged is not None:
                converged = np.asarray(converged, dtype=bool)
                if converged.shape != x.shape:
                    raise ValueError(
                        "result.wall_converged must have one value per wall node"
                    )
                failed = ~converged
                if np.any(failed):
                    residual = getattr(result, "wall_residual_scaled", None)
                    if residual is None:
                        residual = np.full(x.shape, np.nan)
                    residual = np.asarray(residual, dtype=float)
                    if residual.shape != x.shape:
                        raise ValueError(
                            "result.wall_residual_scaled must have one value "
                            "per wall node"
                        )
                    marker_col = T.shape[1] - 1 if plot_hot else cols[-1]
                    self.fig.add_trace(
                        go.Scatter(
                            x=x[failed],
                            y=T[failed, marker_col],
                            mode="markers",
                            name="Non-converged wall nodes",
                            legendgroup=circuit_name,
                            legendgrouptitle=dict(text=circuit_name),
                            showlegend=True,
                            marker=dict(symbol="x", size=11, color="#d62728"),
                            customdata=residual[failed],
                            hovertemplate=(
                                "x = %{x:.6g} m<br>"
                                "T = %{y:.2f} K<br>"
                                "scaled residual = %{customdata:.3e}"
                                "<extra></extra>"
                            ),
                        )
                    )

        self.fig.update_layout(
            title="Wall Temperature",
            xaxis=dict(title="Axial Position (m)"),
            yaxis=dict(title="Temperature (K)"),
            legend=dict(
                title=None,
                tracegroupgap=10,
                # Click a legend entry to toggle that trace alone rather than
                # its whole circuit group.
                groupclick="toggleitem",
            ),
            margin=dict(l=60, r=20, t=60, b=60),
        )

class PlotHeatFlux(PlotBase):
    """Plot wall heat flux from one or more ``RegenResult`` objects."""
    def __init__(self, *regen_results):
        super().__init__(go.Figure())
        self.template("plotly_white")

        for result in regen_results:
            x = _result_grid(result, "x_heat_flux")
            y = np.asarray(result.dQ_dA)
            self.fig.add_trace(
                go.Scatter(x=x, y=y, mode="lines", name=result.circuit_name, showlegend=True)
            )

        self.fig.update_layout(
            title="Heat Flux",
            xaxis=dict(title="Axial Position (m)"),
            yaxis=dict(title="Heat Flux (W/m²)"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )


class PlotVelocity(PlotBase):
    """Plot coolant velocity from one or more ``RegenResult`` objects."""
    def __init__(self, *regen_results):
        super().__init__(go.Figure())
        self.template("plotly_white")

        for result in regen_results:
            x = _result_grid(result, "x_coolant")
            y = np.asarray(result.velocity)
            self.fig.add_trace(
                go.Scatter(x=x, y=y, mode="lines", name=result.circuit_name, showlegend=True)
            )

        self.fig.update_layout(
            title="Coolant Velocity",
            xaxis=dict(title="Axial Position (m)"),
            yaxis=dict(title="Velocity (m/s)"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )


class PlotReynoldsNumber(PlotBase):
    """Plot combustion-gas and coolant Reynolds numbers along the chamber.

    Coolant profiles are evaluated on each result's coolant grid using the
    same bulk-state convention as the regenerative solver.  The laminar and
    turbulent limits from :mod:`pyskyfire.regen.physics` are drawn only when
    an individual profile spans the corresponding limit.

    Parameters
    ----------
    thrust_chamber : ThrustChamber
        Chamber carrying the combustion transport model, contour, and cooling
        circuits.
    *regen_results : RegenResult
        Solved cooling circuits to include.
    combustion, coolant : bool, optional
        Select which sides to display.
    combustion_nodes, combustion_x : optional
        Override the default combustion grid derived from ``regen_results``
        with a uniform node count or explicit axial coordinates.
    show_transition_thresholds : bool, optional
        Label crossed laminar/turbulent limits with horizontal lines.
    """

    def __init__(
        self,
        thrust_chamber,
        *regen_results,
        combustion: bool = True,
        coolant: bool = True,
        combustion_nodes: int | None = None,
        combustion_x=None,
        show_transition_thresholds: bool = True,
        template: str = "plotly_white",
    ):
        super().__init__(go.Figure())
        self.template(template)
        profiles = []

        if combustion:
            transport = thrust_chamber.combustion_transport
            x_gas = _visualization_x(
                transport,
                x=combustion_x,
                results=regen_results or None,
                nodes=combustion_nodes,
            )
            reynolds_gas = np.empty_like(x_gas)
            for index, x_i in enumerate(x_gas):
                state = transport.get_state(float(x_i))
                reynolds_gas[index] = regen_physics.reynolds(
                    float(state.rho),
                    float(state.M * state.a),
                    float(2.0 * thrust_chamber.contour.r(x_i)),
                    float(state.mu),
                )
            profiles.append(reynolds_gas)
            self.fig.add_trace(go.Scatter(
                x=x_gas,
                y=reynolds_gas,
                mode="lines",
                name="Combustion gas",
                showlegend=True,
            ))

        if coolant:
            for result in regen_results:
                x_coolant = _result_grid(result, "x_coolant")
                circuit = thrust_chamber.cooling_circuits[result.circuit_index]
                temperature = np.asarray(result.T_static, dtype=float)
                pressure = np.asarray(result.p_stagnation, dtype=float)
                velocity = np.asarray(result.velocity, dtype=float)
                quality = getattr(result, "coolant_quality", None)
                if quality is not None:
                    quality = np.asarray(quality, dtype=float)

                reynolds_coolant = np.empty_like(x_coolant)
                coolant_state = None
                if quality is not None and np.isfinite(quality).any():
                    from pyskyfire.regen.coolant_state import CoolantState

                    coolant_state = CoolantState(
                        circuit.coolant_transport, float(pressure[0])
                    )

                for index, x_i in enumerate(x_coolant):
                    T_i = float(temperature[index])
                    p_i = float(pressure[index])
                    if quality is not None and np.isfinite(quality[index]):
                        saturation = coolant_state.saturation(p_i)
                    else:
                        saturation = None
                    if saturation is not None:
                        rho_i = saturation.homogeneous_density(quality[index])
                        mu_i = saturation.mu_f
                    else:
                        rho_i = float(
                            circuit.coolant_transport.get_rho(T_i, p_i)
                        )
                        mu_i = float(
                            circuit.coolant_transport.get_mu(T_i, p_i)
                        )
                    reynolds_coolant[index] = regen_physics.reynolds(
                        rho_i,
                        float(velocity[index]),
                        float(circuit.Dh_coolant(x_i)),
                        mu_i,
                    )

                profiles.append(reynolds_coolant)
                self.fig.add_trace(go.Scatter(
                    x=x_coolant,
                    y=reynolds_coolant,
                    mode="lines",
                    name=f"{result.circuit_name} coolant",
                    showlegend=True,
                ))

        if show_transition_thresholds:
            thresholds = (
                (regen_physics.ReDh_laminar, "Laminar limit"),
                (regen_physics.ReDh_turbulent, "Turbulent limit"),
            )
            for threshold, label in thresholds:
                crossed = any(
                    np.isfinite(values).any()
                    and np.nanmin(values) <= threshold <= np.nanmax(values)
                    for values in profiles
                )
                if crossed:
                    self.fig.add_hline(
                        y=threshold,
                        line=dict(color="#666", dash="dash", width=1),
                        annotation_text=f"{label} (Re = {threshold:g})",
                        annotation_position="top left",
                    )

        self.fig.update_layout(
            title="Reynolds Number",
            xaxis=dict(title="Axial Position (m)"),
            yaxis=dict(title="Reynolds Number, Re"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )


# ---------------------------------------------------------------------------
# Coolant phase / boiling
# ---------------------------------------------------------------------------

#: Fixed categorical colour per coolant phase. Assigned in this order and never
#: cycled. Validated all-pairs for CVD on a light surface (worst deutan ΔE 9.2,
#: worst normal-vision ΔE 16.3); the aqua step sits below 3:1 against the
#: surface, which is why the phase strip carries visible labels as well as a
#: legend rather than relying on colour alone.
PHASE_COLORS = {
    "liquid": "#2a78d6",         # blue
    "two_phase": "#eb6834",      # orange -- the one you are looking for
    "vapour": "#1baf7a",         # aqua
    "supercritical": "#4a3aa7",  # violet; never co-occurs with the three above
    "unknown": "#8c8c8c",        # no saturation data (temperature march)
}

#: Legend/axis wording per phase.
PHASE_LABELS = {
    "liquid": "Subcooled liquid",
    "two_phase": "Two-phase",
    "vapour": "Superheated vapour",
    "supercritical": "Supercritical",
    "unknown": "Single-phase (no data)",
}

#: Draw order, so the legend reads the way the coolant actually evolves.
PHASE_ORDER = ("liquid", "two_phase", "vapour", "supercritical", "unknown")

_INK_SECONDARY = "#52514e"


def _phase_array(result) -> Optional[np.ndarray]:
    """Per-station phase labels for a ``RegenResult``, or ``None`` if absent.

    Results produced before the two-phase coolant march carry no phase data, so
    every plot here degrades to an explanatory annotation rather than raising.
    """
    phase = getattr(result, "coolant_phase", None)
    if phase is None:
        return None
    phase = np.asarray(phase, dtype=object).astype(str)
    return phase if phase.size else None


def _phase_segments(phase: np.ndarray):
    """Split a phase label array into contiguous runs.

    Yields
    ------
    tuple
        ``(label, i_start, i_end)`` with both indices inclusive.
    """
    if phase.size == 0:
        return
    start = 0
    for i in range(1, phase.size):
        if phase[i] != phase[start]:
            yield phase[start], start, i - 1
            start = i
    yield phase[start], start, phase.size - 1


def _segment_bounds(x: np.ndarray, i0: int, i1: int):
    """Axial extent of a segment, ordered low→high.

    Edges fall on the midpoints between neighbouring nodes, so adjacent
    segments share a boundary exactly instead of leaving a one-node hole -- the
    phase really does change somewhere inside that interval, and a visible gap
    would misread as a fourth state.

    Cooling circuits may march in either direction, so ``x`` is not guaranteed
    to increase and the raw endpoints cannot be used as ``x0``/``x1``.
    """
    n = x.size
    left = float(x[0]) if i0 == 0 else 0.5 * (float(x[i0 - 1]) + float(x[i0]))
    right = float(x[-1]) if i1 == n - 1 else 0.5 * (float(x[i1]) + float(x[i1 + 1]))
    return (left, right) if left <= right else (right, left)


class PlotCoolantQuality(PlotBase):
    """Plot bulk vapour quality along the circuit from ``RegenResult`` objects.

    Quality is only defined inside the boiling dome, so the solved values are
    drawn as a solid line and the subcooled/superheated stretches either side
    are carried at 0 and 1 by a dashed line. The dashed parts are bookkeeping,
    not model output -- the distinction is the point of the two line styles.

    Parameters
    ----------
    *regen_results : RegenResult
        One or more solved circuits.
    shade_dome : bool, optional
        Shade the axial extent of the two-phase region. Default True.
    template : str, optional
        Plotly template name.

    Notes
    -----
    A quality profile says *where* the coolant boils, not how hot the wall gets
    there: the coolant-side coefficient is still a single-phase correlation
    inside the dome. See :mod:`pyskyfire.regen.coolant_state`.
    """

    def __init__(self, *regen_results, shade_dome: bool = True,
                 template: str = "plotly_white"):
        super().__init__(go.Figure())
        self.template(template)

        multi = len(regen_results) > 1
        any_quality = False
        any_boiling = False
        all_supercritical = True

        for result in regen_results:
            x = _result_grid(result, "x_coolant")
            quality = getattr(result, "coolant_quality", None)
            phase = _phase_array(result)
            if quality is None or phase is None:
                continue
            quality = np.asarray(quality, dtype=float)
            any_quality = True

            name = result.circuit_name
            group = dict(legendgroup=name, legendgrouptitle=dict(text=name)) if multi else {}

            # Carry the curve at 0 through the subcooled run and at 1 past the
            # dome. Supercritical and unknown stations are left blank instead:
            # there is no dome to be on either side of, so drawing 0 there would
            # assert "liquid, not boiling" when the truth is "not a quantity".
            has_dome = np.isin(phase, ("liquid", "two_phase", "vapour"))
            all_supercritical = all_supercritical and not has_dome.any()
            carried = np.where(phase == "vapour", 1.0, 0.0)
            carried = np.where(np.isfinite(quality), quality, carried)
            carried = np.where(has_dome, carried, np.nan)
            in_dome = phase == "two_phase"
            any_boiling = any_boiling or bool(in_dome.any())

            if has_dome.any():
                self.fig.add_trace(go.Scatter(
                    x=x, y=carried,
                    mode="lines",
                    name="Outside dome (carried)" if not multi else f"{name} (carried)",
                    line=dict(color=PHASE_COLORS["liquid"], width=2, dash="dot"),
                    hovertemplate="x = %{x:.4g} m<br>quality = %{y:.3f}<extra></extra>",
                    showlegend=True,
                    **group,
                ))

            if in_dome.any():
                solved = np.where(in_dome, quality, np.nan)
                self.fig.add_trace(go.Scatter(
                    x=x, y=solved,
                    mode="lines",
                    name="Vapour quality" if not multi else f"{name} (dome)",
                    line=dict(color=PHASE_COLORS["two_phase"], width=2),
                    hovertemplate="x = %{x:.4g} m<br>quality = %{y:.3f}<extra></extra>",
                    showlegend=True,
                    **group,
                ))

            if shade_dome:
                for label, i0, i1 in _phase_segments(phase):
                    if label != "two_phase":
                        continue
                    x0, x1 = _segment_bounds(x, i0, i1)
                    self.fig.add_vrect(
                        x0=x0, x1=x1,
                        fillcolor=PHASE_COLORS["two_phase"], opacity=0.10,
                        line_width=0, layer="below",
                    )

            # Direct-label the onset and the exit state rather than every point.
            onset = np.flatnonzero(in_dome)
            if onset.size:
                i0 = int(onset[0])
                # Offset into the subcooled side, which is empty by
                # construction -- centring the label on the onset would sit it
                # on top of the quality curve.
                self.fig.add_annotation(
                    x=float(x[i0]), y=0.0,
                    text=f"boiling onset<br>x = {x[i0]:.4g} m",
                    showarrow=True, arrowhead=0, arrowwidth=1,
                    ax=60, ay=-28,
                    font=dict(color=_INK_SECONDARY, size=11),
                )
                if np.isfinite(carried[-1]):
                    self.fig.add_annotation(
                        x=float(x[-1]), y=float(carried[-1]),
                        text=f"exit {carried[-1]:.3f}",
                        showarrow=False, xanchor="left", yshift=12,
                        font=dict(color=_INK_SECONDARY, size=11),
                    )

        if not any_quality:
            note = ("No quality data on these results.<br>"
                    "The coolant ran on the temperature march "
                    "(mixture or no saturation data).")
        elif all_supercritical:
            note = ("Coolant is supercritical along the whole circuit.<br>"
                    "No boiling dome exists, so quality is undefined.")
        elif not any_boiling:
            note = "Coolant stays single-phase along the whole circuit."
        else:
            note = None

        if note is not None:
            self.fig.add_annotation(
                text=note,
                xref="paper", yref="paper", x=0.5, y=0.5,
                showarrow=False, align="center",
                font=dict(color=_INK_SECONDARY),
            )

        # With no traces at all -- the supercritical case -- Plotly would fall
        # back to a default range and the axis would no longer describe the
        # engine, so pin it to the circuit whatever happened above.
        xaxis = dict(title="Axial Position (m)")
        if not self.fig.data and regen_results:
            spans = [_result_grid(r, "x_coolant") for r in regen_results]
            xaxis["range"] = [min(s.min() for s in spans),
                              max(s.max() for s in spans)]

        self.fig.update_layout(
            title="Coolant Vapour Quality",
            # Deliberately not "x": that is the axial coordinate on this figure.
            xaxis=xaxis,
            yaxis=dict(title="Vapour quality (-)", range=[-0.05, 1.05]),
            legend=dict(title=None, tracegroupgap=10, groupclick="toggleitem"),
            margin=dict(l=60, r=20, t=60, b=60),
        )


class PlotCoolantPhase(PlotBase):
    """Plot the coolant phase regime as a band along the circuit.

    One horizontal strip per circuit, split into contiguous phase regions. This
    is the at-a-glance companion to :class:`PlotCoolantQuality`: it answers
    "does this coolant boil, and where" without reading a curve.

    Parameters
    ----------
    *regen_results : RegenResult
        One or more solved circuits.
    label_bands : bool, optional
        Write the phase name inside bands wide enough to hold it. Default True.
    template : str, optional
        Plotly template name.
    """

    def __init__(self, *regen_results, label_bands: bool = True,
                 template: str = "plotly_white"):
        super().__init__(go.Figure())
        self.template(template)

        rows = []
        legend_done = set()
        x_span = None

        for result in regen_results:
            phase = _phase_array(result)
            if phase is None:
                continue
            rows.append((result, phase))

        # Comparing one circuit under two conditions is the common case, and
        # then every row carries the same circuit name. Disambiguate rather
        # than stacking identical tick labels.
        names = [r.circuit_name for r, _ in rows]
        row_labels = [
            f"{n} ({i + 1})" if names.count(n) > 1 else n
            for i, n in enumerate(names)
        ]

        for row, (result, phase) in enumerate(rows):
            x = _result_grid(result, "x_coolant")
            quality = np.asarray(
                getattr(result, "coolant_quality", np.full(x.size, np.nan)),
                dtype=float,
            )
            y0, y1 = row + 0.15, row + 0.85
            lo, hi = float(np.min(x)), float(np.max(x))
            x_span = (lo, hi) if x_span is None else (min(x_span[0], lo), max(x_span[1], hi))

            for label, i0, i1 in _phase_segments(phase):
                xa, xb = _segment_bounds(x, i0, i1)
                color = PHASE_COLORS.get(label, PHASE_COLORS["unknown"])

                q_seg = quality[i0:i1 + 1]
                q_text = (
                    f"<br>quality {np.nanmin(q_seg):.3f} – {np.nanmax(q_seg):.3f}"
                    if np.isfinite(q_seg).any() else ""
                )

                self.fig.add_trace(go.Scatter(
                    x=[xa, xb, xb, xa, xa],
                    y=[y0, y0, y1, y1, y0],
                    mode="lines",
                    fill="toself",
                    fillcolor=color,
                    # A surface-coloured ring gives the 2 px gap between
                    # touching bands without displacing either boundary.
                    line=dict(color="#ffffff", width=2),
                    name=PHASE_LABELS.get(label, label),
                    legendgroup=label,
                    showlegend=label not in legend_done,
                    # Legend order follows the physical progression, not the
                    # order the circuit happens to be marched in -- a circuit
                    # running -x would otherwise list vapour first.
                    legendrank=(PHASE_ORDER.index(label)
                                if label in PHASE_ORDER else len(PHASE_ORDER)),
                    hoveron="fills",
                    hoverinfo="text",
                    text=(f"{row_labels[row]}<br>"
                          f"{PHASE_LABELS.get(label, label)}<br>"
                          f"x = {xa:.4g} … {xb:.4g} m{q_text}"),
                ))
                legend_done.add(label)

                # Relief for the low-contrast steps: name the band in place
                # whenever it is wide enough, so identity never rests on hue.
                if label_bands and x_span is not None:
                    total = max(x_span[1] - x_span[0], 1e-12)
                    if (xb - xa) / total > 0.10:
                        # The y axis is reversed (row 0 on top), so the band's
                        # low-y edge is its top edge on screen.
                        self.fig.add_annotation(
                            x=0.5 * (xa + xb), y=row + 0.15,
                            text=PHASE_LABELS.get(label, label),
                            showarrow=False, yanchor="bottom",
                            font=dict(color=_INK_SECONDARY, size=11),
                        )

        if not rows:
            self.fig.add_annotation(
                text=("No phase data on these results.<br>"
                      "The coolant ran on the temperature march "
                      "(mixture or no saturation data)."),
                xref="paper", yref="paper", x=0.5, y=0.5,
                showarrow=False, align="center",
                font=dict(color=_INK_SECONDARY),
            )

        self.fig.update_layout(
            title="Coolant Phase",
            xaxis=dict(title="Axial Position (m)"),
            yaxis=dict(
                title=None,
                tickmode="array",
                tickvals=[i + 0.5 for i in range(len(rows))],
                ticktext=row_labels,
                # Reversed so the first result passed reads at the top; the
                # extra headroom at the low end holds the band labels.
                range=[max(len(rows), 1) + 0.05, -0.45],
                showgrid=False, zeroline=False,
            ),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )


class PlotSaturationMargin(PlotBase):
    """Plot bulk coolant temperature against its saturation temperature.

    The subcooling margin ``T_sat - T_bulk`` is what decides whether the coolant
    is anywhere near boiling, and it is the quantity to watch on a circuit that
    never enters the dome -- where :class:`PlotCoolantQuality` has nothing to
    show. Both curves are temperatures, so they share one axis; the gap between
    them is shaded.

    Parameters
    ----------
    *regen_results : RegenResult
        One or more solved circuits.
    template : str, optional
        Plotly template name.

    Notes
    -----
    A positive margin does not mean no boiling at the wall. Subcooled nucleate
    boiling starts once the *wall* passes ``T_sat`` while the bulk is still
    below it, and that is not modelled -- so this plot bounds the bulk
    behaviour only.
    """

    def __init__(self, *regen_results, template: str = "plotly_white"):
        super().__init__(go.Figure())
        self.template(template)

        multi = len(regen_results) > 1
        plotted = False

        for result in regen_results:
            x = _result_grid(result, "x_coolant")
            T_bulk = np.asarray(result.T_static, dtype=float)
            T_sat = getattr(result, "coolant_T_sat", None)
            phase = _phase_array(result)
            name = result.circuit_name
            group = dict(legendgroup=name, legendgrouptitle=dict(text=name)) if multi else {}

            if T_sat is not None and phase is not None:
                T_sat = np.asarray(T_sat, dtype=float)
                # Shade the margin only where the bulk is subcooled liquid.
                # Past saturation the same gap is *superheat*, not margin, and
                # filling it would read as "more headroom" exactly where there
                # is none left.
                liquid = (phase == "liquid") & np.isfinite(T_sat)
                if liquid.any():
                    self.fig.add_trace(go.Scatter(
                        x=x, y=np.where(liquid, T_bulk, np.nan),
                        mode="lines", line=dict(width=0),
                        hoverinfo="skip", showlegend=False, **group,
                    ))
                    self.fig.add_trace(go.Scatter(
                        x=x, y=np.where(liquid, T_sat, np.nan),
                        mode="lines", line=dict(width=0),
                        fill="tonexty", fillcolor="rgba(235,104,52,0.10)",
                        hoverinfo="skip", showlegend=False, **group,
                    ))

            self.fig.add_trace(go.Scatter(
                x=x, y=T_bulk,
                mode="lines",
                name="Bulk coolant" if not multi else f"{name} bulk",
                line=dict(color=PHASE_COLORS["liquid"], width=2),
                hovertemplate="x = %{x:.4g} m<br>T = %{y:.1f} K<extra></extra>",
                showlegend=True,
                **group,
            ))

            if T_sat is None or not np.isfinite(np.asarray(T_sat, dtype=float)).any():
                continue
            T_sat = np.asarray(T_sat, dtype=float)
            plotted = True

            self.fig.add_trace(go.Scatter(
                x=x, y=T_sat,
                mode="lines",
                name="Saturation" if not multi else f"{name} T_sat",
                line=dict(color=PHASE_COLORS["two_phase"], width=2, dash="dash"),
                hovertemplate="x = %{x:.4g} m<br>T_sat = %{y:.1f} K<extra></extra>",
                showlegend=True,
                **group,
            ))

            # The minimum margin is only meaningful over the subcooled run; if
            # the coolant reaches saturation the margin is zero by definition
            # and the useful fact is where that happened.
            liquid = (phase == "liquid") & np.isfinite(T_sat) if phase is not None else np.isfinite(T_sat)
            boils = (phase == "two_phase").any() if phase is not None else False
            if boils:
                i = int(np.flatnonzero(phase == "two_phase")[0])
                text = f"reaches saturation<br>x = {x[i]:.4g} m"
            elif liquid.any():
                margin = np.where(liquid, T_sat - T_bulk, np.nan)
                i = int(np.nanargmin(margin))
                text = f"min subcooling {margin[i]:.1f} K"
            else:
                continue

            self.fig.add_annotation(
                x=float(x[i]), y=float(T_sat[i]),
                text=text,
                showarrow=True, arrowhead=0, arrowwidth=1,
                ax=0, ay=-36,
                font=dict(color=_INK_SECONDARY, size=11),
            )

        if not plotted:
            self.fig.add_annotation(
                text=("No saturation temperature available.<br>"
                      "The coolant is supercritical, or ran on the "
                      "temperature march."),
                xref="paper", yref="paper", x=0.98, y=0.05,
                showarrow=False, align="right", xanchor="right",
                font=dict(color=_INK_SECONDARY),
            )

        self.fig.update_layout(
            title="Coolant Subcooling Margin",
            xaxis=dict(title="Axial Position (m)"),
            yaxis=dict(title="Temperature (K)"),
            legend=dict(title=None, tracegroupgap=10, groupclick="toggleitem"),
            margin=dict(l=60, r=20, t=60, b=60),
        )


IndexLike = Optional[Union[int, Iterable[int]]]

def _normalize_indices(thrust_chamber, circuit_index: IndexLike):
    circuits = thrust_chamber.cooling_circuits
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

        circuits = thrust_chamber.cooling_circuits
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

        circuits = thrust_chamber.cooling_circuits
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

        circuits = thrust_chamber.cooling_circuits
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


class PlotChannelHeight(PlotBase):
    """Plot the raw cooling-channel height profile for selected circuits."""

    def __init__(self, thrust_chamber, circuit_index: IndexLike = None):
        super().__init__(go.Figure())
        self.template("plotly_white")

        circuits = thrust_chamber.cooling_circuits
        indices = _normalize_indices(thrust_chamber, circuit_index)

        for idx in indices:
            circuit = circuits[idx]
            x = np.asarray(circuit.x_domain, dtype=float)
            height = np.asarray(circuit.channel_heights, dtype=float)
            name = getattr(circuit, "name", f"Circuit {idx}")
            self.fig.add_trace(
                go.Scatter(
                    x=x,
                    y=height,
                    mode="lines",
                    name=name,
                    showlegend=True,
                    hovertemplate=(
                        "x = %{x:.4g} m<br>h = %{y:.4g} m<extra></extra>"
                    ),
                )
            )

        self.fig.update_layout(
            title="Cooling Channel Height",
            xaxis=dict(title="Axial Position, x (m)"),
            yaxis=dict(title="Channel Height, h (m)"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )


class PlotChannelWidth(PlotBase):
    """Compare calculated tube OD width with optional reference profiles.

    The rounded-section geometry occupies an angular sector ``theta`` at the
    hot-side tube radius ``r``.  Its measurable tube width there is therefore
    the chord ``2 r sin(theta/2)``.  Reference profiles are expected to expose
    ``x`` and ``tube_width_od`` arrays and may optionally provide ``name``.
    """

    def __init__(
        self,
        thrust_chamber,
        circuit_index: IndexLike = None,
        reference_profiles: Optional[Sequence[Any]] = None,
    ):
        super().__init__(go.Figure())
        self.template("plotly_white")

        circuits = thrust_chamber.cooling_circuits
        indices = _normalize_indices(thrust_chamber, circuit_index)

        for idx in indices:
            circuit = circuits[idx]
            x = np.asarray(circuit.x_domain, dtype=float)
            radius = np.asarray(circuit.centerlines[0][:, 1], dtype=float)
            theta = np.asarray(circuit.channel_local_sector, dtype=float)
            width = 2.0 * radius * np.sin(0.5 * theta)
            name = getattr(circuit, "name", f"Circuit {idx}")
            self.fig.add_trace(
                go.Scatter(
                    x=x,
                    y=width,
                    mode="lines",
                    name=f"{name} (calculated)",
                    showlegend=True,
                    hovertemplate=(
                        "x = %{x:.4g} m<br>tube width OD = %{y:.4g} m"
                        "<extra></extra>"
                    ),
                )
            )

        for index, reference in enumerate(reference_profiles or ()):
            x = np.asarray(reference.x, dtype=float)
            width = np.asarray(reference.tube_width_od, dtype=float)
            if x.shape != width.shape:
                raise ValueError("reference tube-width x and y arrays must match")
            name = getattr(reference, "name", f"Reference {index + 1}")
            self.fig.add_trace(
                go.Scatter(
                    x=x,
                    y=width,
                    mode="markers",
                    name=name,
                    marker=dict(
                        symbol="circle-open",
                        size=7,
                        color="black",
                    ),
                    showlegend=True,
                    hovertemplate=(
                        "x = %{x:.4g} m<br>tube width OD = %{y:.4g} m"
                        "<extra></extra>"
                    ),
                )
            )

        self.fig.update_layout(
            title="Cooling Tube Width",
            xaxis=dict(title="Axial Position, x (m)"),
            yaxis=dict(title="Tube Width OD (m)"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )


# Categorical slots taken in fixed order, coolant first. Adjacent bands in the
# radial stack are the pairs that have to separate, and every sequence this
# assignment can produce clears the CVD and normal-vision floors.
_COOLANT_LAYER_COLOR = "#2a78d6"
_SOLID_LAYER_COLORS = (
    "#eb6834",  # orange
    "#1baf7a",  # aqua
    "#eda100",  # yellow
    "#e87ba4",  # magenta
    "#008300",  # green
    "#4a3aa7",  # violet
    "#e34948",  # red
)


def _layer_label(wall, index: int) -> str:
    """Name a wall layer from its own label, else its material, else its index."""
    name = getattr(wall, "name", None)
    if name:
        return str(name)
    material = getattr(wall, "material", None)
    material_name = getattr(material, "name", None)
    if material_name:
        return str(material_name)
    return f"Layer {index + 1}"


class PlotWallLayers(PlotBase):
    """Draw the wall stack of each cooling circuit to scale in the meridional plane.

    Where :class:`PlotContour` draws the hot-gas surface as a line, this fills
    the physical radial extent of every layer between the hot wall and the
    closeout, so zooming in reveals the real thicknesses. Bands are stacked from
    the hot surface outward exactly as the solver stacks them: the circuit's
    walls in order, then the channel height, then the closeout.

    Two things drawn here are not modelled by pyskyfire and are shown as such:

    * **The closeout wall.** The wall on the atmosphere side of the channel has
      no definition anywhere in the library, so it is drawn as a copy of the
      coolant-facing layer and hatched to mark it as an assumption. Replace it
      once a real definition exists (e.g. for strength calculations).
    * **Ribs.** The channel band is the full circumferential sector, so the
      land between channels is not resolved; the band shows where coolant
      lives radially, not the channel-to-rib split.

    Helical circuits are drawn as if the helix angle were zero -- the
    meridional path rotationally projected onto the plane. A true rendering
    would need a section through the 3-D engine. A warning is emitted when a
    helical circuit is plotted.

    Parameters
    ----------
    thrust_chamber : ThrustChamber
        Chamber whose cooling circuits supply the wall stacks.
    circuit_index : int or iterable of int, optional
        Circuits to draw. Defaults to all of them. Circuits that overlap
        axially are drawn on top of each other, so pick one when the chamber
        interlaces several passes.
    interpolation_points : int, optional
        Axial resolution of the filled bands. The default is generous because
        the point of the figure is that it survives being zoomed into.
    show_contour : bool, optional
        Overlay the full hot-gas contour, including where no circuit runs.
    mirror : bool, optional
        Draw the stack on both sides of the axis.
    """

    def __init__(
        self,
        thrust_chamber,
        circuit_index: IndexLike = None,
        *,
        interpolation_points: int = 1500,
        show_contour: bool = True,
        mirror: bool = True,
        title: str = "Wall Layers (to scale)",
        template: str = "plotly_white",
    ):
        super().__init__(go.Figure())
        self.template(template)
        if interpolation_points < 2:
            raise ValueError("interpolation_points must be at least 2")

        circuits = thrust_chamber.cooling_circuits
        indices = _normalize_indices(thrust_chamber, circuit_index)

        for idx in indices:
            circuit = circuits[idx]
            circuit_name = getattr(circuit, "name", f"Circuit {idx}")

            if getattr(circuit, "is_helical", False):
                warnings.warn(
                    f"Cooling circuit {circuit_name!r} is helical. Its wall "
                    "layers are drawn as the meridional path rotationally "
                    "projected onto the plane, i.e. as if the helix angle were "
                    "zero, so azimuthal extent and true section shape are not "
                    "represented.",
                    RuntimeWarning,
                    stacklevel=2,
                )

            x = self._sample_grid(circuit, interpolation_points)

            # Radial stack, hot surface outward, starting from the same
            # placement radius ``build_channel_centerlines`` builds on, so the
            # drawing agrees with the geometry the solver ran.
            r = np.array(
                [
                    float(
                        circuit.placement.compute_centerline_radius(
                            float(x_i), circuit.contour
                        )
                    )
                    for x_i in x
                ],
                dtype=float,
            )

            walls = list(circuit.walls)
            for layer_index, wall in enumerate(walls):
                thickness = np.asarray(wall.thickness(x), dtype=float)
                r_next = r + thickness
                self._add_band(
                    x,
                    r,
                    r_next,
                    name=_layer_label(wall, layer_index),
                    color=_SOLID_LAYER_COLORS[
                        layer_index % len(_SOLID_LAYER_COLORS)
                    ],
                    group=circuit_name,
                    mirror=mirror,
                )
                r = r_next

            height = np.array(
                [float(circuit.channel_height(float(x_i))) for x_i in x],
                dtype=float,
            )
            r_next = r + height
            self._add_band(
                x,
                r,
                r_next,
                name="Coolant channel",
                color=_COOLANT_LAYER_COLOR,
                group=circuit_name,
                mirror=mirror,
            )
            r = r_next

            if walls:
                closeout = walls[-1]
                closeout_index = len(walls) - 1
                thickness = np.asarray(closeout.thickness(x), dtype=float)
                self._add_band(
                    x,
                    r,
                    r + thickness,
                    name=f"{_layer_label(closeout, closeout_index)} (closeout, assumed)",
                    color=_SOLID_LAYER_COLORS[
                        closeout_index % len(_SOLID_LAYER_COLORS)
                    ],
                    group=circuit_name,
                    mirror=mirror,
                    pattern=True,
                )

        if show_contour:
            contour = thrust_chamber.contour
            knots = np.asarray(contour.xs, dtype=float)
            x_contour = np.linspace(
                float(knots[0]), float(knots[-1]), interpolation_points
            )
            r_contour = np.asarray(contour.r(x_contour), dtype=float)
            self.fig.add_trace(
                go.Scatter(
                    x=x_contour,
                    y=r_contour,
                    mode="lines",
                    name="Hot-gas contour",
                    legendgroup="contour",
                    line=dict(color="#52514e", width=1.5),
                    hovertemplate=(
                        "x = %{x:.4g} m<br>r = %{y:.4g} m<extra></extra>"
                    ),
                )
            )
            if mirror:
                self.fig.add_trace(
                    go.Scatter(
                        x=x_contour,
                        y=-r_contour,
                        mode="lines",
                        name="Hot-gas contour (mirror)",
                        legendgroup="contour",
                        showlegend=False,
                        line=dict(color="#52514e", width=1.5),
                        hoverinfo="skip",
                    )
                )

        self.fig.update_layout(
            title=title or None,
            xaxis=dict(title="Axial Position, x (m)"),
            yaxis=dict(title="Radius, r (m)"),
            legend=dict(title=None, tracegroupgap=10, groupclick="toggleitem"),
            margin=dict(l=60, r=20, t=60, b=60),
        )
        # To scale in both directions, so the layers keep their real proportions
        # as the reader zooms.
        self.fig.update_yaxes(scaleanchor="x", scaleratio=1)

    @staticmethod
    def _sample_grid(circuit, interpolation_points: int) -> np.ndarray:
        """Refine the circuit's own domain so thin layers stay smooth."""
        x_domain = np.asarray(circuit.x_domain, dtype=float)
        x_start = float(np.min(x_domain))
        x_end = float(np.max(x_domain))
        return np.unique(
            np.concatenate(
                [np.linspace(x_start, x_end, interpolation_points), x_domain]
            )
        )

    def _add_band(
        self,
        x: np.ndarray,
        r_inner: np.ndarray,
        r_outer: np.ndarray,
        *,
        name: str,
        color: str,
        group: str,
        mirror: bool,
        pattern: bool = False,
    ) -> None:
        """Fill one layer as a closed polygon, mirrored about the axis."""
        # ``None`` splits a single toself trace into two polygons, which keeps
        # both sides of the axis on one legend entry.
        polygon_x = np.concatenate([x, x[::-1]])
        polygon_y = np.concatenate([r_inner, r_outer[::-1]])
        if mirror:
            polygon_x = np.concatenate([polygon_x, [np.nan], polygon_x])
            polygon_y = np.concatenate([polygon_y, [np.nan], -polygon_y])

        trace = go.Scatter(
            x=polygon_x,
            y=polygon_y,
            mode="lines",
            fill="toself",
            fillcolor=color,
            # No stroke and no inter-band spacer: this is a to-scale drawing,
            # and a line of fixed pixel width would swallow the thin layers it
            # is meant to separate. Hue carries the split instead.
            line=dict(width=0),
            name=name,
            legendgroup=group,
            legendgrouptitle=dict(text=group),
            hoveron="fills",
            hoverinfo="text",
            text=f"{group} — {name}",
        )
        if pattern:
            trace.fillpattern = dict(
                shape="/", fgcolor="#fcfcfb", bgcolor=color, size=7, solidity=0.25
            )
        self.fig.add_trace(trace)


class PlotdAdxCoolantArea(PlotBase):
    def __init__(self, thrust_chamber, circuit_index: IndexLike = None):
        super().__init__(go.Figure())
        self.template("plotly_white")

        circuits = thrust_chamber.cooling_circuits
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

        circuits = thrust_chamber.cooling_circuits
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

        circuits = thrust_chamber.cooling_circuits
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
      - result: ``RegenResult`` containing ``x``, ``T`` (N×n_layers), ``T_static``,
        ``T_aw_hot``, ``h_hot``, ``h_cold``, and ``p_static``
      - thrust_chamber: provides combustion_transport and cooling circuits
      - circuit_index: which cooling circuit to use for coolant props
      - x_query: axial location (m) to sample
    """

    def __init__(self, result, thrust_chamber, circuit_index: int, x_query: float, n_bl: int = 1000):
        super().__init__(go.Figure())
        self.template("plotly_white")

        # ---- 1) nearest axial node ----
        x_arr = np.asarray(result.x, dtype=float)
        i = int(np.argmin(np.abs(x_arr - x_query)))
        x0 = float(x_arr[i])

        # ---- 2) extract temps & heat flux ----
        combustion = thrust_chamber.combustion_transport
        gas = combustion.get_state(x0)
        x_heat = _result_grid(result, "x_heat_flux")
        T_gas_edge = float(
            np.interp(
                x0,
                np.sort(x_heat),
                np.asarray(result.T_aw_hot)[np.argsort(x_heat)],
            )
        )
        T_hw  = float(result.T[i, -1])    # hot-wall
        T_cw  = float(result.T[i,  1])    # coolant-side wall
        x_cool = _result_grid(result, "x_coolant")
        T_c = float(
            np.interp(
                x0,
                np.sort(x_cool),
                np.asarray(result.T_static)[np.argsort(x_cool)],
            )
        )
        h_g = float(
            np.interp(
                x0,
                np.sort(x_heat),
                np.asarray(result.h_hot)[np.argsort(x_heat)],
            )
        )

        # ---- 3) gas-side turbulent thermal boundary layer (Kader) ----
        velocity_g = float(gas.M * gas.a)
        diameter_g = float(2.0 * thrust_chamber.contour.r(x0))
        reynolds_g = regen_physics.reynolds(
            float(gas.rho), velocity_g, diameter_g, float(gas.mu)
        )
        darcy_g = regen_physics.f_darcy_turbulent(
            reynolds_g, diameter_g, x0, roughness=None
        )
        u_tau_g = regen_physics.friction_velocity(velocity_g, darcy_g)
        distance_g, temperature_g, delta_g = regen_physics.kader_temperature_profile(
            T_wall=T_hw,
            T_edge=T_gas_edge,
            h=h_g,
            rho=float(gas.rho),
            cp=float(gas.cp),
            mu=float(gas.mu),
            prandtl=float(gas.Pr),
            u_tau=u_tau_g,
            n_points=n_bl,
        )
        y_g = -distance_g[::-1]
        T_g = temperature_g[::-1]

        # ---- 4) wall layers (hot→coolant) ----
        Ts_rev   = result.T[i, 1:]       # coolant-side → hot-side
        Ts_wall  = Ts_rev[::-1]              # hot-side → coolant-side (interfaces included)
        walls    = thrust_chamber.cooling_circuits[circuit_index].walls
        thicknesses = [float(w.thickness(x0)) for w in walls]
        y_w = np.insert(np.cumsum(thicknesses), 0, 0.0) if thicknesses else np.array([0.0])
        wall_thickness = float(y_w[-1])

        # ---- 5) coolant-side turbulent thermal boundary layer (Kader) ----
        p_static = float(
            np.interp(
                x0,
                np.sort(x_cool),
                np.asarray(result.p_static)[np.argsort(x_cool)],
            )
        )
        T_film   = 0.5 * (T_c + T_cw)
        circuit = thrust_chamber.cooling_circuits[circuit_index]
        coolant = circuit.coolant_transport
        h_c = float(
            np.interp(
                x0,
                np.sort(x_cool),
                np.asarray(result.h_cold)[np.argsort(x_cool)],
            )
        )
        velocity_c = float(
            np.interp(
                x0,
                np.sort(x_cool),
                np.asarray(result.velocity)[np.argsort(x_cool)],
            )
        )
        rho_bulk_c = float(coolant.get_rho(T_c, p_static))
        mu_bulk_c = float(coolant.get_mu(T_c, p_static))
        diameter_c = float(circuit.Dh_coolant(x0))
        reynolds_c = regen_physics.reynolds(
            rho_bulk_c, velocity_c, diameter_c, mu_bulk_c
        )
        darcy_c = regen_physics.f_darcy(
            reynolds_c, diameter_c, x0, circuit.roughness
        )
        darcy_c *= regen_physics.phi_curv_friction(
            reynolds_c, diameter_c, circuit.radius_of_curvature(x0)
        )
        u_tau_c = regen_physics.friction_velocity(velocity_c, darcy_c)
        distance_c, T_cBL, delta_c = regen_physics.kader_temperature_profile(
            T_wall=T_cw,
            T_edge=T_c,
            h=h_c,
            rho=float(coolant.get_rho(T_film, p_static)),
            cp=float(coolant.get_cp(T_film, p_static)),
            mu=float(coolant.get_mu(T_film, p_static)),
            prandtl=float(coolant.get_Pr(T_film, p_static)),
            u_tau=u_tau_c,
            n_points=n_bl,
        )
        y_c = distance_c

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
        material_keys = [
            getattr(wall.material, "name", None)
            or f"{type(wall.material).__name__}-{id(wall.material)}"
            for wall in walls
        ]
        material_labels = [
            getattr(wall, "name", None)
            or getattr(wall.material, "name", None)
            or f"Wall material {index + 1}"
            for index, wall in enumerate(walls)
        ]
        unique_materials = list(dict.fromkeys(material_keys))
        wall_palette = [
            "#7f8c8d", "#d4a72c", "#b87333", "#6c8ebf", "#9b59b6",
            "#5d8a66", "#c96b6b", "#607d8b",
        ]
        material_colors = {
            material: (
                "gray" if len(unique_materials) == 1
                else wall_palette[index % len(wall_palette)]
            )
            for index, material in enumerate(unique_materials)
        }
        for index, material in enumerate(material_keys):
            self.fig.add_shape(
                type="rect", xref="x", yref="paper",
                x0=float(y_w[index]), x1=float(y_w[index + 1]),
                y0=0.0, y1=1.0,
                fillcolor=material_colors[material], opacity=0.8,
                line_width=0, layer="below",
                name=material_labels[index],
                showlegend=(
                    len(unique_materials) > 1
                    and material_keys.index(material) == index
                ),
            )
        self.fig.add_shape(
            type="rect", xref="x", yref="paper",
            x0=wall_thickness, x1=wall_thickness + delta_c, y0=0.0, y1=1.0,
            fillcolor="lightgray", opacity=0.7, line_width=0, layer="below"
        )

        # ---- freestream lines ----
        self.fig.add_trace(go.Scatter(
            x=[x_min, -delta_g], y=[T_gas_edge, T_gas_edge],
            mode="lines",
            line=dict(color="crimson", width=2, dash="dash"),
            name="Gas recovery temperature", showlegend=False
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
        y_min = float(np.nanmin([np.nanmin(T_all), T_gas_edge, T_c]))
        y_max = float(np.nanmax([np.nanmax(T_all), T_gas_edge, T_c]))
        y_mid = 0.5 * (y_min + y_max)

        self.fig.add_annotation(x=x_min, y=y_mid, xref="x", yref="y",
                                text="Gas recovery edge", showarrow=False,
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

class PlotHeatTransferCoefficient(PlotBase):
    """
    Plot returned regenerative-cooling heat transfer coefficients.

    Expects one or more ``RegenResult`` objects.

    ``h_hot`` is the effective hot-side temperature-based coefficient and
    ``h_cold`` is the coolant-side coefficient, both in W/m²/K.
    """

    def __init__(self, *regen_results, hot: bool = True, cold: bool = True):
        super().__init__(go.Figure())
        self.template("plotly_white")

        for result in regen_results:
            name = result.circuit_name

            if hot:
                x = _result_grid(result, "x_heat_flux")
                y_hot = np.asarray(result.h_hot)
                self.fig.add_trace(go.Scatter(
                    x=x, y=y_hot, mode="lines",
                    name=f"{name} — Hot side",
                    showlegend=True,
                ))

            if cold:
                x = _result_grid(result, "x_coolant")
                y_cold = np.asarray(result.h_cold)
                self.fig.add_trace(go.Scatter(
                    x=x, y=y_cold, mode="lines",
                    name=f"{name} — Coolant side",
                    showlegend=True,
                    line=dict(dash="dash"),
                ))

        self.fig.update_layout(
            title="Heat Transfer Coefficient",
            xaxis=dict(title="Axial Position (m)"),
            yaxis=dict(title="h (W/m²/K)"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )


#: Fixed categorical colour per resistance layer, assigned in stack order and
#: never cycled. Slots 2/7/1 of the same reference palette as PHASE_COLORS;
#: validated all-pairs on a light surface (worst adjacent deutan ΔE 13.0,
#: worst normal-vision ΔE 16.3, all three above 3:1 against the surface).
#: The hues are semantic -- warm for the gas film, metal for the wall, cool
#: for the coolant -- so the stack reads without consulting the legend.
RESISTANCE_COLORS = {
    "hot": "#eb6834",    # orange
    "wall": "#4a3aa7",   # violet
    "cold": "#2a78d6",   # blue
}

RESISTANCE_LABELS = {
    "hot": "Hot-gas film",
    "wall": "Wall conduction",
    "cold": "Coolant film",
}


def _interp_sorted(x_target, x_source, y_source) -> np.ndarray:
    """Interpolate onto ``x_target``, tolerating a descending source grid.

    Counter-flow circuits march from the nozzle end, so their native grids run
    the other way and ``np.interp`` would silently return garbage.
    """
    x_source = np.asarray(x_source, dtype=float)
    y_source = np.asarray(y_source, dtype=float)
    order = np.argsort(x_source)
    return np.interp(
        np.asarray(x_target, dtype=float), x_source[order], y_source[order]
    )


class PlotThermalResistance(PlotBase):
    """Stack the hot-side, wall, and coolant-side thermal resistances.

    Every resistance is a temperature drop over the *same* hot-gas-area heat
    flux, so the three sum exactly to the total gas-to-coolant resistance and
    can be compared directly:

    ``R_hot = (T_aw - T_hw) / q''``,
    ``R_wall = (T_hw - T_cw) / q''``,
    ``R_cold = (T_cw - T_coolant) / q''``.

    That referral matters on the coolant side: ``h_cold`` is a bare convective
    coefficient on the wetted channel perimeter, so it is not comparable with
    ``h_hot`` until the channel-to-chamber area ratio is folded in, which this
    form does implicitly.

    Reading the shares tells you what limits the design -- a gas-film-dominated
    stack means a better coolant-side correlation buys nothing.

    Parameters
    ----------
    *regen_results : RegenResult
        Solved circuits. One at a time reads best: stacks from several circuits
        are drawn in separate stack groups and will overlay each other.
    mode : {'share', 'absolute'}, optional
        Percentage of the total resistance, or absolute m²·K/W.
    """

    def __init__(
        self,
        *regen_results,
        mode: str = "share",
        template: str = "plotly_white",
    ):
        if mode not in ("share", "absolute"):
            raise ValueError("mode must be 'share' or 'absolute'")

        super().__init__(go.Figure())
        self.template(template)

        for result in regen_results:
            x = _result_grid(result, "x_wall")
            T = np.asarray(result.T, dtype=float)
            if T.shape[1] < 3:
                raise ValueError(
                    "result.T needs at least a coolant, cold-wall and hot-wall "
                    "column to split the resistance stack"
                )

            # Column layout matches PlotWallTemperature: bulk coolant, cold
            # wall, optional interfaces, hot wall.
            T_cool, T_cw, T_hw = T[:, 0], T[:, 1], T[:, -1]

            x_flux = _result_grid(result, "x_heat_flux")
            q = _interp_sorted(x, x_flux, result.dQ_dA)
            T_aw = _interp_sorted(x, x_flux, result.T_aw_hot)

            # A vanishing flux makes every resistance blow up; the split is
            # meaningless there rather than merely large.
            q_ref = np.nanmax(np.abs(q)) if np.isfinite(q).any() else 0.0
            with np.errstate(divide="ignore", invalid="ignore"):
                usable = np.abs(q) > 1.0e-6 * q_ref
                q_safe = np.where(usable, q, np.nan)
                layers = {
                    "hot": (T_aw - T_hw) / q_safe,
                    "wall": (T_hw - T_cw) / q_safe,
                    "cold": (T_cw - T_cool) / q_safe,
                }

                if mode == "share":
                    total = sum(layers.values())
                    layers = {
                        key: 100.0 * value / total
                        for key, value in layers.items()
                    }

            for key, values in layers.items():
                self.fig.add_trace(go.Scatter(
                    x=x,
                    y=values,
                    mode="lines",
                    name=RESISTANCE_LABELS[key],
                    legendgroup=result.circuit_name,
                    legendgrouptitle=dict(text=result.circuit_name),
                    showlegend=True,
                    stackgroup=result.circuit_name,
                    fillcolor=RESISTANCE_COLORS[key],
                    # 2 px surface-coloured separator between stacked fills.
                    line=dict(color="#ffffff", width=2),
                    hovertemplate=(
                        f"{RESISTANCE_LABELS[key]}<br>"
                        "x = %{x:.4g} m<br>"
                        + (
                            "%{y:.1f}%% of total"
                            if mode == "share"
                            else "%{y:.3e} m²·K/W"
                        )
                        + "<extra></extra>"
                    ),
                ))

        if mode == "share":
            title = "Thermal Resistance Share"
            y_axis = dict(title="Share of total resistance (%)", range=[0, 100])
        else:
            title = "Thermal Resistance"
            y_axis = dict(title="Resistance (m²·K/W)")

        self.fig.update_layout(
            title=title,
            xaxis=dict(title="Axial Position (m)"),
            yaxis=y_axis,
            # Reversed so the legend reads in the same order as the stack.
            legend=dict(
                title=None, traceorder="reversed", groupclick="toggleitem"
            ),
            hovermode="x unified",
            margin=dict(l=60, r=20, t=60, b=60),
        )
