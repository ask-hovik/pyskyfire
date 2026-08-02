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


class PlotCoolantTemperature(PlotBase):
    """Plot bulk static coolant temperature from one or more ``RegenResult`` objects."""
    def __init__(self, *regen_results):
        super().__init__(go.Figure())
        self.template("plotly_white")

        for result in regen_results:
            x = np.asarray(result.x)
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
            x = np.asarray(result.x)
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
            ),
            margin=dict(l=60, r=20, t=60, b=60),
        )

class PlotWallTemperature(PlotBase):
    """
    Build a wall-temperature plot from one or more ``RegenResult`` objects.

    ``T`` must have shape ``(len(x), n_layers)``. The circuit name is used as the
    legend-group title.
    """

    def __init__(
        self,
        *regen_results,
        plot_hot: bool = True,
        plot_interfaces: bool = False,
        plot_coolant_wall: bool = False,
        template: str = "plotly_white",
    ):
        super().__init__(go.Figure())
        self.template(template)

        for result in regen_results:
            x = np.asarray(result.x)
            T = np.asarray(result.T)
            circuit_name = result.circuit_name

            cols = []

            if plot_coolant_wall and T.shape[1] > 1:
                cols.append(1)

            if plot_interfaces and T.shape[1] > 2:
                cols.extend(range(2, T.shape[1] - 1))

            if plot_hot:
                cols.append(T.shape[1] - 1)

            for col in cols:
                is_cold_wall = col == 1 and plot_coolant_wall
                is_hot_wall = col == T.shape[1] - 1 and plot_hot
                is_interface = plot_interfaces and 2 <= col <= T.shape[1] - 2

                if is_cold_wall:
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

        self.fig.update_layout(
            title="Wall Temperature",
            xaxis=dict(title="Axial Position (m)"),
            yaxis=dict(title="Temperature (K)"),
            legend=dict(
                title=None,
                tracegroupgap=10,
            ),
            margin=dict(l=60, r=20, t=60, b=60),
        )

class PlotHeatFlux(PlotBase):
    """Plot wall heat flux from one or more ``RegenResult`` objects."""
    def __init__(self, *regen_results):
        super().__init__(go.Figure())
        self.template("plotly_white")

        for result in regen_results:
            x = np.asarray(result.x)
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
            x = np.asarray(result.x)
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
            x = np.asarray(result.x, dtype=float)
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
            spans = [np.asarray(r.x, dtype=float) for r in regen_results]
            xaxis["range"] = [min(s.min() for s in spans),
                              max(s.max() for s in spans)]

        self.fig.update_layout(
            title="Coolant Vapour Quality",
            # Deliberately not "x": that is the axial coordinate on this figure.
            xaxis=xaxis,
            yaxis=dict(title="Vapour quality (-)", range=[-0.05, 1.05]),
            legend=dict(title=None, tracegroupgap=10),
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
            x = np.asarray(result.x, dtype=float)
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
            x = np.asarray(result.x, dtype=float)
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
            legend=dict(title=None, tracegroupgap=10),
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
        ``dQ_dA``, and ``p_static``
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
        T_inf = float(thrust_chamber.combustion_transport.get_T(x0))
        T_hw  = float(result.T[i, -1])    # hot-wall
        T_cw  = float(result.T[i,  1])    # coolant-side wall
        T_c   = float(result.T_static[i]) # bulk coolant
        qpp   = float(result.dQ_dA[i])

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
        Ts_rev   = result.T[i, 1:]       # coolant-side → hot-side
        Ts_wall  = Ts_rev[::-1]              # hot-side → coolant-side (interfaces included)
        walls    = thrust_chamber.cooling_circuits[circuit_index].walls
        thicknesses = [float(w.thickness(x0)) for w in walls]
        y_w = np.insert(np.cumsum(thicknesses), 0, 0.0) if thicknesses else np.array([0.0])
        wall_thickness = float(y_w[-1])

        # ---- 5) coolant BL (1/7th power law) ----
        p_static = float(result.p_static[i])
        T_film   = 0.5 * (T_c + T_cw)
        coolant  = thrust_chamber.cooling_circuits[circuit_index].coolant_transport
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
            x = np.asarray(result.x)
            name = result.circuit_name

            if hot:
                y_hot = np.asarray(result.h_hot)
                self.fig.add_trace(go.Scatter(
                    x=x, y=y_hot, mode="lines",
                    name=f"{name} — Hot side",
                    showlegend=True,
                ))

            if cold:
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