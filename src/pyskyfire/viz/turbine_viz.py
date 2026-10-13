"""Plotly figures of an axial turbine: triangles, passage, cascade and spanwise angles.

Each function takes one of the turbine objects and returns a ``go.Figure``,
ready for :meth:`Report.add_figure` or ``fig.show()``. Lengths are shown in
millimetres and angles in degrees. Pitchwise distances and angles are positive
in the direction the rotor turns, as in :mod:`pyskyfire.turbine`.
"""

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

__all__ = ["plot_velocity_triangles", "plot_turbine_meridional", "plot_turbine_cascade",
           "plot_turbine_spanwise"]

_STATOR, _ROTOR = "slategray", "tomato"


def _row_color(row):
    return _ROTOR if row.is_rotor else _STATOR


def _arrow(fig, start, end, color, text, col, dash=None):
    fig.add_trace(go.Scatter(x=[start[0], end[0]], y=[start[1], end[1]], mode="lines",
                             line=dict(color=color, width=2.5, dash=dash), name=text,
                             showlegend=False, hoverinfo="skip"), row=1, col=col)
    axes = "" if col == 1 else str(col)
    fig.add_annotation(x=end[0], y=end[1], ax=start[0], ay=start[1], xref=f"x{axes}", yref=f"y{axes}",
                       axref=f"x{axes}", ayref=f"y{axes}", showarrow=True, arrowhead=3,
                       arrowsize=1.2, arrowwidth=2, arrowcolor=color, text="")
    mid = 0.5*(np.asarray(start) + np.asarray(end))
    fig.add_annotation(x=mid[0], y=mid[1], xref=f"x{axes}", yref=f"y{axes}", text=text,
                       showarrow=False, font=dict(color=color), xshift=12)


def plot_velocity_triangles(turbine):
    """Each stage's rotor inlet and exit triangles, drawn from a common apex.

    Velocities point down the page in the flow direction, with swirl across
    it, positive in the rotation direction. At each station the absolute
    velocity ``c`` is the relative velocity ``w`` plus the blade speed ``u``.
    """
    stages = turbine.stages
    fig = make_subplots(rows=1, cols=len(stages),
                        subplot_titles=[f"Stage {k + 1}" for k in range(len(stages))])
    for k, stage in enumerate(stages, start=1):
        for station, label, dash in ((stage.middle, "2", None), (stage.exit, "3", "dot")):
            c = (station.c_u, -station.c_x)
            w = (station.w_u, -station.c_x)
            _arrow(fig, (0, 0), c, "black", f"c{label}", k, dash)
            _arrow(fig, (0, 0), w, _ROTOR, f"w{label}", k, dash)
            _arrow(fig, w, c, "royalblue", "u", k, dash)
        fig.update_xaxes(title_text="Swirl, rotation direction [m/s]", zeroline=True, row=1, col=k)
        fig.update_yaxes(title_text="Axial, flow downwards [m/s]", scaleanchor=f"x{'' if k == 1 else k}",
                         row=1, col=k)
    fig.update_layout(title="Velocity triangles at mean radius: rotor inlet (2, solid) and exit (3, dotted)",
                      template="plotly_white")
    return fig


def plot_turbine_meridional(geometry):
    """The passage in (x, r), with rotating wall stretches, gaps and blade edges.

    ``geometry`` is a :class:`~pyskyfire.turbine.TurbineGeometry`. Thick red
    stretches of wall turn with the shaft; the blade outlines are each row's
    leading and trailing edges and its tip or hub end.
    """
    fig = go.Figure()
    p = geometry.passage
    for wall in ("hub", "shroud"):
        curve = getattr(p, wall)
        fig.add_trace(go.Scatter(x=curve[:, 0]*1e3, y=curve[:, 1]*1e3, name=wall.capitalize(),
                                 line=dict(color="black", width=1.5)))
        turning = p.rotating(wall)
        if turning.any():
            fig.add_trace(go.Scatter(x=np.where(turning, curve[:, 0], np.nan)*1e3,
                                     y=np.where(turning, curve[:, 1], np.nan)*1e3,
                                     name=f"{wall.capitalize()}, rotating",
                                     line=dict(color=_ROTOR, width=6)))
    for rg in geometry.rows:
        side = rg.blade.suction
        outline = np.vstack((side[:, 0], side[-1, :], side[::-1, -1], side[0, ::-1]))
        fig.add_trace(go.Scatter(x=outline[:, 0]*1e3, y=outline[:, 1]*1e3, name=rg.row.name,
                                 fill="toself", line=dict(color=_row_color(rg.row), width=1),
                                 opacity=0.6))
    fig.update_layout(title="Meridional passage", xaxis_title="x [mm]", yaxis_title="r [mm]",
                      yaxis_scaleanchor="x", template="plotly_white")
    return fig


def plot_turbine_cascade(geometry, span_fraction=0.5, passages=3):
    """Every row at one span, unrolled flat: axial against pitchwise distance.

    Suction sides are solid and pressure sides dashed. Each row shows
    ``passages`` neighbouring blades at its own pitch.
    """
    fig = go.Figure()
    sign = geometry.rotation_sign
    for rg in geometry.rows:
        k = int(np.argmin(abs(rg.span_fractions - span_fraction)))
        r, pitch = rg.section_radius[k], rg.pitch[k]
        for name, side, dash in (("suction", rg.blade.suction, "solid"),
                                 ("pressure", rg.blade.pressure, "dash")):
            pts = side[k]
            for i in range(passages):
                fig.add_trace(go.Scatter(
                    x=pts[:, 0]*1e3, y=(sign*pts[:, 2]*r + i*pitch)*1e3,
                    line=dict(color=_row_color(rg.row), dash=dash, width=1.5),
                    name=f"{rg.row.name} {name}", legendgroup=f"{rg.row.name}{name}",
                    showlegend=i == 0))
    span = geometry.rows[0].span_fractions
    shown = span[int(np.argmin(abs(span - span_fraction)))]
    fig.update_layout(title=f"Blade sections at {shown:.0%} span, unrolled",
                      xaxis_title="x [mm]", yaxis_title="r·θ, rotation direction [mm]",
                      yaxis_scaleanchor="x", template="plotly_white")
    return fig


def plot_turbine_spanwise(geometry):
    """Each row's metal angles, measured gauging angle, throat and pitch, hub to tip.

    The gauging angle ``acos(throat/pitch)`` is measured on the built blades;
    beside the outlet metal angle it shows how far the section turns the flow
    from what it was drawn for. Line styles are the same in every row: see
    the subplot titles.
    """
    fig = make_subplots(rows=1, cols=2, subplot_titles=(
        "Angles [deg]: inlet metal dotted, outlet metal solid, gauging dashed",
        "Throat solid, pitch dotted [mm]"))
    for rg in geometry.rows:
        color, name = _row_color(rg.row), rg.row.name
        inlet = [np.degrees(rg.row.section_angles(f, r)[0])
                 for f, r in zip(rg.span_fractions, rg.section_radius)]
        traces = (((inlet, "dot"), (np.degrees(rg.outlet_metal_angle), "solid"),
                   (np.degrees(rg.gauging_angle), "dash")),
                  ((rg.throat*1e3, "solid"), (rg.pitch*1e3, "dot")))
        first = True
        for col, lines in enumerate(traces, start=1):
            for x, dash in lines:
                fig.add_trace(go.Scatter(x=x, y=rg.span_fractions, name=name, legendgroup=name,
                                         showlegend=first, mode="lines",
                                         line=dict(color=color, dash=dash)), row=1, col=col)
                first = False
    fig.update_yaxes(title_text="Span fraction", row=1, col=1)
    fig.update_layout(title="Hub to tip", template="plotly_white")
    return fig
