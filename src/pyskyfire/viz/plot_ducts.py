"""Pressure along a solved duct, against arc length."""

import numpy as np
import plotly.graph_objects as go

from .core import PlotBase

__all__ = ["PlotPipePressure", "PlotPipePressureLoss"]


def _labels(results, labels):
    if labels is None:
        return [f"Pipe {i + 1}" for i in range(len(results))] if len(results) > 1 else ["Pipe"]
    labels = list(labels)
    if len(labels) != len(results):
        raise ValueError("labels must have one entry per result")
    return labels


class PlotPipePressure(PlotBase):
    """Static and/or stagnation pressure against arc length for solved pipes.

    Each input is a :class:`~pyskyfire.ducts.DuctResult`, or any result with
    the same ``s``, ``p_stagnation`` and ``p_static`` station arrays. Markers
    sit on the :class:`~pyskyfire.ducts.PathStation` positions the solver
    evaluated; between them the line is only a straight join, so a lumped
    :class:`~pyskyfire.ducts.RoundedBend` shows as one step across its arc.
    """

    def __init__(self, *results, labels=None, static=True, stagnation=True,
                 markers=True, template="plotly_white"):
        super().__init__(go.Figure())
        self.template(template)
        mode = "lines+markers" if markers else "lines"
        for result, label in zip(results, _labels(results, labels)):
            s = np.asarray(result.s, dtype=float)
            if stagnation:
                self.fig.add_trace(go.Scatter(
                    x=s, y=np.asarray(result.p_stagnation), mode=mode,
                    name="Stagnation pressure", legendgroup=label,
                    legendgrouptitle=dict(text=label), line=dict(dash="dash"),
                ))
            if static:
                self.fig.add_trace(go.Scatter(
                    x=s, y=np.asarray(result.p_static), mode=mode,
                    name="Static pressure", legendgroup=label,
                    legendgrouptitle=dict(text=label), line=dict(dash="solid"),
                ))
        self.layout(
            title="Pipe pressure",
            xaxis=dict(title="Arc length, s (m)"),
            yaxis=dict(title="Pressure (Pa)"),
            legend=dict(title=None, tracegroupgap=10, groupclick="toggleitem"),
            margin=dict(l=60, r=20, t=60, b=60),
        )


class PlotPipePressureLoss(PlotBase):
    """Cumulative irreversible pressure loss against arc length for solved pipes.

    The loss at each station is the inlet stagnation pressure less the
    stagnation pressure there, so every curve starts at zero and ends at the
    result's ``total_pressure_loss``. Inputs are as for
    :class:`PlotPipePressure`.
    """

    def __init__(self, *results, labels=None, markers=True, template="plotly_white"):
        super().__init__(go.Figure())
        self.template(template)
        mode = "lines+markers" if markers else "lines"
        for result, label in zip(results, _labels(results, labels)):
            p = np.asarray(result.p_stagnation, dtype=float)
            self.fig.add_trace(go.Scatter(
                x=np.asarray(result.s, dtype=float), y=p[0] - p, mode=mode, name=label,
            ))
        self.layout(
            title="Cumulative pressure loss",
            xaxis=dict(title="Arc length, s (m)"),
            yaxis=dict(title="Stagnation pressure loss (Pa)"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )
