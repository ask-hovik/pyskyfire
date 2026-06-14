
"""
pyskyfire/viz/plot_film_cooling.py

PlotBase subclasses for visualising contour-based Grisson film cooling results.

These plots are built around the current result containers returned by:
    GrissonFilmCoolingModel.solve(...)

Returned data
-------------
LiquidFilmResults:
    x
    Gamma
    m_vap
    h_conv
    Q_conv
    Q_rad
    T_film
    x_dryout
    Gamma_dryout

GaseousFilmResults:
    x
    Mbl
    T_aw
    T_w
    h_conv
    Q_rad
"""
from __future__ import annotations

import numpy as np
import plotly.graph_objects as go

from .core import PlotBase


def _maybe_array(lst) -> np.ndarray:
    return np.asarray(lst, dtype=float)


class PlotFilmWallTemperature(PlotBase):
    """Combined temperature plot for liquid and gaseous film phases."""

    def __init__(
        self,
        liquid,
        gaseous,
        T_g_ref: float | None = None,
        T_sat_ref: float | None = None,
        template: str = "plotly_white",
    ):
        super().__init__(go.Figure())
        self.template(template)

        if liquid.x:
            x_liq = _maybe_array(liquid.x)
            T_film = _maybe_array(liquid.T_film)
            self.fig.add_trace(go.Scatter(
                x=x_liq, y=T_film,
                mode="lines",
                name="T_film",
                line=dict(width=2),
            ))

        if gaseous.x:
            x_gas = _maybe_array(gaseous.x)
            T_aw = _maybe_array(gaseous.T_aw)
            T_w = _maybe_array(gaseous.T_w)

            self.fig.add_trace(go.Scatter(
                x=x_gas, y=T_aw,
                mode="lines",
                name="T_aw",
                line=dict(width=2),
            ))
            self.fig.add_trace(go.Scatter(
                x=x_gas, y=T_w,
                mode="lines",
                name="T_w",
                line=dict(width=2, dash="dash"),
            ))

        if liquid.x_dryout is not None:
            self.fig.add_vline(
                x=liquid.x_dryout,
                line=dict(dash="dot", width=1.5),
                annotation_text=f"Dryout x={liquid.x_dryout:.4f} m",
                annotation_position="top right",
            )

        if T_g_ref is not None:
            self.fig.add_hline(
                y=T_g_ref,
                line=dict(dash="dot", width=1),
                annotation_text=f"T_g = {T_g_ref:.0f} K",
                annotation_position="bottom right",
            )

        if T_sat_ref is not None:
            self.fig.add_hline(
                y=T_sat_ref,
                line=dict(dash="dot", width=1),
                annotation_text=f"T_sat = {T_sat_ref:.0f} K",
                annotation_position="top right",
            )

        self.fig.update_layout(
            title="Film Cooling — Temperature Profile",
            xaxis=dict(title="Axial position, x (m)"),
            yaxis=dict(title="Temperature (K)"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )


class PlotFilmDryout(PlotBase):
    """Liquid-film loading Γ and dryout point."""

    def __init__(self, liquid, template: str = "plotly_white"):
        super().__init__(go.Figure())
        self.template(template)

        if liquid.x:
            x = _maybe_array(liquid.x)
            Gamma = _maybe_array(liquid.Gamma)

            self.fig.add_trace(go.Scatter(
                x=x, y=Gamma,
                mode="lines",
                name="Γ",
                fill="tozeroy",
                fillcolor="rgba(31,119,180,0.15)",
                line=dict(width=2),
            ))

        if liquid.x_dryout is not None:
            self.fig.add_vline(
                x=liquid.x_dryout,
                line=dict(dash="dot", width=2),
                annotation_text=f"Dryout x={liquid.x_dryout:.4f} m",
                annotation_position="top right",
            )
            self.fig.add_trace(go.Scatter(
                x=[liquid.x_dryout],
                y=[liquid.Gamma_dryout],
                mode="markers",
                name="Dryout point",
                marker=dict(size=10, symbol="x"),
            ))
        else:
            self.fig.add_annotation(
                text="Film survived to end of domain",
                xref="paper", yref="paper",
                x=0.98, y=0.95,
                showarrow=False,
                align="right",
            )

        self.fig.update_layout(
            title="Film Cooling — Film Loading and Dryout",
            xaxis=dict(title="Axial position, x (m)"),
            yaxis=dict(title="Γ (kg m⁻¹ s⁻¹)"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )


class PlotFilmBoundaryLayer(PlotBase):
    """Boundary-layer mass flux in the gaseous film phase."""

    def __init__(self, gaseous, template: str = "plotly_white"):
        super().__init__(go.Figure())
        self.template(template)

        if not gaseous.x:
            self.fig.add_annotation(
                text="No gaseous phase (film did not dry out)",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=14, color="gray"),
            )
        else:
            x = _maybe_array(gaseous.x)
            Mbl = _maybe_array(gaseous.Mbl)

            self.fig.add_trace(go.Scatter(
                x=x, y=Mbl,
                mode="lines",
                name="M_bl",
                line=dict(width=2),
            ))

        self.fig.update_layout(
            title="Film Cooling — Gaseous Boundary-Layer Mass Flux",
            xaxis=dict(title="Axial position, x (m)"),
            yaxis=dict(title="M_bl (kg m⁻² s⁻¹)"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )


class PlotFilmHeatFlux(PlotBase):
    """Convective, radiative, and total heat flux into the liquid film."""

    def __init__(self, liquid, template: str = "plotly_white"):
        super().__init__(go.Figure())
        self.template(template)

        if liquid.x:
            x = _maybe_array(liquid.x)
            Q_conv = _maybe_array(liquid.Q_conv)
            Q_rad = _maybe_array(liquid.Q_rad)
            Q_total = Q_conv + Q_rad

            self.fig.add_trace(go.Scatter(
                x=x, y=Q_conv,
                mode="lines",
                name="Q_conv",
                line=dict(width=2),
            ))
            self.fig.add_trace(go.Scatter(
                x=x, y=Q_rad,
                mode="lines",
                name="Q_rad",
                line=dict(width=2),
            ))
            self.fig.add_trace(go.Scatter(
                x=x, y=Q_total,
                mode="lines",
                name="Q_total",
                line=dict(width=2, dash="dash"),
            ))

        if liquid.x_dryout is not None:
            self.fig.add_vline(
                x=liquid.x_dryout,
                line=dict(dash="dot", width=1.5),
                annotation_text="Dryout",
                annotation_position="top right",
            )

        self.fig.update_layout(
            title="Film Cooling — Heat Flux into Liquid Film",
            xaxis=dict(title="Axial position, x (m)"),
            yaxis=dict(title="Heat flux (W m⁻²)"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )


class PlotFilmEvaporationRate(PlotBase):
    """Local liquid-film evaporation rate."""

    def __init__(self, liquid, template: str = "plotly_white"):
        super().__init__(go.Figure())
        self.template(template)

        if liquid.x:
            x = _maybe_array(liquid.x)
            m_vap = _maybe_array(liquid.m_vap)

            self.fig.add_trace(go.Scatter(
                x=x, y=m_vap,
                mode="lines",
                name="ṁ_vap",
                fill="tozeroy",
                fillcolor="rgba(140,86,75,0.15)",
                line=dict(width=2),
            ))

        if liquid.x_dryout is not None:
            self.fig.add_vline(
                x=liquid.x_dryout,
                line=dict(dash="dot", width=1.5),
                annotation_text="Dryout",
                annotation_position="top right",
            )

        self.fig.update_layout(
            title="Film Cooling — Evaporation Rate",
            xaxis=dict(title="Axial position, x (m)"),
            yaxis=dict(title="ṁ_vap (kg m⁻² s⁻¹)"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )


class PlotFilmHeatTransferCoefficient(PlotBase):
    """Heat-transfer coefficient in both liquid and gaseous film phases."""

    def __init__(self, liquid, gaseous, template: str = "plotly_white"):
        super().__init__(go.Figure())
        self.template(template)

        if liquid.x:
            self.fig.add_trace(go.Scatter(
                x=_maybe_array(liquid.x),
                y=_maybe_array(liquid.h_conv),
                mode="lines",
                name="h_conv (liquid)",
                line=dict(width=2),
            ))

        if gaseous.x:
            self.fig.add_trace(go.Scatter(
                x=_maybe_array(gaseous.x),
                y=_maybe_array(gaseous.h_conv),
                mode="lines",
                name="h_conv (gaseous)",
                line=dict(width=2, dash="dash"),
            ))

        if liquid.x_dryout is not None:
            self.fig.add_vline(
                x=liquid.x_dryout,
                line=dict(dash="dot", width=1.5),
                annotation_text="Dryout",
                annotation_position="top right",
            )

        self.fig.update_layout(
            title="Film Cooling — Convective Heat-Transfer Coefficient",
            xaxis=dict(title="Axial position, x (m)"),
            yaxis=dict(title="h_conv (W m⁻² K⁻¹)"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )


class PlotFilmRadiativeFraction(PlotBase):
    """Radiative contribution to total heat flux in the liquid-film phase."""

    def __init__(self, liquid, template: str = "plotly_white"):
        super().__init__(go.Figure())
        self.template(template)

        if liquid.x:
            x = _maybe_array(liquid.x)
            Q_conv = _maybe_array(liquid.Q_conv)
            Q_rad = _maybe_array(liquid.Q_rad)
            denom = np.maximum(Q_conv + Q_rad, 1e-12)
            frac = Q_rad / denom

            self.fig.add_trace(go.Scatter(
                x=x, y=frac,
                mode="lines",
                name="Q_rad / Q_total",
                fill="tozeroy",
                fillcolor="rgba(148,103,189,0.15)",
                line=dict(width=2),
            ))

        if liquid.x_dryout is not None:
            self.fig.add_vline(
                x=liquid.x_dryout,
                line=dict(dash="dot", width=1.5),
                annotation_text="Dryout",
                annotation_position="top right",
            )

        self.fig.update_layout(
            title="Film Cooling — Radiative Fraction of Liquid-Phase Heat Load",
            xaxis=dict(title="Axial position, x (m)"),
            yaxis=dict(title="Radiative fraction (-)"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )


class PlotGaseousFilmRadiation(PlotBase):
    """Radiative heat flux in the gaseous film phase."""

    def __init__(self, gaseous, template: str = "plotly_white"):
        super().__init__(go.Figure())
        self.template(template)

        if not gaseous.x:
            self.fig.add_annotation(
                text="No gaseous phase (film did not dry out)",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=14, color="gray"),
            )
        else:
            self.fig.add_trace(go.Scatter(
                x=_maybe_array(gaseous.x),
                y=_maybe_array(gaseous.Q_rad),
                mode="lines",
                name="Q_rad (gaseous phase)",
                line=dict(width=2),
            ))

        self.fig.update_layout(
            title="Film Cooling — Gaseous-Phase Radiative Heat Flux",
            xaxis=dict(title="Axial position, x (m)"),
            yaxis=dict(title="Q_rad (W m⁻²)"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )


class PlotGaseousFilmTemperatureGap(PlotBase):
    """Difference between actual and adiabatic wall temperature in gaseous phase."""

    def __init__(self, gaseous, template: str = "plotly_white"):
        super().__init__(go.Figure())
        self.template(template)

        if not gaseous.x:
            self.fig.add_annotation(
                text="No gaseous phase (film did not dry out)",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=14, color="gray"),
            )
        else:
            x = _maybe_array(gaseous.x)
            gap = _maybe_array(gaseous.T_w) - _maybe_array(gaseous.T_aw)

            self.fig.add_trace(go.Scatter(
                x=x, y=gap,
                mode="lines",
                name="T_w - T_aw",
                line=dict(width=2),
            ))

        self.fig.update_layout(
            title="Film Cooling — Gaseous-Phase Temperature Gap",
            xaxis=dict(title="Axial position, x (m)"),
            yaxis=dict(title="T_w - T_aw (K)"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )


class PlotFilmEffectiveWallHTC(PlotBase):
    """Effective wall-facing heat transfer coefficient for liquid and gaseous film phases."""

    def __init__(self, liquid, gaseous, template: str = "plotly_white"):
        super().__init__(go.Figure())
        self.template(template)

        if liquid.x and hasattr(liquid, "h_wall_eff"):
            self.fig.add_trace(go.Scatter(
                x=_maybe_array(liquid.x),
                y=_maybe_array(liquid.h_wall_eff),
                mode="lines",
                name="h_wall_eff (liquid)",
                line=dict(width=2),
            ))

        if gaseous.x and hasattr(gaseous, "h_wall_eff"):
            self.fig.add_trace(go.Scatter(
                x=_maybe_array(gaseous.x),
                y=_maybe_array(gaseous.h_wall_eff),
                mode="lines",
                name="h_wall_eff (gaseous)",
                line=dict(width=2, dash="dash"),
            ))

        if liquid.x_dryout is not None:
            self.fig.add_vline(
                x=liquid.x_dryout,
                line=dict(dash="dot", width=1.5),
                annotation_text="Dryout",
                annotation_position="top right",
            )

        self.fig.update_layout(
            title="Film Cooling — Effective Wall Heat-Transfer Coefficient",
            xaxis=dict(title="Axial position, x (m)"),
            yaxis=dict(title="h_wall_eff (W m⁻² K⁻¹)"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )