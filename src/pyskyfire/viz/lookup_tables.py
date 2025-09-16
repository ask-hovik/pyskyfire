# src/pyskyfire/viz/theta_vs_epsilon.py
import os, json, numpy as np, math
import plotly.graph_objects as go
from .core import PlotBase
import re
from pyskyfire.regen import f_darcy
from typing import Iterable, Optional, Sequence


class PlotThetaVsEpsilon(PlotBase):
    def __init__(self, template: str = "plotly_white"):
        super().__init__(go.Figure())
        self.template(template)

        # --- load data relative to this file ---
        base = os.path.dirname(os.path.abspath(__file__))
        data_dir = os.path.join(base, "..", "regen", "data")
        with open(os.path.join(data_dir, "theta_e.json"), encoding="utf-8") as fe:
            data_e = json.load(fe)
        with open(os.path.join(data_dir, "theta_n.json"), encoding="utf-8") as fn:
            data_n = json.load(fn)

        # --- muted color ramps (pick from qualitative lists) ---
        buGn   = ["#EDF8FB","#CCECE6","#99D8C9","#66C2A4","#41AE76","#238B45","#006D2C","#00441B"]
        ylOrRd = ["#FFFFCC","#FFEDA0","#FED976","#FEB24C","#FD8D3C","#FC4E2A","#E31A1C","#BD0026","#800026"]
        def pick(scale, lo=0.35, hi=0.75, n=5):
            idxs = np.linspace(int(lo*(len(scale)-1)), int(hi*(len(scale)-1)), n).round().astype(int)
            return [scale[i] for i in idxs]
        cols_e = dict(zip([60,70,80,90,100], pick(buGn, 0.45, 0.75, 5)))
        cols_n = dict(zip([60,70,80,90,100], pick(ylOrRd, 0.35, 0.75, 5)))

        # per-percent vertical nudges (points) for endpoint labels
        y_tweak_e = {60:0, 70:0, 80:3, 90:0, 100:-2}
        y_tweak_n = {60:0, 70:0, 80:4, 90:0, 100:-6}

        def pct_from_key(k: str) -> int:
            try: return int(k.split("_")[2])
            except: return 0

        all_eps = []

        # ---- θₑ curves ----
        for k, c in data_e.items():
            eps = np.asarray(c["epsilon"], float); th = np.asarray(c["theta"], float)
            pct = pct_from_key(k); col = cols_e.get(pct, "#238B45")
            all_eps.append(eps)

            # θ_e traces
            self.fig.add_trace(go.Scatter(
                x=eps, y=th, mode="lines",
                line=dict(color=col, width=2),
                name=f"{pct}%",                 # <- only the percent
                legendgroup="theta_e",          # <- family tag
                meta="e",                       # <- alt tag (redundant but handy)
                showlegend=False
            ))

            # centered θₑ label on 60%
            if pct == 60:
                mid = len(eps)//2
                self.fig.add_annotation(
                    x=float(eps[mid]), y=float(th[mid]),
                    xref="x", yref="y",
                    text="θₑ", xanchor="center", yanchor="bottom",
                    xshift=0, yshift=10, showarrow=False,
                    font=dict(size=13, color="black")
                )

        # ---- θₙ curves ----
        for k, c in data_n.items():
            eps = np.asarray(c["epsilon"], float); th = np.asarray(c["theta"], float)
            pct = pct_from_key(k); col = cols_n.get(pct, "#FD8D3C")
            all_eps.append(eps)

            # θ_n traces
            self.fig.add_trace(go.Scatter(
                x=eps, y=th, mode="lines",
                line=dict(color=col, width=2),
                name=f"{pct}%",
                legendgroup="theta_n",
                meta="n",
                showlegend=False
            ))

            # centered θₙ label on 60%
            if pct == 60:
                mid = len(eps)//2
                self.fig.add_annotation(
                    x=float(eps[mid]), y=float(th[mid]),
                    xref="x", yref="y",
                    text="θₙ", xanchor="center", yanchor="bottom",
                    xshift=0, yshift=10, showarrow=False,
                    font=dict(size=13, color="black")
                )


            for tr in self.fig.data:
                xs = np.asarray(tr.x, float); ys = np.asarray(tr.y, float)
                m = np.isfinite(xs) & np.isfinite(ys)
                if not m.any():
                    continue
                i = np.where(m)[0][-1]

                # family: prefer legendgroup, fall back to meta
                lg = (getattr(tr, "legendgroup", "") or "").lower()
                fam = "e" if lg == "theta_e" else "n" if lg == "theta_n" else (tr.meta if tr.meta in ("e","n") else None)

                # percent from name like "80%"
                txt = tr.name or ""
                m_pct = re.search(r"(\d+)%", txt)
                pct = int(m_pct.group(1)) if m_pct else None

                yshift = (y_tweak_e.get(pct, 0) if fam == "e"
                        else y_tweak_n.get(pct, 0) if fam == "n"
                        else 0)
                color = (getattr(getattr(tr, "line", None), "color", None) or "black")

                self.fig.add_annotation(
                    x=float(xs[i]), y=float(ys[i]),
                    xref="x", yref="y",
                    text=txt,                       # e.g., "80%"
                    xanchor="left", yanchor="middle",
                    xshift=6, yshift=yshift,
                    showarrow=False,
                    font=dict(size=11, color=color),
                )

        self.fig.update_layout(
            title=r"$\theta \ \text{vs}\ \varepsilon$",
            xaxis=dict(title=r"$\varepsilon$"),
            yaxis=dict(title=r"$\theta\ (^{\circ})$")
        )

# TODO: the moody diagram log with annotations doesn't work well. So have to redo it. But the bones are there. 
class PlotMoodyDiagram(PlotBase):
    """
    Moody-style plot of Darcy friction factor vs Reynolds number, using
    the correlations implemented in pyskyfire.physics.f_darcy.

    Each curve corresponds to a chosen relative roughness ε/D.
    """

    def __init__(
        self,
        rel_rough_list: Optional[Sequence[float]] = None,  # ε/D values
        Re_min: float = 7e2,
        Re_max: float = 1e8,
        n_pts: int = 400,
        template: str = "plotly_white",
    ):
        super().__init__(go.Figure())
        self.template(template)

        # Defaults (same spirit as your matplotlib version)
        if rel_rough_list is None:
            rel_rough_list = [
                     3e-5,
                1e-4, 3e-4,
                1e-3, 3e-3,
                1e-2, 3e-2,
            ]

        # Reynolds sample (log-spaced)
        Re_vals = np.logspace(np.log10(Re_min), np.log10(Re_max), n_pts)

        # Use Dh=1 so absolute roughness = ε (since ε/D * Dh = ε)
        Dh = 1.0
        x_dummy = 0.0  # physics.f_darcy signature requires x; not used for constants

        # Plot a curve per ε/D
        for eps_over_D in rel_rough_list:
            if eps_over_D == 0:
                roughness = None
            else:
                # constant absolute roughness ε = (ε/D) * Dh
                eps_abs = eps_over_D * Dh
                roughness = (lambda _x, e=eps_abs: e)

            f_vals = np.array([f_darcy(ReDh=Re, Dh=Dh, x=x_dummy, roughness=roughness) for Re in Re_vals])

            label = f"{eps_over_D:.6f}".rstrip("0").rstrip(".")
            self.fig.add_trace(go.Scatter(
                x=Re_vals, y=f_vals,
                mode="lines",
                name=label,              # legend/hover name (ε/D)
                legendgroup=label,
                showlegend=False,        # we'll label on-canvas instead
            ))

        # Axis setup: log Re; friction factor axis commonly shown log as well
        # (you can switch yaxis type to "linear" if you prefer)
        self.fig.update_layout(
            title="Darcy friction factor vs Reynolds number",
            xaxis=dict(title="Re", type="log", showgrid=True),
            yaxis=dict(title="f",  type="log", showgrid=True),
            margin=dict(l=70, r=70, t=50, b=60),
            legend=dict(title=None),
        )

        # Optional: set nice y ticks at 0.01 steps (works on log too, like your MPL)
        # Compute a reasonable range from data first
        all_f = np.concatenate([np.asarray(tr.y, float) for tr in self.fig.data])
        fmin = float(np.nanmin(all_f))
        fmax = float(np.nanmax(all_f))
        # keep existing auto-range but provide explicit tickvals in 0.01 increments
        # between the current y-limits intersection
        lo = max(0.005, 10 ** self.fig.layout.yaxis.range[0] if self.fig.layout.yaxis.range else fmin)
        hi = 10 ** self.fig.layout.yaxis.range[1] if self.fig.layout.yaxis.range else fmax
        ticks = np.arange(np.ceil(lo * 100) / 100, np.floor(hi * 100) / 100 + 1e-9, 0.01)
        if ticks.size > 0 and ticks.size <= 60:  # avoid too many tick labels
            self.fig.update_yaxes(tickvals=ticks, tickformat=".2f")

        # Add endpoint labels near the last valid point of each curve
        # and pad Re-range on the right so labels aren’t clipped on a log axis
        xs_all = []
        for tr in self.fig.data:
            x = np.asarray(tr.x, float)
            y = np.asarray(tr.y, float)
            m = np.isfinite(x) & np.isfinite(y)
            if not m.any():
                continue
            i = np.where(m)[0][-1]
            xs_all.append(x)

            self.fig.add_annotation(
                x=float(x[i]), y=float(y[i]),
                xref="x", yref="y",
                text=tr.name,  # ε/D text
                xanchor="left", yanchor="bottom",
                xshift=6, yshift=0,
                showarrow=False,
                font=dict(size=10, color="black"),
            )

        # Right padding (~15%) in log units so labels have breathing room
        if xs_all:
            xcat = np.concatenate(xs_all)
            xmin = float(np.nanmin(xcat))
            xmax = float(np.nanmax(xcat))
            lo = np.log10(max(xmin, 1e-12))
            hi = np.log10(max(xmax, xmin * 1.0001))
            self.fig.update_xaxes(range=[lo, hi + np.log10(1.15)])

