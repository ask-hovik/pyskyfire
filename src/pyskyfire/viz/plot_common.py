from __future__ import annotations

import numpy as np
import plotly.graph_objects as go

from .core import PlotBase


class PlotResidualHistory(PlotBase):
    def __init__(self, residuals, name="Residual", title="Residual Convergence History", template="plotly_white"):
        super().__init__(go.Figure()); self.template(template)
        r = np.asarray(residuals, float); r = np.where(r > 0, r, np.nan)
        x = np.arange(1, r.size + 1)
        self.fig.add_trace(go.Scatter(x=x, y=r, mode="lines+markers", name=name))
        self.layout(
            title=title or None,
            xaxis=dict(title="Iteration"),
            yaxis=dict(title="Residual (dimensionless)", type="log"),
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=60),
        )


class PlotStationProperty(PlotBase):
    def __init__(self,
                 station_dicts,
                 station_list,
                 property_name,
                 labels=None,
                 title=True,
                 ylabel=None,
                 template="plotly_white",
                 ylim=None):
        super().__init__(go.Figure()); self.template(template)

        # normalize inputs
        from collections.abc import Iterable
        if not isinstance(station_dicts, Iterable) or isinstance(station_dicts, dict):
            station_dicts = [station_dicts]
        n = len(station_dicts)

        # labels / legend
        if labels is False:
            show_legend = False
            labels_list = [f"Case {i+1}" for i in range(n)]
        elif labels is True or labels is None:
            show_legend = True
            labels_list = [f"Case {i+1}" for i in range(n)]
        else:
            if len(labels) != n: raise ValueError("len(labels) must match number of station_dicts")
            show_legend = True
            labels_list = labels

        def get_val(st, pname, sname):
            if isinstance(st, dict):
                if pname not in st: raise KeyError(f"Property '{pname}' missing for station '{sname}'")
                return st[pname]
            if hasattr(st, pname): return getattr(st, pname)
            raise TypeError(f"Unsupported station entry type {type(st)} for '{sname}'")

        x = np.arange(1, len(station_list)+1)
        ticktext = [s.replace("_"," ").title() for s in station_list]

        for d, lab in zip(station_dicts, labels_list):
            y = [get_val(d[name], property_name, name) for name in station_list]
            self.fig.add_trace(go.Scatter(x=x, y=y, mode="lines+markers",
                                          name=lab, showlegend=show_legend))

        ttl = (f"{property_name} vs Station" if title is True else (title or None))
        ya = dict(title=ylabel or property_name)
        if ylim is not None: ya["range"] = list(ylim)

        self.layout(
            title=ttl,
            xaxis=dict(title=None, tickmode="array", tickvals=x, ticktext=ticktext, tickangle=90),
            yaxis=ya,
            legend=dict(title=None),
            margin=dict(l=60, r=20, t=60, b=80),
        )


class PlotPTDiagram(PlotBase):
    def __init__(self,
                 station_dicts,
                 station_list,
                 fluid_name=None,
                 title=None,
                 sat_points=200,
                 labels=None,
                 template="plotly_white",
                 annotate_ha=None,
                 annotate_va=None,
                 scale="log"):
        super().__init__(go.Figure()); self.template(template)

        # normalize input
        from collections.abc import Iterable
        if not isinstance(station_dicts, Iterable) or isinstance(station_dicts, dict):
            station_dicts = [station_dicts]
        n_cases = len(station_dicts)

        # labels
        labels = labels or [f"Case {i+1}" for i in range(n_cases)]
        if len(labels) != n_cases:
            raise ValueError("len(labels) must equal number of station dictionaries.")

        # extractor
        def get_TP(obj, name):
            if hasattr(obj, "T") and hasattr(obj, "p"):
                return float(obj.T), float(obj.p)
            if isinstance(obj, dict):
                return float(obj["T"]), float(obj["p"])
            raise TypeError(f"Unsupported station entry type {type(obj)} for '{name}'.")

        # gather data
        all_T, all_P = [], []
        for d in station_dicts:
            Ts, Ps = [], []
            for name in station_list:
                if name not in d:
                    raise KeyError(f"Station '{name}' not found in station_dict.")
                t, p = get_TP(d[name], name)
                Ts.append(t); Ps.append(p)
            all_T.append(Ts); all_P.append(Ps)

        # traces
        for Ts, Ps, lab in zip(all_T, all_P, labels):
            self.fig.add_trace(go.Scatter(x=Ts, y=Ps, mode="lines+markers", name=lab))

        # label annotations (pixel offsets via ha/va)
        def _broadcast(arg, default):
            if arg is None or isinstance(arg, str): return [arg or default]*len(station_list)
            if len(arg) != len(station_list): raise ValueError("annotate_* must match station_list length.")
            return arg
        has = _broadcast(annotate_ha, "left")
        vas = _broadcast(annotate_va, "top")
        def _offset(ha, va, d=8):
            dx = {"left": d, "center": 0, "right": -d}.get(ha, 0)
            dy = {"bottom": d, "middle": 0, "center": 0, "top": -d}.get(va, 0)
            return dx, dy

        for Ts, Ps in zip(all_T, all_P):
            for name, T, P, ha, va in zip(station_list, Ts, Ps, has, vas):
                dx, dy = _offset(ha, va)
                self.fig.add_annotation(
                    x=T, y=P,
                    text=name.replace("_", " ").title(),
                    showarrow=False,
                    xref="x", yref="y",
                    xanchor=ha or "left", yanchor=va or "top",
                    ax=dx, ay=dy
                )

        # saturation curve + critical point (CoolProp optional)
        T_cr = P_cr = None
        if fluid_name:
            try:
                from CoolProp.CoolProp import PropsSI
                T_tr = float(PropsSI("Ttriple", fluid_name))
                T_cr = float(PropsSI("Tcrit",   fluid_name))
                P_cr = float(PropsSI("Pcrit",   fluid_name))
                T_sat = np.linspace(T_tr*1.01, T_cr*0.99, int(sat_points))
                P_sat = [float(PropsSI("P", "T", T, "Q", 0, fluid_name)) for T in T_sat]
                self.fig.add_trace(go.Scatter(
                    x=T_sat, y=P_sat, mode="lines",
                    name=f"{fluid_name} sat. line",
                    line=dict(dash="dash")
                ))
                self.fig.add_trace(go.Scatter(
                    x=[T_cr], y=[P_cr], mode="markers",
                    name="Critical pt", marker=dict(symbol="x", size=10)
                ))
            except Exception:
                # CoolProp unavailable or bad fluid name → skip saturation overlay
                pass

        # scales + zoom
        all_T_flat = [t for Ts in all_T for t in Ts]
        all_P_flat = [p for Ps in all_P for p in Ps]
        tmin, tmax = min(all_T_flat), max(all_T_flat)
        pmin, pmax = min(all_P_flat), max(all_P_flat)

        xaxis = dict(title="Temperature (K)")
        yaxis = dict(title="Pressure (Pa)")

        if scale == "log":
            fx, fy = 1.1, 1.4
            xmin = max(min(tmin, T_cr or tmin)/fx, 1e-6)
            xmax = (max(tmax, T_cr or tmax))*fx
            ymin = max(min(pmin, P_cr or pmin)/fy, 1e-6)
            ymax = (max(pmax, P_cr or pmax))*fy
            xaxis.update(type="log", range=[np.log10(xmin), np.log10(xmax)])
            yaxis.update(type="log", range=[np.log10(ymin), np.log10(ymax)])
        elif scale == "linear":
            dt = (tmax - tmin)*0.05 if tmax > tmin else max(tmax, 1.0)*0.05
            dp = (pmax - pmin)*0.05 if pmax > pmin else max(pmax, 1.0)*0.05
            xmin = min(tmin, T_cr or tmin) - dt
            xmax = max(tmax, T_cr or tmax) + dt
            ymin = min(pmin, P_cr or pmin) - dp
            ymax = max(pmax, P_cr or pmax) + dp
            xaxis.update(range=[xmin, xmax])
            yaxis.update(range=[ymin, ymax])
        else:
            raise ValueError("scale must be 'log' or 'linear'.")

        self.layout(
            title=title or None,
            xaxis=xaxis,
            yaxis=yaxis,
            legend=dict(title=None),
            margin=dict(l=80, r=30, t=60, b=60),
        )
