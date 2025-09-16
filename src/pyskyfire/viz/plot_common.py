from __future__ import annotations
import plotly.graph_objects as go
from .core import PlotBase
from pyvis.network import Network
import os, numpy as np, html as _html
from typing import Dict, List, Literal, Any, Optional


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

import plotly.graph_objects as go
import numpy as np

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

import plotly.graph_objects as go
import numpy as np

class PlotEngineNetwork(PlotBase):
    """
    Plot an EngineNetwork interactively with Plotly.

    Features:
      • Station label modes: "name" | "values" | "both" | "hidden"
      • Mass-flow–scaled arrow widths (optional)
      • Dashed signal edges
      • Lightweight layout: networkx.spring_layout if available, else circular
    """
    def __init__(self,
                 engine_network,
                 title: str | None = "Engine Network",
                 station_mode: str = "values",
                 mass_flow_based_arrows: bool = False,
                 edge_length: int = 200,
                 template: str = "plotly_white",
                 height: int | None = 800,
                 width: int | None = None):
        super().__init__(go.Figure())
        self.template(template)

        # ---------- helpers ----------
        def _wrap_label(text: str, max_len: int = 25, max_words_line: int = 3) -> str:
            words = text.split()
            if len(text) <= max_len and len(words) <= max_words_line:
                return text
            mid = len(words) // 2
            return "{}<br>{}".format(" ".join(words[:mid]), " ".join(words[mid:]))

        def _nice(name: str) -> str:
            return _wrap_label(name.replace("_", " ").title())

        # colors
        type_colors = {
            "TransmissionBlock": "#FFE066",
            "PumpBlock":         "#C9B3E6",
            "RegenBlock":        "#D46A6A",
            "TurbineBlock":      "#66CDAA",
        }
        default_block_color = "#AFC6E0"
        station_fill = "#FFFFFF"
        edge_color = "#888"

        # width scaling from mdot
        stations = engine_network.stations
        mdot_max = max((s.mdot for s in stations.values() if not np.isnan(s.mdot)), default=1.0)
        def _w(mdot):
            if (not mass_flow_based_arrows) or mdot is None or np.isnan(mdot):
                return 1.5
            # ~1.5–6 px
            return 1.5 + 4.5*np.sqrt(max(mdot, 0.0)/mdot_max)

        # ---------- nodes ----------
        nodes = {}  # name -> dict(x=?, y=?, kind='block'|'station', color=?)
        # blocks
        for blk in engine_network.blocks:
            nodes[blk.name] = {
                "kind": "block",
                "label": _nice(blk.name),
                "color": type_colors.get(blk.__class__.__name__, default_block_color)
            }
        # stations (unless hidden)
        if station_mode != "hidden":
            for name, st in stations.items():
                if station_mode == "name":
                    label = _nice(name)
                elif station_mode == "values":
                    label = f"p={st.p:0.3e} Pa<br>T={st.T:0.1f} K<br>ṁ={st.mdot:0.3f} kg/s"
                else:  # both
                    label = f"{_nice(name)}<br>p={st.p:0.3e} Pa<br>T={st.T:0.1f} K<br>ṁ={st.mdot:0.3f} kg/s"
                nodes[name] = {"kind": "station", "label": label, "color": station_fill}

        # ---------- edges ----------
        flow_edges = []   # (u, v, width, title)
        signal_edges = [] # (u, v)

        # station connectivity maps
        prod_map = {k: [] for k in stations}
        cons_map = {k: [] for k in stations}
        for blk in engine_network.blocks:
            for s in getattr(blk, "station_outputs", []) or []:
                prod_map.setdefault(s, []).append(blk.name)
            for s in getattr(blk, "station_inputs", []) or []:
                cons_map.setdefault(s, []).append(blk.name)

        if station_mode == "hidden":
            seen = set()
            for s_key, producers in prod_map.items():
                consumers = cons_map.get(s_key, [])
                if not producers or not consumers:
                    continue
                md = stations.get(s_key).mdot if s_key in stations else np.nan
                for u in producers:
                    for v in consumers:
                        sig = (u, v, s_key)
                        if sig in seen: continue
                        seen.add(sig)
                        flow_edges.append((u, v, _w(md), f"via {s_key}"))
        else:
            for blk in engine_network.blocks:
                for s_in in getattr(blk, "station_inputs", []) or []:
                    md = stations.get(s_in).mdot if s_in in stations else np.nan
                    flow_edges.append((s_in, blk.name, _w(md), None))
                for s_out in getattr(blk, "station_outputs", []) or []:
                    md = stations.get(s_out).mdot if s_out in stations else np.nan
                    flow_edges.append((blk.name, s_out, _w(md), None))

        # signals (dashed, no arrows)
        sig_src = {}
        for blk in engine_network.blocks:
            for s in getattr(blk, "signal_outputs", []) or []:
                sig_src[s] = blk.name
        for blk in engine_network.blocks:
            for s in getattr(blk, "signal_inputs", []) or []:
                u = sig_src.get(s)
                if u:
                    signal_edges.append((u, blk.name))

        # ---------- layout (spring if available, else circular) ----------
        names = list(nodes.keys())
        idx = {n: i for i, n in enumerate(names)}
        N = len(names)
        pos = np.zeros((N, 2), float)

        try:
            import networkx as nx  # optional
            G = nx.DiGraph()
            G.add_nodes_from(names)
            G.add_edges_from([(u, v) for (u, v, *_ ) in flow_edges] + signal_edges)
            # k scaled by "edge_length" to give breathing room
            pos_dict = nx.spring_layout(G, k=max(edge_length, 80)/800.0, iterations=300, seed=1)
            for n, (x, y) in pos_dict.items():
                pos[idx[n], :] = (x, y)
        except Exception:
            # fallback: deterministic circle
            theta = np.linspace(0, 2*np.pi, N, endpoint=False)
            pos[:, 0] = np.cos(theta)
            pos[:, 1] = np.sin(theta)

        # scale/center
        # (Plotly auto-scales; keep coordinates as-is)
        for n, i in idx.items():
            nodes[n]["x"], nodes[n]["y"] = pos[i, 0], pos[i, 1]

        # ---------- traces: nodes ----------
        # blocks
        bx, by, btxt, bcolors = [], [], [], []
        # stations
        sx, sy, stxt = [], [], []
        for n, data in nodes.items():
            if data["kind"] == "block":
                bx.append(data["x"]); by.append(data["y"])
                btxt.append(data["label"]); bcolors.append(data["color"])
            else:
                sx.append(data["x"]); sy.append(data["y"])
                stxt.append(data["label"])

        if bx:
            self.fig.add_trace(go.Scatter(
                x=bx, y=by, mode="markers+text",
                text=btxt, textposition="middle center",
                hoverinfo="skip",
                marker=dict(symbol="square", size=28, color=bcolors, line=dict(color="#444", width=1.5)),
                name="Blocks",
                showlegend=False
            ))
        if sx:
            self.fig.add_trace(go.Scatter(
                x=sx, y=sy, mode="markers+text",
                text=stxt, textposition="middle center",
                hoverinfo="skip",
                marker=dict(symbol="square", size=26, color=station_fill, line=dict(color="#555", width=1.3)),
                name="Stations",
                showlegend=False
            ))

        # ---------- edges: flow (arrows via annotations) ----------
        for (u, v, w, ttl) in flow_edges:
            xu, yu = nodes[u]["x"], nodes[u]["y"]
            xv, yv = nodes[v]["x"], nodes[v]["y"]

            ann_kw = dict(
                x=xv, y=yv, ax=xu, ay=yu,
                xref="x", yref="y", axref="x", ayref="y",
                showarrow=True, standoff=8,
                arrowhead=3, arrowsize=1, arrowwidth=w, arrowcolor=edge_color,
                # text is blank (we only want a tooltip)
                text=""
            )
            if ttl:
                # Annotations support hovertext; don't set hoverlabel.namelength here
                ann_kw["hovertext"] = ttl

            self.fig.add_annotation(**ann_kw)

            # thin underlying line for visibility / export
            self.fig.add_shape(
                type="line",
                x0=xu, y0=yu, x1=xv, y1=yv,
                line=dict(width=max(1.0, w*0.5), color=edge_color)
            )
        # ---------- edges: signals (dashed, no arrows) ----------
        for (u, v) in signal_edges:
            xu, yu = nodes[u]["x"], nodes[u]["y"]
            xv, yv = nodes[v]["x"], nodes[v]["y"]
            self.fig.add_shape(
                type="line",
                x0=xu, y0=yu, x1=xv, y1=yv,
                line=dict(width=1.2, color=edge_color, dash="dot")
            )

        # ---------- layout ----------
        self.layout(
            title=title or None,
            xaxis=dict(visible=False),
            yaxis=dict(visible=False, scaleanchor="x", scaleratio=1),
            showlegend=False,
            height=height,
            width=width,
            margin=dict(l=10, r=10, t=50 if title else 10, b=10),
        )


# ----------------------------------------------------------------------------
# Interactive Network Visualization
# ----------------------------------------------------------------------------

def render_engine_network(
    engine_network,
    *,
    height: str = "900px",
    width: str = "100%",
    edge_length: int = 200,
    physics_settings: Dict | str | None = "default",
    mass_flow_based_arrows: bool = False,
    station_mode: Literal["name", "values", "both", "hidden"] = "values",
) -> str:
    """
    Build a PyVis network and return an <iframe srcdoc="..."> HTML snippet
    suited for Report.add_raw_html(...). Produces a single-file embed.
    """
    # ---------------------- PyVis canvas & physics ----------------------
    net = Network(height=height, width=width, directed=True, bgcolor="#ffffff", font_color="#222")
    if isinstance(physics_settings, dict):
        net.barnes_hut(**physics_settings)
    elif physics_settings == "default":
        net.barnes_hut(
            gravity=-80000,
            central_gravity=0.3,
            spring_length=edge_length,
            spring_strength=1.0e-3,
            damping=0.09,
            overlap=0,
        )



    # ---------------------- Styles / colors -----------------------------
    type_colours: Dict[str, str] = {
        "TransmissionBlock": "#FFE066",
        "PumpBlock":         "#C9B3E6",
        "RegenBlock":        "#D46A6A",
        "TurbineBlock":      "#66CDAA",
    }
    default_block_colour = "#AFC6E0"
    station_colour = "#FFFFFF"
    edge_colour = "#888"

    font_block   = {"color": "#222", "face": "Latin Modern Roman"}
    font_station = {"align": "left", "color": "#333", "face": "Latin Modern Roman"}

    # ---------------------- Helpers ------------------------------------
    def _wrap_label(text: str, max_len: int = 25, max_words_line: int = 3) -> str:
        words = text.split()
        if len(text) <= max_len and len(words) <= max_words_line:
            return text
        mid = len(words) // 2
        return "{}\n{}".format(" ".join(words[:mid]), " ".join(words[mid:]))

    def _nice(name: str) -> str:
        return _wrap_label(name.replace("_", " ").title())

    stations = engine_network.stations
    mdot_max = max((s.mdot for s in stations.values() if not np.isnan(s.mdot)), default=1.0)

    def _width(mdot: float | int | None) -> int:
        if not mass_flow_based_arrows or mdot is None or np.isnan(mdot):
            return 1
        return int(1 + 5 * np.sqrt(max(mdot, 0.0) / mdot_max))  # 1–6 px

    # ---------------------- Station maps (for hidden mode) --------------
    prod_map: Dict[str, List[str]] = {key: [] for key in stations}
    cons_map: Dict[str, List[str]] = {key: [] for key in stations}
    for blk in engine_network.blocks:
        for st_key in getattr(blk, "station_outputs", []) or []:
            prod_map.setdefault(st_key, []).append(blk.name)
        for st_key in getattr(blk, "station_inputs", []) or []:
            cons_map.setdefault(st_key, []).append(blk.name)

    # ---------------------- Add STATION nodes ---------------------------
    if station_mode != "hidden":
        for name, st in stations.items():
            if station_mode == "name":
                label = _nice(name)
            elif station_mode == "values":
                label = f"p = {st.p:0.3e} Pa\nT = {st.T:0.1f} K\nṁ = {st.mdot:0.3f} kg/s"
            else:  # both
                label = f"{_nice(name)}\n" \
                        f"p = {st.p:0.3e} Pa\nT = {st.T:0.1f} K\nṁ = {st.mdot:0.3f} kg/s"

            net.add_node(
                name,
                label=label,
                shape="box",
                color={"background": station_colour, "border": "#555"},
                font=font_station,
            )

    # ---------------------- Add BLOCK nodes -----------------------------
    for blk in engine_network.blocks:
        colour = type_colours.get(blk.__class__.__name__, default_block_colour)
        net.add_node(
            blk.name,
            label=_nice(blk.name),
            shape="box",
            shape_properties={"borderRadius": 10},
            color={"background": colour, "border": "#444"},
            font=font_block,
        )

    # ---------------------- Add FLOW edges ------------------------------
    if station_mode == "hidden":
        # block→block edges via common station
        seen: set[tuple[str, str, str]] = set()
        for st_key, producers in prod_map.items():
            consumers = cons_map.get(st_key, [])
            if not producers or not consumers:
                continue
            mdot = stations.get(st_key).mdot if st_key in stations else np.nan
            for p in producers:
                for c in consumers:
                    sig = (p, c, st_key)
                    if sig in seen: 
                        continue
                    seen.add(sig)
                    net.add_edge(
                        p, c,
                        arrows="to",
                        length=edge_length // 2,
                        color=edge_colour,
                        width=_width(mdot),
                        title=f"via {st_key}",
                    )
    else:
        # station↔block edges
        for blk in engine_network.blocks:
            for st_key in getattr(blk, "station_inputs", []) or []:
                mdot = stations.get(st_key).mdot if st_key in stations else np.nan
                net.add_edge(st_key, blk.name, arrows="to", length=edge_length // 2,
                             color=edge_colour, width=_width(mdot))
            for st_key in getattr(blk, "station_outputs", []) or []:
                mdot = stations.get(st_key).mdot if st_key in stations else np.nan
                net.add_edge(blk.name, st_key, arrows="to", length=edge_length // 2,
                             color=edge_colour, width=_width(mdot))

    # ---------------------- Dashed SIGNAL edges -------------------------
    sig_src: Dict[str, str] = {}
    for blk in engine_network.blocks:
        for sig in getattr(blk, "signal_outputs", []) or []:
            sig_src[sig] = blk.name
    for blk in engine_network.blocks:
        for sig in getattr(blk, "signal_inputs", []) or []:
            src = sig_src.get(sig)
            if src:
                net.add_edge(src, blk.name, arrows="to", length=edge_length // 2,
                             color=edge_colour, width=1, dashes=True)

    # ---------------------- Generate standalone HTML --------------------
    # Prefer inline HTML for single-file embedding
    try:
        html_full = net.generate_html(notebook=False)
    except Exception:
        # Fallback for older pyvis: write then read
        import tempfile
        fd, tmp = None, None
        try:
            import tempfile, io
            with tempfile.NamedTemporaryFile("w+", suffix=".html", delete=False, encoding="utf-8") as f:
                tmp = f.name
                net.write_html(tmp, open_browser=False)
                f.flush()
            with open(tmp, "r", encoding="utf-8") as f:
                html_full = f.read()
        finally:
            if tmp and os.path.exists(tmp):
                try: os.remove(tmp)
                except Exception: pass

    # Wrap into an iframe via srcdoc so it nests cleanly inside the Report
    srcdoc = _html.escape(html_full, quote=True)
    iframe = (
        f"<iframe "
        f"style='width:100%;height:{_html.escape(height)};border:0;' "
        f"sandbox='allow-scripts allow-same-origin' "
        f"srcdoc=\"{srcdoc}\"></iframe>"
    )
    return iframe
