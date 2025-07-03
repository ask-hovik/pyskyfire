from __future__ import annotations
import matplotlib.pyplot as plt
import pyskyfire as psf
import os
import matplotlib.pyplot as plt
from collections.abc import Iterable
from CoolProp.CoolProp import PropsSI
from pyvis.network import Network
from typing import Literal, Dict, List
import numpy as np
import itertools


# Apply project-wide professional plotting style
plt.style.use('ggplot')
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Latin Modern Roman'],
    'font.size': 10,
    'savefig.dpi': 600
})

def plot_residual_history(residuals,
                          title="Residual Convergence History",
                          color='C0',
                          linestyle='-',
                          marker='o'):
    """
    Plot the convergence history of residuals on a semilogarithmic y-axis,
    styled consistently with the project's plotting conventions.

    Parameters:
        residuals (list of float): Residual values per iteration.
        title (str): Title for the plot.
        color (str): Matplotlib color spec for the line.
        linestyle (str): Matplotlib linestyle for the line.
        marker (str): Matplotlib marker style for the data points.

    Returns:
        fig, ax: Matplotlib Figure and Axes objects for further customization or saving.
    """
    # Set up figure and axis
    fig, ax = plt.subplots(figsize=(7, 5))

    # X = iteration count, Y = residual
    iterations = list(range(1, len(residuals) + 1))
    ax.semilogy(iterations,
                residuals,
                color=color,
                linestyle=linestyle,
                marker=marker,
                linewidth=2,
                label='Residual')

    # Labels and title
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Residual (dimensionless)")
    ax.set_title(title, pad=12)

    # Grid for readability
    ax.grid(True, which='both', linestyle='--', alpha=0.7)

    # Hide top and right spines for cleaner look
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Thicker ticks
    ax.tick_params(width=1.2)

    # Legend box
    ax.legend(frameon=True, framealpha=0.9, edgecolor='black')

    plt.tight_layout()
    return fig, ax


def plot_station_property(
        station_dicts,
        station_list,
        property_name,
        labels=None,
        title=True,
        ax=None,
        ylabel=None,
        color_cycle=None,
        style_cycle=None,
        marker_cycle=None,
        ylim=None):
    """
    Plot *property_name* for an ordered list of stations taken from one
    **or many** station dictionaries.

    title parameter behavior:
      • None  → no title
      • True  → default title: "{property_name} vs Station"
      • str   → use the given string

    Other params as before.
    """
    # Normalize station_dicts input
    from collections.abc import Iterable
    if not isinstance(station_dicts, Iterable) or isinstance(station_dicts, dict):
        station_dicts = [station_dicts]
    n_cases = len(station_dicts)

    # Labels
    if labels is False:
        do_legend   = False
        labels_list = [None] * n_cases
    elif labels is True or labels is None:
        do_legend   = True
        labels_list = [f"Case {i}" for i in range(n_cases)]
    else:
        # assume it's a list
        if len(labels) != n_cases:
            raise ValueError("len(labels) must match number of station_dicts")
        do_legend   = True
        labels_list = labels

    # Color cycle: default from rcParams if not provided
    default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    color_cycle = list(color_cycle or default_colors)
    if len(color_cycle) < n_cases:
        color_cycle = list(itertools.islice(itertools.cycle(color_cycle), n_cases))

    # Line style cycle
    style_cycle = list(style_cycle or ["-"] * n_cases)
    if len(style_cycle) < n_cases:
        raise ValueError(f"Need at least {n_cases} linestyles, got {len(style_cycle)}")

    # Marker cycle
    default_markers = ["o", "s", "^", "d", "v", "P", "X"]
    marker_cycle = list(marker_cycle or default_markers)
    if len(marker_cycle) < n_cases:
        marker_cycle = list(itertools.islice(itertools.cycle(marker_cycle), n_cases))

    # Extract y-values
    def _extract(st, name):
        if isinstance(st, psf.common.Station):
            if not hasattr(st, property_name):
                raise KeyError(f"Station '{name}' lacks attribute '{property_name}'.")
            return getattr(st, property_name)
        elif isinstance(st, dict):
            try:
                return st[property_name]
            except KeyError as exc:
                raise KeyError(f"Property '{property_name}' missing for station '{name}'.") from exc
        else:
            raise TypeError(f"Unsupported station entry type {type(st)} for '{name}'.")

    all_y = [
        [_extract(d[name], name) for name in station_list]
        for d in station_dicts
    ]

    # Begin plotting
    x_vals = list(range(1, len(station_list) + 1))

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5))
    else:
        fig = ax.figure

    for y, lab, c, ls, mk in zip(all_y,
                                 labels_list,
                                 color_cycle,
                                 style_cycle,
                                 marker_cycle):
        if do_legend:
            ax.plot(x_vals, y,
                    color=c, linestyle=ls, marker=mk,
                    linewidth=2, label=lab)
        else:
            ax.plot(x_vals, y,
                    color=c, linestyle=ls, marker=mk,
                    linewidth=2)

    # Axes formatting
    ax.set_xticks(x_vals)
    # Format station names: replace underscores, capitalize words
    formatted_labels = [name.replace('_', ' ').title() for name in station_list]
    ax.set_xticklabels(formatted_labels, rotation=90)
    #ax.set_xlabel('Station')
    ax.set_ylabel(ylabel or property_name)

    # Title logic
    if title is True:
        ax.set_title(f"{property_name} vs Station", pad=12)
    elif isinstance(title, str):
        ax.set_title(title, pad=12)
    # else title is None → skip setting any title

    # Grid, spines, ticks, legend
    ax.grid(True, which='both', linestyle='--', alpha=0.7)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(width=1.2)
    ax.legend(frameon=True, framealpha=0.9, edgecolor='black')
    if ylim is not None: 
        ax.set_ylim(ylim)

    plt.tight_layout()
    return fig, ax



def plot_PT_diagram(
        station_dicts,
        station_list,
        fluid_name=None,
        title=None,
        sat_points=200,
        labels=None,
        color_cycle=None,
        style_cycle=None,
        marker_cycle=None,
        annotate_ha=None,
        annotate_va=None,
        legend_loc='best',
        scale='log'):
    """
    Plot P–T diagram for one or multiple station runs,
    styled consistently with project conventions.

    Parameters:
    -----------
    station_dicts : dict or iterable of dicts
        Single station mapping or list of mappings:
        {station_name: {'T': K, 'p': Pa, ...}, ...}
    station_list : list[str]
        Ordered list of station names to plot.
    fluid_name : str, optional
        Pure fluid name for CoolProp saturation curve.
    title : str or None
        Plot title.
    sat_points : int, optional
        Number of points for saturation curve.
    labels : list[str], optional
        Legend labels for each case (defaults to "Case 0", ...).
    color_cycle : list[str], optional
        Color cycle for each case; defaults to matplotlib rcParams cycle.
    style_cycle : list[str], optional
        Line style cycle per case.
    marker_cycle : list[str], optional
        Marker cycle per case.
    annotate_ha : list[str] or str, optional
        Horizontal alignments for station labels.
    annotate_va : list[str] or str, optional
        Vertical alignments for station labels.
    legend_loc : str, optional
        Legend location code (Matplotlib loc string).
    scale : {'log', 'linear'}
        Axis scale for both x and y.

    Returns:
    --------
    fig, ax : matplotlib Figure and Axes
    """
    from collections.abc import Iterable
    import numpy as np
    from CoolProp.CoolProp import PropsSI

    # Normalize input
    if not isinstance(station_dicts, Iterable) or isinstance(station_dicts, dict):
        station_dicts = [station_dicts]
    n_cases = len(station_dicts)

    # Legend labels
    labels = labels or [f"Case {i}" for i in range(n_cases)]
    if len(labels) != n_cases:
        raise ValueError("len(labels) must equal number of station dictionaries.")

    # Color cycle
    default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    colors = list(color_cycle or default_colors)
    if len(colors) < n_cases:
        colors = list(itertools.islice(itertools.cycle(colors), n_cases))

    # Style and marker cycles
    style_cycle = style_cycle or ["-"] * n_cases
    marker_cycle = marker_cycle or ["o", "s", "^", "d", "v", "P", "X"] * ((n_cases + 6) // 7)

    # Extract T and p for each case
    all_T, all_P = [], []
    for d in station_dicts:
        Ts, Ps = [], []
        for name in station_list:
            if name not in d:
                raise KeyError(f"Station '{name}' not found in station_dict.")
            val = d[name]
            if hasattr(val, 'T') and hasattr(val, 'p'):
                Ts.append(val.T)
                Ps.append(val.p)
            elif isinstance(val, dict):
                Ts.append(val['T'])
                Ps.append(val['p'])
            else:
                raise TypeError(f"Unsupported station entry type '{type(val)}' for '{name}'.")
        all_T.append(Ts)
        all_P.append(Ps)

    # Prepare annotation alignments
    def _broadcast(arg, name):
        if isinstance(arg, str) or arg is None:
            return [arg] * len(station_list)
        if len(arg) != len(station_list):
            raise ValueError(f"{name} must be a string or list same length as station_list.")
        return arg

    has = _broadcast(annotate_ha or 'left', 'annotate_ha')
    vas = _broadcast(annotate_va or 'top', 'annotate_va')

    # Helper to compute offset points based on ha/va
    def _offset(ha, va, delta=5):
        dx = {'left': delta, 'center': 0, 'right': -delta}.get(ha, 0)
        dy = {'bottom': delta, 'center': 0, 'top': -delta}.get(va, 0)
        return (dx, dy)

    # Plot
    fig, ax = plt.subplots(figsize=(7, 5))
    for Ts, Ps, lab, c, ls, mk in zip(all_T, all_P, labels, colors, style_cycle, marker_cycle):
        ax.plot(Ts, Ps, color=c, linestyle=ls, marker=mk, linewidth=2, label=lab)
        for name, T, P, ha, va in zip(station_list, Ts, Ps, has, vas):
            off = _offset(ha, va)
            ax.annotate(name.replace('_', ' ').title(), xy=(T, P), xytext=off, textcoords='offset points', ha=ha, va=va, fontsize=8)

    # Saturation curve in black
    T_cr = P_cr = None
    if fluid_name:
        T_tr = PropsSI('Ttriple', fluid_name)
        T_cr = PropsSI('Tcrit', fluid_name)
        T_sat = np.linspace(T_tr * 1.01, T_cr * 0.99, sat_points)
        P_sat = [PropsSI('P', 'T', T, 'Q', 0, fluid_name) for T in T_sat]

        ax.plot(T_sat, P_sat, linestyle='--', color='black', linewidth=2, label=f'{fluid_name} sat. line')
        P_cr = PropsSI('Pcrit', fluid_name)
        ax.scatter([T_cr], [P_cr], marker='x', s=50, color='black')
        # Annotate critical point
        ax.annotate('Critical pt', xy=(T_cr, P_cr), xytext=(5, -5), textcoords='offset points', ha='left', va='top', fontsize=8)

    # scales
    if scale == 'log':
        ax.set_xscale('log')
        ax.set_yscale('log')
        from matplotlib.ticker import FixedLocator, NullFormatter
        # major ticks every 100 K
        from matplotlib.ticker import FixedLocator, NullFormatter, ScalarFormatter
        xmin_lim, xmax_lim = ax.get_xlim()
        major_ticks = np.arange(
            np.ceil(xmin_lim/100)*100,
            np.floor(xmax_lim/100)*100 + 100,
            100)
        ax.xaxis.set_major_locator(FixedLocator(major_ticks))
        # format major ticks as plain numbers
        major_formatter = ScalarFormatter()
        major_formatter.set_scientific(False)
        major_formatter.set_useOffset(False)
        ax.xaxis.set_major_formatter(major_formatter)
        # minor ticks every 10 K
        minor_ticks = np.arange(
            np.ceil(xmin_lim/10)*10,
            np.floor(xmax_lim/10)*10 + 10,
            10)
        ax.xaxis.set_minor_locator(FixedLocator(minor_ticks))
        ax.xaxis.set_minor_formatter(NullFormatter())
    elif scale == 'linear':
        ax.set_xscale('linear')
        ax.set_yscale('linear')
    else:
        raise ValueError("`scale` must be 'log' or 'linear'.")

    # Zoom to station path
    all_T_flat = [t for Ts in all_T for t in Ts]
    all_P_flat = [p for Ps in all_P for p in Ps]
    tmin, tmax = min(all_T_flat), max(all_T_flat)
    pmin, pmax = min(all_P_flat), max(all_P_flat)
    if scale == 'linear':
        dt = (tmax - tmin) * 0.05 if tmax > tmin else tmax * 0.05
        dp = (pmax - pmin) * 0.05 if pmax > pmin else pmax * 0.05
        xmin, xmax = tmin - dt, tmax + dt
        ymin, ymax = pmin - dp, pmax + dp
    else:  # log scale: multiplicative margins
        factor_x = 1.1
        factor_y = 1.4
        xmin, xmax = tmin / factor_x, tmax * factor_x
        ymin, ymax = pmin / factor_y, pmax * factor_y
    # include critical point
    if T_cr is not None and P_cr is not None:
        if scale == 'linear':
            xmin = min(xmin, T_cr - dt)
            xmax = max(xmax, T_cr + dt)
            ymin = min(ymin, P_cr - dp)
            ymax = max(ymax, P_cr + dp)
        else:
            xmin = min(xmin, T_cr / factor_x)
            xmax = max(xmax, T_cr * factor_x)
            ymin = min(ymin, P_cr / factor_y)
            ymax = max(ymax, P_cr * factor_y)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    # labels, title, grid, spines, ticks
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Pressure (Pa)')
    if title:
        ax.set_title(title, pad=12)

    ax.grid(True, which='both', linestyle='--', alpha=0.7)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(width=1.2)

    # legend + layout
    ax.legend(loc=legend_loc, frameon=True, framealpha=0.9, edgecolor='black')
    plt.tight_layout()
    return fig, ax






# ----------------------------------------------------------------------------
# Interactive Network Visualization
# ----------------------------------------------------------------------------



def _wrap_label(text: str, max_len: int = 25, max_words_line: int = 3) -> str:
    """Wrap *text* onto (at most) two lines if it is long.

    Splits approximately in the middle of the word list when either the total
    character length exceeds *max_len* *or* the number of words exceeds
    *max_words_line*.
    """
    words = text.split()
    if len(text) <= max_len and len(words) <= max_words_line:
        return text

    # Split after half the words (second half may be the same size or shorter).
    mid = len(words) // 2
    return "{}\n{}".format(" ".join(words[:mid]), " ".join(words[mid:]))


def _nice(name: str) -> str:
    """Return *name* prettified for display (spaces + title‑case + wrapping)."""
    prettified = name.replace("_", " ").title()
    return _wrap_label(prettified)


def plot_engine_network(
    engine_network,
    output_path: str,
    output_file: str = "engine_network2.html",
    notebook: bool = False,
    height: str = "1080px",
    width: str = "100%",
    edge_length: int = 200,
    physics_settings: Dict | str | None = None,
    mass_flow_based_arrows: bool = False,
    station_mode: Literal["name", "values", "both", "hidden"] = "values",
):
    """Interactive network diagram of an *EngineNetwork* with extensive customisation.

    The function prettifies both **block** and **station** labels by
    1) replacing underscores with spaces, 2) applying *Title Case*, and 3) wrapping
    long labels over two lines for better readability (using
    :func:`_wrap_label`).  Latin Modern Roman is used for all labels.
    """

    # ------------------------------------------------------------------
    #  0. Output file / directory --------------------------------------
    # ------------------------------------------------------------------
    os.makedirs(output_path, exist_ok=True)
    full_output = os.path.join(output_path, output_file)

    # ------------------------------------------------------------------
    #  1. PyVis canvas --------------------------------------------------
    # ------------------------------------------------------------------
    net = Network(height=height, width=width, directed=True)
    if notebook:
        net.toggle_notebook()

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

    # ------------------------------------------------------------------
    #  2. Colour palette & fonts ---------------------------------------
    # ------------------------------------------------------------------
    type_colours: Dict[str, str] = {
        "TransmissionBlock": "#FFE066",
        "PumpBlock": "#C9B3E6",
        "RegenBlock": "#D46A6A",
        "TurbineBlock": "#66CDAA",
    }
    default_block_colour = "#AFC6E0"
    station_colour = "#FFFFFF"
    edge_colour = "#888"

    font_block = {"color": "#222", "face": "Latin Modern Roman"}
    font_station = {"align": "left", "color": "#333", "face": "Latin Modern Roman"}

    # ------------------------------------------------------------------
    #  3. Utility: edge‑width from mass‑flow ---------------------------
    # ------------------------------------------------------------------
    stations = engine_network.stations  # shorthand
    mdot_max = max((s.mdot for s in stations.values() if not np.isnan(s.mdot)), default=1.0)

    def _width(mdot: float | int | None) -> int:
        if not mass_flow_based_arrows or mdot is None or np.isnan(mdot):
            return 1
        return int(1 + 5 * np.sqrt(max(mdot, 0.0) / mdot_max))  # 1–6 px

    # ------------------------------------------------------------------
    #  4. Add STATION nodes & lookup tables ---------------------------
    # ------------------------------------------------------------------
    # Map station → producer blocks & consumer blocks for the *hidden* mode.
    prod_map: Dict[str, List[str]] = {key: [] for key in stations}
    cons_map: Dict[str, List[str]] = {key: [] for key in stations}

    for blk in engine_network.blocks:
        for st_key in getattr(blk, "station_outputs", []) or []:
            prod_map.setdefault(st_key, []).append(blk.name)
        for st_key in getattr(blk, "station_inputs", []) or []:
            cons_map.setdefault(st_key, []).append(blk.name)

    if station_mode != "hidden":
        for name, st in stations.items():
            if station_mode == "name":
                label = _nice(name)
            elif station_mode == "values":
                label = (
                    f"p = {st.p:0.3e} Pa\nT = {st.T:0.1f} K\nṁ = {st.mdot:0.3f} kg/s"
                )
            else:  # both
                label = (
                    f"{_nice(name)}\n"
                    f"p = {st.p:0.3e} Pa\nT = {st.T:0.1f} K\nṁ  = {st.mdot:0.3f} kg/s"
                )

            net.add_node(
                name,  # identifier remains untouched
                label=label,
                shape="box",
                color={"background": station_colour, "border": "#555"},
                font=font_station,
            )

    # ------------------------------------------------------------------
    #  5. Add BLOCK nodes ---------------------------------------------
    # ------------------------------------------------------------------
    for blk in engine_network.blocks:
        colour = type_colours.get(blk.__class__.__name__, default_block_colour)
        net.add_node(
            blk.name,  # identifier (key) unchanged
            label=_nice(blk.name),  # prettified + wrapped label
            shape="box",
            shape_properties={"borderRadius": 10},
            color={"background": colour, "border": "#444"},
            font=font_block,
        )

    # ------------------------------------------------------------------
    #  6. Add FLOW edges ----------------------------------------------
    # ------------------------------------------------------------------
    if station_mode == "hidden":
        # Direct block‑to‑block edges through the common station.
        seen: set[tuple[str, str, str]] = set()  # (prod, cons, station) to avoid dupes
        for st_key, producers in prod_map.items():
            consumers = cons_map.get(st_key, [])
            if not producers or not consumers:
                continue  # avoid dangling edges
            mdot = stations.get(st_key).mdot if st_key in stations else np.nan
            for p in producers:
                for c in consumers:
                    sig = (p, c, st_key)
                    if sig in seen:
                        continue
                    seen.add(sig)
                    net.add_edge(
                        p,
                        c,
                        arrows="to",
                        length=edge_length // 2,
                        color={"color": edge_colour},
                        width=_width(mdot),
                        title=f"via {st_key}",
                    )
    else:
        # Standard: station → block (inlets) and block → station (outlets)
        for blk in engine_network.blocks:
            # inlet
            for st_key in getattr(blk, "station_inputs", []) or []:
                mdot = stations.get(st_key).mdot if st_key in stations else np.nan
                net.add_edge(
                    st_key,
                    blk.name,
                    arrows="to",
                    length=edge_length // 2,
                    color={"color": edge_colour},
                    width=_width(mdot),
                )
            # outlet
            for st_key in getattr(blk, "station_outputs", []) or []:
                mdot = stations.get(st_key).mdot if st_key in stations else np.nan
                net.add_edge(
                    blk.name,
                    st_key,
                    arrows="to",
                    length=edge_length // 2,
                    color={"color": edge_colour},
                    width=_width(mdot),
                )

    # ------------------------------------------------------------------
    #  7. Dashed SIGNAL edges -----------------------------------------
    # ------------------------------------------------------------------
    sig_src: Dict[str, str] = {}
    for blk in engine_network.blocks:
        for sig in getattr(blk, "signal_outputs", []) or []:
            sig_src[sig] = blk.name

    for blk in engine_network.blocks:
        for sig in getattr(blk, "signal_inputs", []) or []:
            src = sig_src.get(sig)
            if src:
                net.add_edge(
                    src,
                    blk.name,
                    arrows="to",
                    dashes=True,
                    length=edge_length // 2,
                    color={"color": edge_colour},
                    width=1,
                )

    # ------------------------------------------------------------------
    #  8. Export -------------------------------------------------------
    # ------------------------------------------------------------------
    net.write_html(full_output, open_browser=False, notebook=notebook)
    return net
