"""Define, solve and draw pipes with rounded and continuous bends.

Holds the code of the "Define and solve a pipe" how-to guide and writes the
figures and 3D scenes embedded in it.

Run from the repository with:
    python examples/pipe/pipe_paths.py
Requires pyskyfire[viz]. No browser or display is needed for HTML output.
"""
import argparse
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.interpolate import CubicSpline

import pyskyfire as psf

NODE_STYLE = dict(color="grey", size=12, symbol="circle-open", line=dict(width=2))
STATION_STYLE = dict(color="#327fa8", size=8, symbol="x")
BEND_COLOR = "#e58b27"


def _bend_arc(pipe):
    """Display points of the pipe's lumped rounded bend, from tangent to tangent."""
    display = pipe.display_path(segments=48)
    bend = next(i for i, c in enumerate(pipe.connections)
                if isinstance(c, psf.ducts.RoundedBend))
    segments = np.flatnonzero(display.connection == bend)
    return display.points[np.r_[segments[0], segments + 1]]


def _walls(pipe):
    """Both bore walls in the y-z plane the example paths lie in."""
    display = pipe.display_path(segments=48)
    yz = display.points[:, 1:]
    legs = np.diff(yz, axis=0)
    legs /= np.linalg.norm(legs, axis=1)[:, None]
    tangent = np.vstack((legs[0], legs[:-1] + legs[1:], legs[-1]))
    tangent /= np.linalg.norm(tangent, axis=1)[:, None]
    offset = display.radius[:, None] * np.column_stack((-tangent[:, 1], tangent[:, 0]))
    return yz + offset, yz - offset


def plot_path_construction(nodes, pipe):
    """The PathNodes given, the PathStations generated from them, and the bend arc."""
    nodes = np.asarray([node.position for node in nodes])
    stations = pipe.centerline
    arc = _bend_arc(pipe)

    # The path lies in the y-z plane, so a 2D view shows it without distortion.
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=nodes[:, 1], y=nodes[:, 2], mode="lines+markers", name="PathNode",
        line=dict(color="grey", width=1, dash="dot"), marker=NODE_STYLE,
    ))
    fig.add_trace(go.Scatter(
        x=arc[:, 1], y=arc[:, 2], mode="lines",
        name="RoundedBend (one connection)",
        line=dict(color=BEND_COLOR, width=3, dash="dash"),
    ))
    fig.add_trace(go.Scatter(
        x=stations[:, 1], y=stations[:, 2], mode="markers", name="PathStation",
        marker=STATION_STYLE,
    ))
    fig.update_layout(
        title="PathNodes and the PathStations generated from them",
        xaxis_title="y [m]",
        yaxis=dict(title="z [m]", scaleanchor="x", scaleratio=1),
    )
    return fig


def plot_bend_comparison(nodes, rounded, continuous):
    """The same corner as one lumped bend and as a resolved, flaring arc."""
    nodes = np.asarray([node.position for node in nodes])
    fig = make_subplots(
        rows=1, cols=2, horizontal_spacing=0.08,
        subplot_titles=(
            f"RoundedBend: {len(rounded.stations)} PathStations",
            f"Continuous: {len(continuous.stations)} PathStations",
        ),
    )
    for col, pipe in enumerate((rounded, continuous), start=1):
        fig.add_trace(go.Scatter(
            x=nodes[:, 1], y=nodes[:, 2], mode="lines+markers", name="PathNode",
            legendgroup="nodes", showlegend=col == 1,
            line=dict(color="lightgrey", width=1, dash="dot"),
            marker={**NODE_STYLE, "color": "lightgrey"},
        ), row=1, col=col)
        for side, wall in enumerate(_walls(pipe)):
            fig.add_trace(go.Scatter(
                x=wall[:, 0], y=wall[:, 1], mode="lines", name="Bore wall",
                legendgroup="wall", showlegend=col == 1 and side == 0,
                line=dict(color="grey", width=1), hoverinfo="skip",
            ), row=1, col=col)
        stations = pipe.centerline
        if col == 1:
            arc = _bend_arc(pipe)
            fig.add_trace(go.Scatter(
                x=arc[:, 1], y=arc[:, 2], mode="lines",
                name="RoundedBend (one connection)",
                line=dict(color=BEND_COLOR, width=3, dash="dash"),
            ), row=1, col=col)
        else:
            bend = [i for i, c in enumerate(pipe.connections)
                    if isinstance(c, psf.ducts.Continuous)]
            arc = stations[bend[0]:bend[-1] + 2]
            fig.add_trace(go.Scatter(
                x=arc[:, 1], y=arc[:, 2], mode="lines",
                name="Continuous (one connection per station pair)",
                line=dict(color=BEND_COLOR, width=3),
            ), row=1, col=col)
        fig.add_trace(go.Scatter(
            x=stations[:, 1], y=stations[:, 2], mode="markers", name="PathStation",
            legendgroup="stations", showlegend=col == 1,
            marker={**STATION_STYLE, "size": 6},
            text=[f"PathStation {i}, D = {d * 1e3:.1f} mm"
                  for i, d in enumerate(pipe.diameter)],
            hovertemplate="%{text}<extra></extra>",
        ), row=1, col=col)
    # Same ranges in both panels, so the bores compare at one scale.
    fig.update_xaxes(title_text="y [m]", range=[-0.15, 0.65])
    fig.update_yaxes(range=[-0.05, 0.85])
    fig.update_yaxes(title_text="z [m]", scaleanchor="x", scaleratio=1, row=1, col=1)
    fig.update_yaxes(scaleanchor="x2", scaleratio=1, row=1, col=2)
    fig.update_layout(legend=dict(orientation="h", y=-0.2))
    return fig


def plot_spline_curvature(spline, pipe):
    """The curvature of each fitted connection, against the spline's own."""
    t = np.linspace(spline.x[0], spline.x[-1], 2001)
    first, second = spline(t, 1), spline(t, 2)
    exact = (np.linalg.norm(np.cross(first, second), axis=1)
             / np.linalg.norm(first, axis=1)**3)
    # The spline parameter is chord length, so measure true arc length here.
    speed = np.linalg.norm(first, axis=1)
    s = np.concatenate(([0.0], np.cumsum((speed[1:] + speed[:-1]) / 2 * np.diff(t))))
    path = pipe.path
    midpoints = (path.s[:-1] + path.s[1:]) / 2

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=s, y=exact, mode="lines", name="CubicSpline curvature",
                             line=dict(color="grey")))
    fig.add_trace(go.Scatter(x=midpoints, y=1 / path.curvature_radius, mode="markers",
                             name="Continuous connections", marker=STATION_STYLE))
    fig.update_layout(
        title="Curvature along the spline pipe",
        xaxis_title="Arc length, s [m]",
        yaxis=dict(title="Curvature, 1/R [1/m]", rangemode="tozero"),
        template="plotly_white",
    )
    return fig


def plot_manual_path(path, radius):
    """Hand-placed PathStations, the connections joining them, and the true arc."""
    stations = path.centerline
    angle = np.linspace(0.0, np.pi / 2, 100)
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=radius * (1 - np.cos(angle)), y=0.3 + radius * np.sin(angle), mode="lines",
        name=f"True arc, R = {radius} m", line=dict(color="lightgrey", width=6),
        hoverinfo="skip",
    ))
    for kind, color in ((psf.ducts.Straight, "grey"), (psf.ducts.Continuous, BEND_COLOR)):
        shown = False
        for i, connection in enumerate(path.connections):
            if not isinstance(connection, kind):
                continue
            fig.add_trace(go.Scatter(
                x=stations[i:i + 2, 1], y=stations[i:i + 2, 2], mode="lines",
                name=kind.__name__, legendgroup=kind.__name__, showlegend=not shown,
                line=dict(color=color, width=3), hoverinfo="skip",
            ))
            shown = True
    fig.add_trace(go.Scatter(
        x=stations[:, 1], y=stations[:, 2], mode="markers", name="PathStation",
        marker=STATION_STYLE, text=[f"PathStation {i}" for i in range(len(stations))],
        hovertemplate="%{text}<extra></extra>",
    ))
    fig.update_layout(
        title="PathStations placed by hand, joined by the connections given",
        xaxis_title="y [m]",
        yaxis=dict(title="z [m]", scaleanchor="x", scaleratio=1),
    )
    return fig


def plot_manual_curvature(estimated, given, radius):
    """Each connection's curvature, estimated from its neighbours and given outright."""
    fig = go.Figure()
    s = estimated.s
    exact = [1 / radius if isinstance(c, psf.ducts.Continuous) else 0.0
             for c in estimated.connections]
    fig.add_trace(go.Scatter(
        x=np.repeat(s, 2)[1:-1], y=np.repeat(exact, 2), mode="lines",
        name="True curvature", line=dict(color="grey"),
    ))
    midpoints = (s[:-1] + s[1:]) / 2
    for path, name, marker in (
        (given, "Continuous(radius=0.2)", dict(color="#2a9d8f", size=11, symbol="circle-open",
                                               line=dict(width=2))),
        (estimated, "Continuous(), estimated", STATION_STYLE),
    ):
        fig.add_trace(go.Scatter(
            x=midpoints, y=1 / path.curvature_radius, mode="markers", name=name, marker=marker,
        ))
    fig.update_layout(
        title="Curvature of each connection",
        xaxis_title="Arc length, s [m]",
        yaxis_title="Curvature, 1/R [1/m]",
        template="plotly_white",
    )
    return fig


def plot_tapered_stations(supplied, pipe):
    """The two stations given and the stations the pipe placed along its taper."""
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=pipe.path.s, y=pipe.diameter * 1e3, mode="lines+markers",
        name=f"Resolved: {len(pipe.stations)} PathStations",
        line=dict(color="grey"), marker={**STATION_STYLE, "size": 7},
    ))
    fig.add_trace(go.Scatter(
        x=supplied.s, y=np.array([0.010, 0.012]) * 1e3, mode="markers",
        name=f"Supplied: {len(supplied.stations)} PathStations", marker=NODE_STYLE,
    ))
    fig.update_layout(
        title="Stations added where the bore changes",
        xaxis_title="Arc length, s [m]",
        yaxis_title="Bore diameter [mm]",
        template="plotly_white",
    )
    return fig


def main(output_dir: Path | None = None) -> None:
    if output_dir is None:
        output_dir = Path(__file__).parent / "output"
    output_dir.mkdir(parents=True, exist_ok=True)

    # tutorial:start:rounded
    nodes = [
        psf.ducts.PathNode((0.0, 0.0, 0.0)),
        psf.ducts.PathNode((0.0, 0.0, 0.5), bend=psf.ducts.RoundedBend(radius=0.25)),
        psf.ducts.PathNode((0.0, 0.5, 0.75)),
    ]

    path = psf.ducts.Path.from_nodes(nodes)
    rounded_pipe = psf.ducts.Pipe(path, diameter=0.025, roughness=15e-6)
    # tutorial:end:rounded

    plot_path_construction(nodes, rounded_pipe).write_html(
        output_dir / "rounded-bend-path.html", include_plotlyjs="cdn")

    # tutorial:start:scene-rounded
    scene = psf.viz.Scene()
    scene.add(psf.viz.PipeComponent(rounded_pipe), name="rounded_bend")
    scene.save_html(output_dir / "rounded-bend-pipe.html")
    # tutorial:end:scene-rounded

    # tutorial:start:flared
    nodes = [
        psf.ducts.PathNode((0.0, 0.0, 0.0)),
        psf.ducts.PathNode((0.0, 0.0, 0.5), bend=psf.ducts.Continuous(radius=0.25)),
        psf.ducts.PathNode((0.0, 0.5, 0.75)),
    ]

    path = psf.ducts.Path.from_nodes(nodes)
    length = path.s[-1]

    def flared_diameter(s):
        """Bore flaring linearly from 25 mm at the inlet to 90 mm at the outlet."""
        return 0.025 + 0.065 * s / length

    flared_pipe = psf.ducts.Pipe(path, diameter=flared_diameter, roughness=15e-6)
    # tutorial:end:flared

    plot_bend_comparison(nodes, rounded_pipe, flared_pipe).write_html(
        output_dir / "bend-comparison.html", include_plotlyjs="cdn")

    scene = psf.viz.Scene()
    scene.add(psf.viz.PipeComponent(flared_pipe), name="flared_bend")
    scene.save_html(output_dir / "flared-bend-pipe.html")

    # tutorial:start:stations
    radius = 0.2
    angles = np.linspace(0.0, np.pi / 2, 7)
    arc = [psf.ducts.PathStation((0.0, radius * (1 - np.cos(a)), 0.3 + radius * np.sin(a)))
           for a in angles]
    stations = [psf.ducts.PathStation((0.0, 0.0, 0.0)), *arc,
                psf.ducts.PathStation((0.0, 0.6, 0.5))]

    connections = [psf.ducts.Straight(), *[psf.ducts.Continuous()] * 6, psf.ducts.Straight()]
    manual_pipe = psf.ducts.Pipe(psf.ducts.Path(stations, connections),
                                 diameter=0.025, roughness=15e-6)
    # tutorial:end:stations

    # tutorial:start:stations-radius
    connections = [psf.ducts.Straight(), *[psf.ducts.Continuous(radius=radius)] * 6,
                   psf.ducts.Straight()]
    known_pipe = psf.ducts.Pipe(psf.ducts.Path(stations, connections),
                                diameter=0.025, roughness=15e-6)
    # tutorial:end:stations-radius

    plot_manual_path(manual_pipe.path, radius).write_html(
        output_dir / "manual-path.html", include_plotlyjs="cdn")
    plot_manual_curvature(manual_pipe.path, known_pipe.path, radius).write_html(
        output_dir / "manual-curvature.html", include_plotlyjs="cdn")

    # tutorial:start:tapered
    supplied = psf.ducts.Path(
        [psf.ducts.PathStation((0.0, 0.0, 0.0)), psf.ducts.PathStation((0.0, 0.0, 0.5))],
        [psf.ducts.Straight()],
    )
    tapered_pipe = psf.ducts.Pipe(supplied, diameter=[0.010, 0.012], roughness=15e-6)
    # tutorial:end:tapered

    plot_tapered_stations(supplied, tapered_pipe).write_html(
        output_dir / "tapered-stations.html", include_plotlyjs="cdn")

    # tutorial:start:spline
    # From a flange at the origin facing up (+z) to one at (0.5, 0.4, 0.6)
    # facing along +x, through one waypoint between them.
    points = np.array([
        [0.00, 0.00, 0.00],
        [0.15, 0.25, 0.40],
        [0.50, 0.40, 0.60],
    ])

    # Parametrise by chord length, so that the parameter runs roughly with
    # distance along the curve, and fix the direction at each flange.
    chord = np.linalg.norm(np.diff(points, axis=0), axis=1)
    t = np.concatenate(([0.0], np.cumsum(chord)))
    spline = CubicSpline(t, points, bc_type=((1, [0.0, 0.0, 1.0]), (1, [1.0, 0.0, 0.0])))

    path = psf.ducts.Path.from_curve(spline, t_span=(t[0], t[-1]))
    spline_pipe = psf.ducts.Pipe(path, diameter=0.025, roughness=15e-6)
    # tutorial:end:spline

    # tutorial:start:helix
    def helix(t):
        return np.column_stack((0.1 * np.cos(t), 0.1 * np.sin(t), 0.02 * t))

    helix_path = psf.ducts.Path.from_curve(helix, t_span=(0.0, 4 * np.pi))
    # tutorial:end:helix
    print(f"Helix: {len(helix_path.stations)} stations, length {helix_path.s[-1]:.4f} m")

    plot_spline_curvature(spline, spline_pipe).write_html(
        output_dir / "spline-curvature.html", include_plotlyjs="cdn")

    scene = psf.viz.Scene()
    scene.add(psf.viz.PipeComponent(spline_pipe), name="spline")
    scene.save_html(output_dir / "spline-pipe.html")

    # tutorial:start:solve
    fluid = psf.thermo.Fluid("fuel", ["Water"], [1.0])
    transport = psf.thermo.FluidTransport(fluid)
    resolver = psf.thermo.FluidStateResolver(transport, p_ref=2.0e6)
    inlet = resolver.from_tp(300.0, 2.0e6)

    rounded_result = rounded_pipe.solve(
        inlet_state=inlet,
        mass_flow=2.0,
        transport=transport,
        friction_model="fluids",
    )

    print(rounded_result.total_pressure_loss)
    print(rounded_result.states[-1])
    # tutorial:end:solve

    # tutorial:start:continuous
    nodes = [
        psf.ducts.PathNode((0.0, 0.0, 0.0)),
        psf.ducts.PathNode((0.0, 0.0, 0.5), bend=psf.ducts.Continuous(radius=0.25)),
        psf.ducts.PathNode((0.0, 0.5, 0.75)),
    ]
    path = psf.ducts.Path.from_nodes(nodes)
    continuous_pipe = psf.ducts.Pipe(path, diameter=0.025, roughness=15e-6)

    continuous_result = continuous_pipe.solve(
        inlet_state=inlet,
        mass_flow=2.0,
        transport=transport,
        friction_model="fluids",
    )
    # tutorial:end:continuous

    # tutorial:start:plot-loss
    plot = psf.viz.PlotPipePressureLoss(
        rounded_result,
        continuous_result,
        labels=["RoundedBend", "Continuous"],
    )
    plot.save_html(output_dir / "pressure-loss-comparison.html")
    # tutorial:end:plot-loss

    for name, result in (("RoundedBend", rounded_result),
                         ("Continuous", continuous_result)):
        bend = [i for i, c in enumerate(result.connections)
                if not isinstance(c, psf.ducts.Straight)]
        print(f"{name}: total {result.total_pressure_loss:.0f} Pa, "
              f"bend {result.pressure_loss[bend].sum():.0f} Pa, "
              f"bend K {result.loss_coefficient[bend].sum():.3f}")

    # tutorial:start:plot-pressure
    flared_result = flared_pipe.solve(
        inlet_state=inlet,
        mass_flow=2.0,
        transport=transport,
        friction_model="fluids",
    )

    plot = psf.viz.PlotPipePressure(flared_result, labels=["Flared pipe"])
    plot.save_html(output_dir / "flared-pressure.html")
    # tutorial:end:plot-pressure

    print(f"Flared: {len(flared_pipe.stations)} stations, "
          f"loss {flared_result.total_pressure_loss:.0f} Pa, "
          f"stagnation {flared_result.p_stagnation[[0, -1]]} Pa, "
          f"static {flared_result.p_static[[0, -1]]} Pa")
    print(f"Tapered: {len(tapered_pipe.stations)} stations")
    for name, pipe in (("Estimated", manual_pipe), ("Given", known_pipe)):
        result = pipe.solve(inlet_state=inlet, mass_flow=2.0, transport=transport)
        print(f"Hand-placed, radius {name.lower()}: radii {np.round(pipe.path.curvature_radius, 3)}, "
              f"loss {result.total_pressure_loss:.1f} Pa, length {pipe.path.s[-1]:.5f} m")

    # tutorial:start:network
    block = psf.network.PipeBlock(
        "fuel_feed",
        "tank_out",
        "pump_in",
        rounded_pipe,
        transport,
    )
    # tutorial:end:network

    # tutorial:start:scene-pressure
    scene = psf.viz.Scene(color_by="p_static")
    scene.add(psf.viz.PipeComponent(flared_pipe, flared_result), name="flared_bend")
    scene.save_html(output_dir / "flared-pressure-3d.html")
    # tutorial:end:scene-pressure

    # How finely the spline is followed: the default 5 degrees per connection
    # against much finer limits.
    for degrees in (5.0, 1.0, 0.25):
        sampled = psf.ducts.Pipe(
            psf.ducts.Path.from_curve(spline, (t[0], t[-1]), max_turn=np.radians(degrees)),
            diameter=0.025, roughness=15e-6)
        result = sampled.solve(inlet_state=inlet, mass_flow=2.0,
                               transport=transport, friction_model="fluids")
        print(f"Spline, max_turn {degrees} deg: {len(sampled.stations)} stations, "
              f"loss {result.total_pressure_loss:.1f} Pa, length {sampled.path.s[-1]:.5f} m")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Directory for generated HTML outputs.",
    )

    args = parser.parse_args()
    main(output_dir=args.output_dir)
