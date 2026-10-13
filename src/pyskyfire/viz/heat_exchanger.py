"""Heat-exchanger swept geometry and wall/flow report panels."""

import numpy as np
import plotly.graph_objects as go
from .report import Report

__all__ = ["plot_heat_exchanger_geometry", "heat_exchanger_report"]


def plot_heat_exchanger_geometry(exchanger, *, installation=None):
    """Draw the actual exposed tube envelopes, including transported sections."""
    fig = go.Figure()
    colors = ["#e58b27", "#327fa8", "#73a857", "#a665a8", "#df615c"]
    for color, c in zip(
        (colors[i % len(colors)] for i in range(len(exchanger.circuits))),
        exchanger.circuits.values(),
    ):
        vertices = []
        triangles = []
        for j in range(c.n_channels):
            points = c.section_points(j, samples=8)
            if installation is not None:
                points = (
                    points @ installation.gas.rotation.T + installation.gas.position
                )
            n, m, _ = points.shape
            offset = len(vertices)
            vertices.extend(points.reshape(-1, 3))
            for i in range(n - 1):
                for k in range(m):
                    a = offset + i * m + k
                    b = offset + i * m + (k + 1) % m
                    d = offset + (i + 1) * m + k
                    e = offset + (i + 1) * m + (k + 1) % m
                    triangles.extend(((a, b, e), (a, e, d)))
        v = np.asarray(vertices)
        t = np.asarray(triangles)
        fig.add_trace(
            go.Mesh3d(
                x=v[:, 0],
                y=v[:, 1],
                z=v[:, 2],
                i=t[:, 0],
                j=t[:, 1],
                k=t[:, 2],
                color=color,
                name=c.name,
                showlegend=True,
                opacity=1.0,
            )
        )
        points = c.placed_paths[0].centerline
        if installation is not None:
            points = points @ installation.gas.rotation.T + installation.gas.position
        fig.add_trace(
            go.Scatter3d(
                x=points[[0, -1], 0],
                y=points[[0, -1], 1],
                z=points[[0, -1], 2],
                mode="markers",
                text=["inlet", "outlet"],
                textposition="top center",
                marker=dict(size=3, color=color),
                showlegend=False,
                name=c.name,
            )
        )
    fig.update_layout(
        title="Heat-exchanger channel envelopes",
        scene=dict(
            aspectmode="data",
            xaxis_title="Axial x [m]",
            yaxis_title="y [m]",
            zaxis_title="z [m]",
        ),
        height=720,
    )
    return fig


def heat_exchanger_report(
    exchanger, result, *, title="Heat exchanger", installation=None, scene=None
):
    """Create a Report with geometry, fluid states, wall temperatures, and balance."""
    report = Report(title)
    tab = report.add_tab("Summary")
    tab.add_table(
        dict(
            iterations=result.iterations,
            relative_coupling_residual=result.relative_residual,
            energy_balance_error_W=result.energy_balance_error,
            prescribed_heat_W=result.prescribed_heat,
            displaced_volume_m3=result.displaced_volume,
        )
    )
    if installation is not None:
        tab.add_table(vars(installation.volume_change), caption="Volume change [m³]")
        tab.add_markdown(
            "Gas temperature and chemistry are prescribed. Volume is summed over nonoverlapping channel envelopes, including the coolant passages. Ideal series headers are not included in the swept volume."
        )
    geom = report.add_tab("Geometry")
    if scene is None:
        geom.add_figure(
            plot_heat_exchanger_geometry(exchanger, installation=installation)
        )
    else:
        geom.add_iframe(scene.data_url, caption="Chamber and heat exchanger")
    for name, r in result.circuits.items():
        tab = report.add_tab(name)
        c = r.circuit
        wall_values = [
            T
            for w in result.walls
            for T, owner in (
                (w.temperature_a, w.circuit_a),
                (w.temperature_b, w.circuit_b),
            )
            if owner == name
        ]
        tab.add_table(
            dict(
                copies=c.n_channels,
                total_mass_flow_kg_s=r.mass_flow,
                per_channel_mass_flow_kg_s=r.mass_flow / c.n_channels,
                inlet_total_area_m2=c.area[0] * c.n_channels,
                inlet_temperature_K=r.T[0],
                outlet_temperature_K=r.T[-1],
                static_pressure_drop_Pa=r.p[0] - r.p[-1],
                heat_to_fluid_W=r.total_heat,
                max_wall_temperature_K=max(wall_values) if wall_values else None,
            )
        )
        for values, label in (
            (r.T, "Fluid temperature [K]"),
            (r.p / 1e5, "Static pressure [bar]"),
            (r.velocity, "Velocity [m/s]"),
        ):
            fig = go.Figure(go.Scatter(x=r.s, y=values, mode="lines", name=name))
            fig.update_layout(
                xaxis_title="Distance along flow path [m]",
                yaxis_title=label,
                title=label,
            )
            tab.add_figure(fig)
        fig = go.Figure()
        for surface in c.cross_section.surface_names:
            inner = []
            outer = []
            for w in result.walls:
                if w.circuit_a == name and w.surface_a == surface:
                    s = (c.s[w.segment_a] + c.s[w.segment_a + 1]) / 2
                    inner.append((s, w.temperature_a))
                    outer.append((s, w.temperature_b))
                if w.circuit_b == name and w.surface_b == surface:
                    s = (c.s[w.segment_b] + c.s[w.segment_b + 1]) / 2
                    inner.append((s, w.temperature_b))
                    outer.append((s, w.temperature_a))
            for points, side in (
                (inner, "fluid interface"),
                (outer, "opposite interface"),
            ):
                if points:
                    x, y = np.asarray(sorted(points)).T
                    fig.add_trace(
                        go.Scatter(x=x, y=y, name=f"{surface}: {side}", mode="lines")
                    )
        fig.update_layout(
            title="Wall temperatures by exposed surface",
            xaxis_title="Distance along flow path [m]",
            yaxis_title="Wall temperature [K]",
        )
        tab.add_figure(fig)
    return report
