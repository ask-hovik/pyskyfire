"""Injector scenes and face views.

:class:`InjectorComponent` draws the wetted surface of every hole, its bore
behind the face, one colour per feed, optionally carrying each hole's flow from
an :class:`~pyskyfire.injector.InjectorResult` as scene fields. Arrows from
each hole to its impingement point can be added; they are annotation, not
hardware.
:func:`plot_injector_face` is the same injector seen from the chamber: hole
outlines, the jets projected onto the face, and where they meet.

The injector frame is the scene frame: the face is the plane x = 0 and the
chamber lies towards +x, as for a thrust chamber whose contour starts at the
face.
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import math

import numpy as np
import plotly.graph_objects as go

from .components import _carrying, _fan, _merge, _surface
from .scene import Anchor, FieldData, SceneNode

__all__ = ["InjectorComponent", "plot_injector_face", "FEED_COLORS"]

#: Categorical slots for feeds, in the injector's feed order. Validated for
#: colour-vision deficiency against a light surface; identity is also carried
#: by legends and node names.
FEED_COLORS = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300")


def _feed_colors(feeds, colors):
    if colors is None:
        if len(feeds) > len(FEED_COLORS):
            raise ValueError(f"more than {len(FEED_COLORS)} feeds need colors given explicitly")
        return {feed: FEED_COLORS[i] for i, feed in enumerate(feeds)}
    missing = [feed for feed in feeds if feed not in colors]
    if missing:
        raise ValueError(f"colors has no entry for feeds {missing}")
    return dict(colors)


def _feed_labels(feeds, labels):
    labels = {} if labels is None else dict(labels)
    return {feed: labels.get(feed, feed if isinstance(feed, str) else f"feed {i + 1}")
            for i, feed in enumerate(feeds)}


def _diameter(orifice):
    if orifice.hydraulic_diameter is not None:
        return orifice.hydraulic_diameter
    return math.sqrt(4.0 * orifice.area / math.pi)


def _ring(axis, radius, n):
    """``n`` points on a circle of ``radius`` about the origin, normal to ``axis``."""
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    helper = np.eye(3)[np.argmin(np.abs(axis))]
    u = np.cross(axis, helper)
    u /= np.linalg.norm(u)
    v = np.cross(axis, u)
    angle = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    return radius * (np.cos(angle)[:, None] * u + np.sin(angle)[:, None] * v)


def _bore(exit_, direction, length, radius, n):
    """A hole's wetted wall, from ``length`` upstream of its exit to the face.

    The upstream end is square to the axis. The downstream end is where the
    bore meets the face plane x = 0, an ellipse for an inclined hole; no
    counterbore or chamfer is drawn.
    """
    ring = _ring(direction, radius, n)
    inlet = exit_ - length * direction + ring
    at_exit = exit_ + ring
    face = at_exit - (at_exit[:, 0] / direction[0])[:, None] * direction
    return _surface(np.stack((inlet, face)), closed=True)


def _arrow(start, end, radius, n):
    """A shaft and conical head from ``start`` to ``end``."""
    axis = np.asarray(end, float) - np.asarray(start, float)
    length = np.linalg.norm(axis)
    head = min(0.3 * length, 6.0 * radius)
    base = end - head * axis / length
    shaft = _surface(np.stack((start + _ring(axis, radius, n), base + _ring(axis, radius, n))),
                     closed=True)
    cone = _ring(axis, 2.5 * radius, n)
    return _merge([shaft, _fan(base + cone, end, outward=True),
                   _fan(base + cone, base, outward=False)])


def _disc(radius, n):
    angle = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    rim = np.column_stack((np.zeros(n), radius * np.cos(angle), radius * np.sin(angle)))
    return _surface(np.stack((np.zeros_like(rim), rim)), closed=True)


def _hole_flows(result, holes):
    """Mass flow and jet velocity of each hole from an InjectorResult."""
    mass_flow = np.array([result.feeds[f].hole_mass_flow(o)
                          for f, o in zip(holes.feeds, holes.orifices)])
    velocity = np.array([result.feeds[f].velocity for f in holes.feeds])
    return mass_flow, velocity


@dataclass
class InjectorComponent:
    """Display the wetted bores of an injector's holes, coloured by feed.

    Parameters
    ----------
    injector : pyskyfire.injector.Injector
        The injector.
    result : pyskyfire.injector.InjectorResult, optional
        An operating point. Each hole then carries ``hole_mass_flow`` [kg/s]
        and ``jet_velocity`` [m/s] as point fields, and each feed's bores
        record its CoolProp fluid name as ``metadata["fluid"]``.
    show_impingement : bool, optional
        Add an arrow from each hole's exit to its target, the impingement
        point for an impinging element.
    face_radius : float, optional
        Draw a translucent face of this radius [m], such as the chamber
        radius. The injector itself has no face outline.
    colors, labels : mapping of feed to str, optional
        Colour and display name of each feed. Names default to the feed
        itself when it is a string.
    n_theta : int, optional
        Facets round each bore and jet.

    Notes
    -----
    Only wetted surface is drawn as hardware. A hole's bore runs back from the
    face along its axis for the orifice's ``length``, or two diameters where
    no length is given, and ends where it meets the face plane. The arrows
    are drawn in grey, without fields, as annotation. The ``face`` anchor sits
    on the axis at the face, pointing into the chamber.
    """

    injector: object
    result: object = None
    show_impingement: bool = False
    face_radius: float | None = None
    colors: Mapping | None = None
    labels: Mapping | None = None
    n_theta: int = 16

    def build(self):
        if not isinstance(self.n_theta, int) or self.n_theta < 3:
            raise ValueError("n_theta must be an integer of at least 3")
        feeds = self.injector.feeds
        colors = _feed_colors(feeds, self.colors)
        labels = _feed_labels(feeds, self.labels)
        holes = self.injector.holes()
        flows = None if self.result is None else _hole_flows(self.result, holes)
        n = self.n_theta

        root = SceneNode("injector", metadata={"kind": "injector"})
        if self.face_radius is not None:
            root.children.append(SceneNode("face", mesh=_disc(float(self.face_radius), 4*n),
                                           color="silver", opacity=0.25,
                                           metadata={"kind": "injector_face"}))
        for s in self.injector.element_sets:
            group = SceneNode(s.name, metadata={"kind": "element_set", "count": s.count})
            for feed in s.element.feeds:
                index = [i for i, (name, f) in enumerate(zip(holes.set_name, holes.feeds))
                         if name == s.name and f == feed]
                bores, arrows = [], []
                for i in index:
                    orifice = holes.orifices[i]
                    d = _diameter(orifice)
                    exit_, direction = holes.exit[i], holes.direction[i]
                    length = orifice.length if orifice.length > 0.0 else 2.0 * d
                    bores.append(_bore(exit_, direction, length, 0.5 * d, n))
                    arrows.append(_arrow(exit_, holes.target[i], 0.1 * d, n))
                fields = {}
                if flows is not None:
                    per_point = 2 * n
                    fields = {
                        "hole_mass_flow": FieldData(np.repeat(flows[0][index], per_point),
                                                    "kg/s", "Hole mass flow"),
                        "jet_velocity": FieldData(np.repeat(flows[1][index], per_point),
                                                  "m/s", "Jet velocity"),
                    }
                meta = {"kind": "injector_holes", "feed": labels[feed], "holes": len(index)}
                fluid = None if self.result is None else self.result.feeds[feed].fluid
                group.children.append(SceneNode(f"{labels[feed]} bores", mesh=_merge(bores),
                                                fields=dict(fields), color=colors[feed],
                                                metadata=_carrying(dict(meta), fluid)))
                if self.show_impingement:
                    group.children.append(SceneNode(
                        f"{labels[feed]} impingement", mesh=_merge(arrows), color="dimgray",
                        metadata={**meta, "kind": "impingement_vectors"}))
            root.children.append(group)
        root.anchors = {"face": Anchor.from_direction((0.0, 0.0, 0.0), (1.0, 0.0, 0.0))}
        return root


def plot_injector_face(injector, *, face_radius=None, colors=None, labels=None, title=None):
    """The injector seen from the chamber, in millimetres.

    Each hole's outline on the face is drawn in its feed's colour, with its
    jet projected onto the face as a line to the impingement point, where a
    grey dot marks where the jets meet. ``face_radius`` adds the chamber wall
    as a grey circle. Hovering a hole shows its set, feed and diameter.

    Parameters are as for :class:`InjectorComponent`.
    """
    feeds = injector.feeds
    colors = _feed_colors(feeds, colors)
    labels = _feed_labels(feeds, labels)
    holes = injector.holes()
    y, z = holes.exit[:, 1] * 1e3, holes.exit[:, 2] * 1e3
    circle = np.linspace(0.0, 2.0 * np.pi, 33)
    ink, grid = "#3d3d3a", "#e4e4e0"

    fig = go.Figure()
    if face_radius is not None:
        r = float(face_radius) * 1e3
        wall = np.linspace(0.0, 2.0 * np.pi, 257)
        fig.add_trace(go.Scatter(x=r*np.cos(wall), y=r*np.sin(wall), mode="lines",
                                 line=dict(color="#8a8a86", width=1), name="Chamber wall",
                                 hoverinfo="skip"))
    for feed in feeds:
        index = [i for i, f in enumerate(holes.feeds) if f == feed]
        outline_x, outline_y, jet_x, jet_y, text = [], [], [], [], []
        for i in index:
            radius = holes.footprint[i] * 1e3
            outline_x += list(y[i] + radius*np.cos(circle)) + [None]
            outline_y += list(z[i] + radius*np.sin(circle)) + [None]
            jet_x += [y[i], holes.target[i, 1]*1e3, None]
            jet_y += [z[i], holes.target[i, 2]*1e3, None]
            text.append(f"{holes.set_name[i]}, element {holes.element_index[i] + 1}<br>"
                        f"{labels[feed]}<br>d = {_diameter(holes.orifices[i])*1e3:.3g} mm")
        group = str(labels[feed])
        fig.add_trace(go.Scatter(x=jet_x, y=jet_y, mode="lines", legendgroup=group,
                                 line=dict(color=colors[feed], width=1), opacity=0.6,
                                 showlegend=False, hoverinfo="skip"))
        fig.add_trace(go.Scatter(x=outline_x, y=outline_y, mode="lines", legendgroup=group,
                                 line=dict(color=colors[feed], width=2), name=group,
                                 hoverinfo="skip"))
        # Invisible hit targets, larger than the smallest holes.
        fig.add_trace(go.Scatter(x=y[index], y=z[index], mode="markers", legendgroup=group,
                                 marker=dict(size=14, color="rgba(0,0,0,0)"), text=text,
                                 hovertemplate="%{text}<extra></extra>", showlegend=False))
    targets = np.unique(np.round(holes.target[:, 1:] * 1e3, 9), axis=0)
    fig.add_trace(go.Scatter(x=targets[:, 0], y=targets[:, 1], mode="markers",
                             marker=dict(size=4, color="#8a8a86"), name="Impingement points",
                             hoverinfo="skip"))
    fig.update_xaxes(title_text="y [mm]", gridcolor=grid, zeroline=False, color=ink)
    fig.update_yaxes(title_text="z [mm]", gridcolor=grid, zeroline=False, color=ink,
                     scaleanchor="x", scaleratio=1)
    fig.update_layout(title=title or "Injector face, seen from the chamber",
                      template="plotly_white", hovermode="closest",
                      legend=dict(orientation="h", y=-0.15), height=720)
    return fig
