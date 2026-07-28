"""Mass-flow Sankey diagrams for converged PyskyFire engine networks.

The diagram uses fluid stations as Sankey nodes and fluid blocks as the links
between them.  This gives splitters and mergers their expected Sankey meaning:

* a splitter's branch-link widths are the outlet-station mass flows;
* a merger's branch-link widths are the inlet-station mass flows; and
* an ordinary one-in/one-out block uses its outlet-station mass flow.

The only required input is a converged :class:`~pyskyfire.common.EngineNetwork`.
"""

from __future__ import annotations

from collections import defaultdict
import colorsys
import hashlib
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
import plotly.graph_objects as go

from .core import PlotBase


@dataclass(frozen=True)
class _FlowLink:
    """One directed mass-flow connection in the rendered Sankey graph."""

    source: str
    target: str
    mdot: float
    block_name: str
    role: str


def _flatten_station_keys(value: Any) -> list[str]:
    """Return station keys from scalar or nested-list block metadata."""
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    if isinstance(value, Iterable):
        keys: list[str] = []
        for item in value:
            keys.extend(_flatten_station_keys(item))
        return keys
    return [str(value)]


def _finite_positive(value: Any, min_mdot: float) -> float | None:
    """Return a usable mass flow, or ``None`` when it cannot form a Sankey link."""
    try:
        mdot = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(mdot) or mdot <= min_mdot:
        return None
    return mdot


def _format_fraction(value: Any) -> str:
    """Return a compact, deterministic representation of a mixture fraction."""
    try:
        fraction = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not np.isfinite(fraction):
        return str(value)
    return f"{fraction:.6g}"


def _medium_identity(medium: Any) -> str:
    """Return a stable identity string for a pure fluid or fluid mixture.

    Strings are treated as pure-fluid identifiers.  Objects exposing
    ``propellants`` and ``fractions`` are canonicalised as sorted
    species/fraction pairs, so the same mixture receives the same colour even
    when components were supplied in a different order.
    """
    if medium is None:
        return "Unspecified"

    propellants = getattr(medium, "propellants", None)
    fractions = getattr(medium, "fractions", None)
    if propellants is not None and fractions is not None:
        try:
            pairs = sorted(
                (
                    (str(propellant).strip(), _format_fraction(fraction))
                    for propellant, fraction in zip(propellants, fractions)
                ),
                key=lambda pair: pair[0].casefold(),
            )
        except TypeError:
            pairs = []

        if pairs:
            return " + ".join(
                f"{species} ({fraction})"
                for species, fraction in pairs
            )

    return str(medium).strip() or "Unspecified"


def _medium_color(medium_identity: str) -> str:
    """Return a deterministic, medium-specific colour.

    The digest gives each canonical medium identity a stable pseudo-random hue.
    Saturation and value are also varied slightly so common media remain
    visually distinguishable without requiring a user-managed colour table.
    """
    digest = hashlib.blake2b(
        medium_identity.casefold().encode("utf-8"),
        digest_size=8,
    ).digest()

    hue = int.from_bytes(digest[:2], "big") / 65536.0
    saturation = 0.55 + digest[2] / 255.0 * 0.17
    value = 0.76 + digest[3] / 255.0 * 0.14
    red, green, blue = colorsys.hsv_to_rgb(hue, saturation, value)

    return "#{:02x}{:02x}{:02x}".format(
        round(255 * red),
        round(255 * green),
        round(255 * blue),
    )


def _rgba(hex_color: str, alpha: float) -> str:
    """Convert a ``#RRGGBB`` colour to Plotly-compatible rgba text."""
    if not 0.0 <= alpha <= 1.0:
        raise ValueError("alpha must be in the interval [0, 1]")

    color = hex_color.lstrip("#")
    if len(color) != 6:
        raise ValueError(f"Expected a #RRGGBB colour, got {hex_color!r}")

    red, green, blue = (
        int(color[0:2], 16),
        int(color[2:4], 16),
        int(color[4:6], 16),
    )
    return f"rgba({red}, {green}, {blue}, {alpha:.3g})"


def _darken(hex_color: str, factor: float = 0.72) -> str:
    """Return a darker version of a ``#RRGGBB`` colour for node borders."""
    if not 0.0 < factor <= 1.0:
        raise ValueError("factor must be in the interval (0, 1]")

    color = hex_color.lstrip("#")
    if len(color) != 6:
        raise ValueError(f"Expected a #RRGGBB colour, got {hex_color!r}")

    red, green, blue = (
        int(color[0:2], 16),
        int(color[2:4], 16),
        int(color[4:6], 16),
    )
    return "#{:02x}{:02x}{:02x}".format(
        round(red * factor),
        round(green * factor),
        round(blue * factor),
    )


def _pretty_station_name(name: str) -> str:
    """Make underscore-separated station keys easier to read without changing them."""
    return str(name).replace("_", " ")


def _format_station_hover(name: str, station: Any, medium: str) -> str:
    """Build a compact HTML hover label for a station node."""
    def number(attribute: str, unit: str, fmt: str) -> str:
        value = getattr(station, attribute, None)
        try:
            value = float(value)
        except (TypeError, ValueError):
            return f"{attribute} = unavailable"
        if not np.isfinite(value):
            return f"{attribute} = unavailable"
        return f"{attribute} = {value:{fmt}} {unit}"

    return "<br>".join(
        [
            f"<b>{_pretty_station_name(name)}</b>",
            f"Medium: {medium}",
            number("p", "Pa", ".4g"),
            number("T", "K", ".4g"),
            number("mdot", "kg/s", ".5g"),
        ]
    ) + "<extra></extra>"


class PlotMassFlowSankey(PlotBase):
    """Render a converged engine network as a mass-flow Sankey diagram.

    Parameters
    ----------
    engine_network:
        A converged :class:`~pyskyfire.common.EngineNetwork`.  It must expose
        ``stations`` and ``blocks``.  The station map supplies the converged
        mass flows; block ``station_inputs`` / ``station_outputs`` metadata
        supplies connectivity.
    title:
        Figure title.  Set to ``None`` to omit it.
    template:
        Plotly template name.
    min_mdot:
        Positive links at or below this threshold [kg/s] are omitted.  Their
        station nodes remain visible, which helps reveal unconnected or
        effectively closed branches.
    node_thickness, node_pad:
        Plotly Sankey node geometry in pixels.
    orientation:
        ``"v"`` (default) draws the engine flow from top to bottom.  This gives
        Plotly native horizontal node bars and horizontal station labels, which
        remain attached during every Plotly redraw.  Use ``"h"`` for the
        conventional left-to-right Sankey orientation.
    arrangement:
        Plotly node-placement mode.  ``"freeform"`` (default) lets the user drag
        stations anywhere in the Sankey canvas after the automatic layout has
        been generated.  Use ``"perpendicular"`` to restrict dragging across
        the flow direction, ``"snap"`` for Plotly-managed spacing, or
        ``"fixed"`` to disable node dragging.
    height, width:
        Optional figure dimensions in pixels.

    Notes
    -----
    Stations are the Sankey nodes, while fluid blocks are the links.  The
    layout uses the first appearance of each station in the ordered block list
    to form the main flow progression.  Links that return to an earlier
    station are retained as Plotly cyclic links, allowing pump recirculation
    and other closed engine flows to be displayed.

    A literal ``station -> same station`` block cannot be represented directly
    by Plotly.  In that rare case the plot adds a small virtual return node so
    the physical loop remains visible without creating a broken self-link.

    ``"freeform"`` changes the browser-side figure only.  Regenerating the
    plot from Python recreates the automatic initial layout unless the edited
    node coordinates are later saved and supplied back to the plot.

    Station nodes and flow links are coloured by fluid medium.  A pure fluid
    such as ``CH4`` or ``O2`` receives a deterministic colour.  A mixture is
    keyed by its sorted species/fraction pairs, so for example
    ``O2 + H2O + CO2`` is distinct from pure ``O2`` while receiving the same
    colour in every generated Sankey diagram.
    """

    def __init__(
        self,
        engine_network: Any,
        *,
        title: str | None = "Engine Mass-Flow Sankey",
        template: str = "plotly_white",
        min_mdot: float = 1.0e-12,
        node_thickness: int = 20,
        node_pad: int = 18,
        orientation: Literal["h", "v"] = "v",
        arrangement: Literal["snap", "perpendicular", "freeform", "fixed"] = "freeform",
        height: int | None = 800,
        width: int | None = None,
    ):
        if min_mdot < 0.0:
            raise ValueError("min_mdot must be non-negative")
        if node_thickness <= 0 or node_pad < 0:
            raise ValueError("node_thickness must be positive and node_pad non-negative")
        if orientation not in {"h", "v"}:
            raise ValueError("orientation must be either 'h' or 'v'")
        if arrangement not in {"snap", "perpendicular", "freeform", "fixed"}:
            raise ValueError(
                "arrangement must be one of 'snap', 'perpendicular', "
                "'freeform', or 'fixed'"
            )
        if not hasattr(engine_network, "stations") or not hasattr(engine_network, "blocks"):
            raise TypeError(
                "engine_network must expose 'stations' and 'blocks'; "
                "pass the converged EngineNetwork, e.g. results['net']."
            )

        super().__init__(go.Figure())
        self.template(template)
        self.engine_network = engine_network
        self.min_mdot = float(min_mdot)
        self.orientation = orientation
        self.arrangement = arrangement

        stations = dict(engine_network.stations)
        blocks = list(engine_network.blocks)
        if not stations:
            raise ValueError("engine_network.stations is empty")

        station_media = self._station_media(stations, blocks)
        station_order = self._station_order(stations, blocks)
        node_labels = {name: _pretty_station_name(name) for name in stations}
        virtual_nodes: dict[str, tuple[str, str]] = {}
        links = self._build_links(
            stations=stations,
            blocks=blocks,
            virtual_nodes=virtual_nodes,
        )

        # Include virtual nodes only where a component has a literal self-link
        # or an otherwise ambiguous many-input/many-output mapping.
        for node_id, (label, medium) in virtual_nodes.items():
            node_labels[node_id] = label
            station_media[node_id] = medium
            station_order.append(node_id)

        node_names = list(station_order)
        node_index = {name: index for index, name in enumerate(node_names)}
        x, y = self._layout_nodes(
            node_names,
            links,
            station_order,
            station_media,
            orientation=orientation,
        )

        medium_colors = {
            medium: _medium_color(medium)
            for medium in dict.fromkeys(station_media[name] for name in node_names)
        }
        node_colors = [
            medium_colors[station_media[name]]
            for name in node_names
        ]
        node_line_colors = [
            _darken(medium_colors[station_media[name]])
            for name in node_names
        ]

        visible_links = [
            link
            for link in links
            if _finite_positive(link.mdot, self.min_mdot) is not None
        ]

        node_hover = []
        for name in node_names:
            if name in stations:
                node_hover.append(
                    _format_station_hover(name, stations[name], station_media[name])
                )
            else:
                node_hover.append(
                    f"<b>{node_labels[name]}</b><br>Virtual routing node"
                    f"<br>Medium: {station_media[name]}<extra></extra>"
                )

        self.fig.add_trace(
            go.Sankey(
                arrangement=arrangement,
                orientation=orientation,
                node=dict(
                    # Keep Plotly's own labels.  Unlike plot annotations, they
                    # are part of the Sankey trace and remain attached to the
                    # nodes during all interactive redraws.
                    label=[node_labels[name] for name in node_names],
                    x=x,
                    y=y,
                    pad=node_pad,
                    thickness=node_thickness,
                    color=node_colors,
                    customdata=node_hover,
                    hovertemplate="%{customdata}",
                    line=dict(color=node_line_colors, width=0.7),
                ),
                link=dict(
                    source=[node_index[link.source] for link in visible_links],
                    target=[node_index[link.target] for link in visible_links],
                    value=[link.mdot for link in visible_links],
                    color=[
                        _rgba(medium_colors[station_media[link.source]], 0.42)
                        for link in visible_links
                    ],
                    # The source and target labels are stored explicitly so
                    # link hover remains informative when native node labels
                    # are hidden in vertical-label mode.
                    customdata=[
                        [
                            link.block_name,
                            link.role,
                            node_labels[link.source],
                            node_labels[link.target],
                        ]
                        for link in visible_links
                    ],
                    hovertemplate=(
                        "<b>%{customdata[2]} → %{customdata[3]}</b>"
                        "<br>Block: %{customdata[0]}"
                        "<br>Role: %{customdata[1]}"
                        "<br>Mass flow: %{value:.5g} kg/s"
                        "<extra></extra>"
                    ),
                ),
            )
        )

        self.station_nodes = tuple(stations)
        self.flow_links = tuple(visible_links)
        self.medium_colors = dict(medium_colors)
        self.node_positions = {
            name: (x[index], y[index])
            for index, name in enumerate(node_names)
        }

        self.layout(
            title=title,
            height=height,
            width=width,
            font=dict(size=12),
            margin=dict(l=20, r=20, t=60 if title else 20, b=20),
        )

    @staticmethod
    def _station_media(stations: dict[str, Any], blocks: list[Any]) -> dict[str, str]:
        """Associate every station with the first fluid medium that touches it."""
        media: dict[str, str] = {}
        for block in blocks:
            block_medium = getattr(block, "medium", None)
            if block_medium is None:
                continue
            medium_name = _medium_identity(block_medium)
            keys = (
                _flatten_station_keys(getattr(block, "station_inputs", []))
                + _flatten_station_keys(getattr(block, "station_outputs", []))
            )
            for key in keys:
                if key in stations:
                    media.setdefault(key, medium_name)
        return {name: media.get(name, "Unspecified") for name in stations}

    @staticmethod
    def _station_order(stations: dict[str, Any], blocks: list[Any]) -> list[str]:
        """Return deterministic station order led by fluid-block declaration order."""
        order: list[str] = []
        seen: set[str] = set()
        for block in blocks:
            inputs = _flatten_station_keys(getattr(block, "station_inputs", []))
            outputs = _flatten_station_keys(getattr(block, "station_outputs", []))
            if not inputs and not outputs:
                continue
            for key in [*inputs, *outputs]:
                if key in stations and key not in seen:
                    order.append(key)
                    seen.add(key)
        for key in stations:
            if key not in seen:
                order.append(key)
        return order

    def _build_links(
        self,
        *,
        stations: dict[str, Any],
        blocks: list[Any],
        virtual_nodes: dict[str, tuple[str, str]],
    ) -> list[_FlowLink]:
        """Translate declared fluid-block topology into station-to-station links."""
        links: list[_FlowLink] = []
        virtual_index = 0

        def station_mdot(name: str) -> float | None:
            if name not in stations:
                raise KeyError(
                    f"Station {name!r} is declared by a block but is absent from "
                    "engine_network.stations. Pass a fully constructed, converged network."
                )
            return _finite_positive(getattr(stations[name], "mdot", None), self.min_mdot)

        def block_name(block: Any, index: int) -> str:
            return str(getattr(block, "name", f"{type(block).__name__}_{index + 1}"))

        def make_virtual(label: str, medium: str) -> str:
            nonlocal virtual_index
            virtual_index += 1
            node_id = f"__pyskyfire_sankey_virtual_{virtual_index}"
            virtual_nodes[node_id] = (label, medium)
            return node_id

        def add_link(source: str, target: str, mdot: float | None, name: str, role: str, medium: str) -> None:
            if mdot is None:
                return
            if source == target:
                # Plotly does not reliably render a literal self-link.  Retain
                # the mass flow by routing through one virtual loop node instead.
                loop_node = make_virtual(f"{name} (return)", medium)
                links.append(_FlowLink(source, loop_node, mdot, name, role))
                links.append(_FlowLink(loop_node, target, mdot, name, role))
            else:
                links.append(_FlowLink(source, target, mdot, name, role))

        for index, block in enumerate(blocks):
            inputs = _flatten_station_keys(getattr(block, "station_inputs", []))
            outputs = _flatten_station_keys(getattr(block, "station_outputs", []))
            if not inputs and not outputs:
                continue  # Signal-only block.
            if not inputs or not outputs:
                # A fluid link requires a transformation from at least one
                # declared station to at least one declared station.
                continue

            name = block_name(block, index)
            medium = _medium_identity(getattr(block, "medium", None))

            if len(inputs) == 1:
                # One inlet, potentially multiple splitter outlets.  The
                # output station determines the branch mass flow.
                source = inputs[0]
                for target in outputs:
                    add_link(
                        source,
                        target,
                        station_mdot(target),
                        name,
                        "outlet branch" if len(outputs) > 1 else "throughflow",
                        medium,
                    )
            elif len(outputs) == 1:
                # Multiple merger inlets, one outlet.  Each input station
                # carries its own branch mass flow.
                target = outputs[0]
                for source in inputs:
                    add_link(
                        source,
                        target,
                        station_mdot(source),
                        name,
                        "inlet branch" if len(inputs) > 1 else "throughflow",
                        medium,
                    )
            else:
                # A general many-to-many fluid component has no unambiguous
                # station-to-station allocation in the current block metadata.
                # Preserve all known mass flows with one labelled virtual node.
                component_node = make_virtual(name, medium)
                for source in inputs:
                    add_link(
                        source,
                        component_node,
                        station_mdot(source),
                        name,
                        "inlet branch",
                        medium,
                    )
                for target in outputs:
                    add_link(
                        component_node,
                        target,
                        station_mdot(target),
                        name,
                        "outlet branch",
                        medium,
                    )

        return links

    @staticmethod
    def _layout_nodes(
        node_names: list[str],
        links: list[_FlowLink],
        node_order: list[str],
        node_media: dict[str, str],
        *,
        orientation: Literal["h", "v"],
    ) -> tuple[list[float], list[float]]:
        """Assign deterministic node positions while retaining cycles.

        Stations are initially ordered by first appearance in the block list.
        A link to an earlier station becomes a feedback link.  Removing only
        these feedback links leaves an acyclic graph for progression assignment;
        the original feedback links are still passed to Plotly and are rendered
        as recirculation loops.

        The internal coordinate layout always uses progress on x and fluid
        lanes on y.  Plotly rotates that internal layout when
        ``orientation="v"`` is selected, producing a top-to-bottom flow with
        side-by-side fluid lanes on screen.
        """
        rank = {name: index for index, name in enumerate(node_order)}
        forward_links = [
            link for link in links
            if rank[link.source] < rank[link.target]
        ]

        successors: dict[str, list[str]] = defaultdict(list)
        predecessors: dict[str, list[str]] = defaultdict(list)
        for link in forward_links:
            successors[link.source].append(link.target)
            predecessors[link.target].append(link.source)

        layers = {name: 0 for name in node_names}
        for source in node_order:
            for target in successors[source]:
                layers[target] = max(layers[target], layers[source] + 1)

        max_layer = max(layers.values(), default=0)
        if max_layer == 0:
            x_by_name = {name: 0.5 for name in node_names}
        else:
            x_by_name = {
                name: 0.03 + 0.94 * layers[name] / max_layer
                for name in node_names
            }

        media_order: list[str] = []
        for name in node_order:
            medium = node_media.get(name, "Unspecified")
            if medium not in media_order:
                media_order.append(medium)
        if not media_order:
            media_order = ["Unspecified"]

        nodes_by_medium: dict[str, list[str]] = defaultdict(list)
        for name in node_names:
            nodes_by_medium[node_media.get(name, "Unspecified")].append(name)

        # Allocate a vertical band to each fluid.  This keeps, for example,
        # fuel and oxidizer circuits distinct while preserving common columns.
        band_gap = 0.055 if len(media_order) > 1 else 0.0
        usable_height = 1.0 - band_gap * (len(media_order) - 1)
        band_height = usable_height / len(media_order)
        y_by_name: dict[str, float] = {}

        for medium_index, medium in enumerate(media_order):
            medium_nodes = nodes_by_medium[medium]
            lower = medium_index * (band_height + band_gap)
            upper = lower + band_height

            nodes_by_layer: dict[int, list[str]] = defaultdict(list)
            for name in medium_nodes:
                nodes_by_layer[layers[name]].append(name)
            for layer_nodes in nodes_by_layer.values():
                layer_nodes.sort(key=rank.__getitem__)

            # A few forward/backward barycentre passes reduce gratuitous link
            # crossings without requiring networkx or a separate layout engine.
            normalized_y = {
                name: (index + 1) / (len(medium_nodes) + 1)
                for index, name in enumerate(sorted(medium_nodes, key=rank.__getitem__))
            }
            medium_forward = [
                link for link in forward_links
                if node_media.get(link.source) == medium
                and node_media.get(link.target) == medium
            ]
            local_predecessors: dict[str, list[str]] = defaultdict(list)
            local_successors: dict[str, list[str]] = defaultdict(list)
            for link in medium_forward:
                local_successors[link.source].append(link.target)
                local_predecessors[link.target].append(link.source)

            for _ in range(3):
                for layer in range(1, max_layer + 1):
                    layer_nodes = nodes_by_layer.get(layer, [])
                    layer_nodes.sort(
                        key=lambda name: (
                            float(np.mean([normalized_y[p] for p in local_predecessors[name]]))
                            if local_predecessors[name]
                            else normalized_y[name],
                            rank[name],
                        )
                    )
                    for index, name in enumerate(layer_nodes):
                        normalized_y[name] = (index + 1) / (len(layer_nodes) + 1)

                for layer in range(max_layer - 1, -1, -1):
                    layer_nodes = nodes_by_layer.get(layer, [])
                    layer_nodes.sort(
                        key=lambda name: (
                            float(np.mean([normalized_y[s] for s in local_successors[name]]))
                            if local_successors[name]
                            else normalized_y[name],
                            rank[name],
                        )
                    )
                    for index, name in enumerate(layer_nodes):
                        normalized_y[name] = (index + 1) / (len(layer_nodes) + 1)

            # Leave a little room at the top and bottom so high-width return
            # loops do not sit exactly against the figure border.
            for name in medium_nodes:
                y_by_name[name] = lower + band_height * (
                    0.06 + 0.88 * normalized_y[name]
                )

        # ``node.x`` and ``node.y`` are specified in Plotly's *internal*
        # Sankey coordinates.  For ``orientation="v"``, Plotly rotates that
        # internal canvas when rendering, so passing swapped coordinates here
        # would rotate the placement twice.  Keep the same internal layout for
        # both modes: progression on x and medium lanes on y.  Plotly performs
        # the final rotation for the vertical presentation.
        return (
            [x_by_name[name] for name in node_names],
            [y_by_name[name] for name in node_names],
        )
