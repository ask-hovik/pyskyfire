"""Interactive HTML visualisation for :class:`EngineNetwork` objects.

The visualiser deliberately treats the network topology as model data and the
node positions / edge bend points as presentation data.  This makes it useful
both as a script output today and as the basis for a future graphical editor.
"""

from __future__ import annotations

import base64
import hashlib
import html
import json
import re
import webbrowser
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping
from importlib.resources import files


_X6_ASSET_NAME = "x6.min.js"


def _read_packaged_asset(filename: str) -> str:
    """Return a UTF-8 JavaScript asset bundled with pyskyfire.viz."""
    try:
        return (
            files("pyskyfire.viz")
            .joinpath("assets", filename)
            .read_text(encoding="utf-8")
        )
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"Missing bundled visualisation asset: {filename}. "
            "Expected it in pyskyfire/viz/assets/."
        ) from exc

_LAYOUT_VERSION = 1


_BLOCK_TITLES = {
    "PumpBlock": "Pump",
    "TurbineBlock": "Turbine",
    "RegenBlock": "Regen cooling",
    "SimpleDuctBlock": "Duct",
    "MassFlowSplitterBlock": "Flow splitter",
    "MassFlowMergerBlock": "Flow merger",
    "TransmissionBlock": "Shaft transmission",
}


def _flatten_keys(value: Any) -> list[str]:
    """Return metadata keys, accepting the nested lists used by older blocks."""
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    if isinstance(value, Iterable):
        keys: list[str] = []
        for item in value:
            keys.extend(_flatten_keys(item))
        return keys
    return [str(value)]


def _slug(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_-]+", "-", value).strip("-")
    return value or "item"


def _stable_id(prefix: str, *parts: Any) -> str:
    """Return an X6-safe ID stable across runs if network topology is unchanged."""
    source = "\x1f".join(str(part) for part in parts)
    digest = hashlib.sha1(source.encode("utf-8")).hexdigest()[:12]
    label = _slug(str(parts[0])) if parts else prefix
    return f"{prefix}:{label}:{digest}"


def _json_for_script(value: Any) -> str:
    """Serialize JSON safely for an inline ``<script>`` element."""
    return json.dumps(value, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")


def _read_layout(layout: Mapping[str, Any] | str | Path | None) -> dict[str, Any]:
    if layout is None:
        return {"version": _LAYOUT_VERSION, "nodes": {}, "edges": {}}

    if isinstance(layout, (str, Path)):
        layout_path = Path(layout)
        try:
            value = json.loads(layout_path.read_text(encoding="utf-8"))
        except OSError as exc:
            raise OSError(f"Could not read network layout from {layout_path}") from exc
        except json.JSONDecodeError as exc:
            raise ValueError(f"Network layout file is not valid JSON: {layout_path}") from exc
    elif isinstance(layout, Mapping):
        value = dict(layout)
    else:
        raise TypeError("layout must be a mapping, path, or None")

    if not isinstance(value, Mapping):
        raise ValueError("Network layout must be a JSON object")

    nodes = value.get("nodes", {})
    edges = value.get("edges", {})
    if not isinstance(nodes, Mapping) or not isinstance(edges, Mapping):
        raise ValueError("Network layout requires object-valued 'nodes' and 'edges' fields")

    return {
        "version": value.get("version", _LAYOUT_VERSION),
        "nodes": dict(nodes),
        "edges": dict(edges),
    }


def _block_title(block: Any) -> str:
    return _BLOCK_TITLES.get(type(block).__name__, type(block).__name__.removesuffix("Block"))


def _block_name(block: Any, index: int) -> str:
    return str(getattr(block, "name", f"{type(block).__name__}_{index + 1}"))


def _block_signal_inputs(block: Any) -> list[str]:
    """Return declared scalar inputs, excluding internal pressure-drop bookkeeping."""
    return [
        key
        for key in _flatten_keys(getattr(block, "signal_inputs", []))
        if not key.startswith("dp_")
    ]

def _block_colors(spec: _BlockSpec) -> tuple[str, str]:
    """Return (fill, stroke) colors for a block."""
    if spec.is_signal_only:
        return "#fff8e1", "#f9a825"   # yellow

    if spec.class_name == "RegenBlock":
        return "#fdecea", "#c62828"   # red
    if spec.class_name == "SimpleDuctBlock":
        return "#e3f2fd", "#1565c0"   # blue
    if spec.class_name == "PumpBlock":
        return "#e8f5e9", "#2e7d32"   # green
    if spec.class_name == "TurbineBlock":
        return "#f3e5f5", "#7b1fa2"   # purple

    return "#eef5fb", "#607d8b"       # default / unspecified


@dataclass(frozen=True)
class _BlockSpec:
    id: str
    class_name: str
    name: str
    title: str
    medium: str | None
    station_inputs: tuple[str, ...]
    station_outputs: tuple[str, ...]
    signal_inputs: tuple[str, ...]
    signal_outputs: tuple[str, ...]

    @property
    def is_signal_only(self) -> bool:
        return not self.station_inputs and not self.station_outputs

@dataclass
class NetworkVisualizer:
    """Build an editable X6 HTML view of an engine network.

    Parameters
    ----------
    network:
        An object exposing ``blocks`` with PyskyFire's station/signal metadata.
    title:
        Browser title and heading displayed above the canvas.
    layout:
        Optional layout mapping or JSON file previously downloaded from the
        viewer.  It contains only presentation data, never solver data.
    Notes
    -----
    The HTML visualiser is editable in the browser.  A browser cannot safely
    overwrite an arbitrary local file, so the toolbar downloads a sidecar
    layout JSON file.  Pass that file back through ``layout=`` on a later run.
    """

    network: Any
    title: str = "PyskyFire Engine Network"
    layout: Mapping[str, Any] | str | Path | None = None

    def __post_init__(self) -> None:
        if not hasattr(self.network, "blocks"):
            raise TypeError("network must expose a 'blocks' attribute")
        self._layout = _read_layout(self.layout)

    @property
    def data_url(self) -> str:
        """Base64 ``data:text/html`` URL suitable for :meth:`Tab.add_iframe`."""
        encoded = base64.b64encode(self.to_html().encode("utf-8")).decode("ascii")
        return f"data:text/html;base64,{encoded}"

    def save_html(self, path: str | Path) -> Path:
        """Write the viewer to an HTML file and return the resolved path."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.to_html(), encoding="utf-8")
        return path

    def show(self) -> None:
        """Open a temporary standalone visualisation in the default browser."""
        import tempfile

        with tempfile.NamedTemporaryFile(prefix="pyskyfire-network-", suffix=".html", delete=False) as handle:
            path = Path(handle.name)
        self.save_html(path)
        webbrowser.open(path.as_uri())

    def to_html(self) -> str:
        """Return a standalone HTML document containing the X6 viewer."""
        payload = self._make_payload()
        template = _HTML_TEMPLATE
        return (
            template.replace("__TITLE__", html.escape(self.title))
            .replace("__X6_SOURCE__", _read_packaged_asset(_X6_ASSET_NAME))
            .replace("__PAYLOAD__", _json_for_script(payload))
        )

    def _make_payload(self) -> dict[str, Any]:
        specs = self._block_specs()
        nodes = self._make_nodes(specs)
        edges, external_nodes = self._make_edges(specs, nodes)
        nodes.extend(external_nodes)
        return {
            "version": 1,
            "title": self.title,
            "nodes": nodes,
            "edges": edges,
            "layout": self._layout,
        }

    def _block_specs(self) -> list[_BlockSpec]:
        specs: list[_BlockSpec] = []
        for index, block in enumerate(self.network.blocks):
            station_inputs = tuple(_flatten_keys(getattr(block, "station_inputs", [])))
            station_outputs = tuple(_flatten_keys(getattr(block, "station_outputs", [])))
            signal_inputs = tuple(_block_signal_inputs(block))
            signal_outputs = tuple(
                key
                for key in _flatten_keys(getattr(block, "signal_outputs", []))
                if not key.startswith("dp_")
            )

            class_name = type(block).__name__
            name = _block_name(block, index)
            block_id = _stable_id(
                "block",
                class_name,
                name,
                ",".join(station_inputs),
                ",".join(station_outputs),
            )

            medium = getattr(block, "medium", None)
            specs.append(
                _BlockSpec(
                    id=block_id,
                    class_name=class_name,
                    name=name,
                    title=_block_title(block),
                    medium=str(medium) if medium is not None else None,
                    station_inputs=station_inputs,
                    station_outputs=station_outputs,
                    signal_inputs=signal_inputs,
                    signal_outputs=signal_outputs,
                )
            )
        return specs

    def _make_nodes(self, specs: list[_BlockSpec]) -> list[dict[str, Any]]:
        """Make a top-to-bottom starting layout in one column per fluid medium."""
        media: list[str] = []
        for spec in specs:
            if spec.medium and spec.medium not in media:
                media.append(spec.medium)

        lane_x = {medium: 110 + 275 * index for index, medium in enumerate(media)}
        signal_lane_x = 110 + 275 * max(len(media), 1)
        lane_counts: dict[str, int] = {medium: 0 for medium in media}
        signal_count = 0
        nodes: list[dict[str, Any]] = []

        for spec in specs:
            max_ports = max(
                len(spec.station_inputs),
                len(spec.station_outputs),
                len(spec.signal_inputs),
                len(spec.signal_outputs),
                1,
            )
            height = max(70, 46 + 16 * max_ports)
            width = 188

            if spec.medium:
                x = lane_x[spec.medium]
                y = 100 + 165 * lane_counts[spec.medium]
                lane_counts[spec.medium] += 1
            else:
                x = signal_lane_x
                y = 100 + 165 * signal_count
                signal_count += 1

            ports = []
            # The browser later reassigns these to whichever side faces the
            # connected block. These are only sensible initial orientations.
            ports.extend(
                {"id": f"station-in:{key}", "group": "station-top"}
                for key in spec.station_inputs
            )
            ports.extend(
                {"id": f"station-out:{key}", "group": "station-bottom"}
                for key in spec.station_outputs
            )
            ports.extend(
                {"id": f"signal-in:{key}", "group": "signal-left"}
                for key in spec.signal_inputs
            )
            ports.extend(
                {"id": f"signal-out:{key}", "group": "signal-right"}
                for key in spec.signal_outputs
            )

            fill_color, stroke_color = _block_colors(spec)

            nodes.append(
                {
                    "id": spec.id,
                    "kind": "block",
                    "label": f"{spec.name}",
                    "x": x,
                    "y": y,
                    "width": width,
                    "height": height,
                    "ports": ports,
                    "data": {
                        "name": spec.name,
                        "type": spec.title,
                        "medium": spec.medium,
                        "fill": fill_color,
                        "stroke": stroke_color,
                    },
                }
            )

        return nodes

    def _make_edges(
        self,
        specs: list[_BlockSpec],
        nodes: list[dict[str, Any]],
    ) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
        """Connect producer and consumer ports, adding sources for boundary data."""
        node_by_id = {node["id"]: node for node in nodes}
        station_producers: dict[str, list[_BlockSpec]] = {}
        signal_producers: dict[str, list[_BlockSpec]] = {}
        for spec in specs:
            for key in spec.station_outputs:
                station_producers.setdefault(key, []).append(spec)
            for key in spec.signal_outputs:
                signal_producers.setdefault(key, []).append(spec)

        edges: list[dict[str, Any]] = []
        external_nodes: list[dict[str, Any]] = []
        external_node_ids: set[str] = set()

        def add_external_node(key: str, kind: str, consumer: _BlockSpec) -> str:
            node_id = _stable_id("source", kind, key)
            if node_id in external_node_ids:
                return node_id
            external_node_ids.add(node_id)

            consumer_node = node_by_id[consumer.id]
            if kind == "station":
                # Boundary fluid states normally enter the cycle from above.
                x = int(consumer_node["x"]) + 6
                y = max(20, int(consumer_node["y"]) - 105)
                port_group = "station-bottom"
                port_id = f"station-out:{key}"
                label = f"Boundary station\n{key}"
            else:
                # Scalar boundary conditions are conventionally kept beside
                # their consumer rather than in the fluid-flow columns.
                x = max(20, int(consumer_node["x"]) - 195)
                y = int(consumer_node["y"]) + 10
                port_group = "signal-right"
                port_id = f"signal-out:{key}"
                label = f"Boundary signal\n{key}"

            external_nodes.append(
                {
                    "id": node_id,
                    "kind": "source",
                    "label": label,
                    "x": x,
                    "y": y,
                    "width": 175,
                    "height": 48,
                    "ports": [{"id": port_id, "group": port_group}],
                    "data": {"key": key, "kind": kind},
                }
            )
            return node_id

        def add_edges(kind: str, producers: Mapping[str, list[_BlockSpec]]) -> None:
            for consumer in specs:
                keys = consumer.station_inputs if kind == "station" else consumer.signal_inputs
                for key in keys:
                    source_specs = producers.get(key, [])
                    if source_specs:
                        for producer in source_specs:
                            edge_id = _stable_id("edge", kind, key, producer.id, consumer.id)
                            edges.append(
                                {
                                    "id": edge_id,
                                    "kind": kind,
                                    "label": key,
                                    "source": {
                                        "cell": producer.id,
                                        "port": f"{kind}-out:{key}",
                                    },
                                    "target": {
                                        "cell": consumer.id,
                                        "port": f"{kind}-in:{key}",
                                    },
                                }
                            )
                    else:
                        source_id = add_external_node(key, kind, consumer)
                        edge_id = _stable_id("edge", kind, key, source_id, consumer.id)
                        edges.append(
                            {
                                "id": edge_id,
                                "kind": kind,
                                "label": key,
                                "source": {
                                    "cell": source_id,
                                    "port": f"{kind}-out:{key}",
                                },
                                "target": {
                                    "cell": consumer.id,
                                    "port": f"{kind}-in:{key}",
                                },
                            }
                        )

        add_edges("station", station_producers)
        add_edges("signal", signal_producers)
        return edges, external_nodes


def make_network_viz(
    network: Any,
    *,
    title: str = "PyskyFire Engine Network",
    layout: Mapping[str, Any] | str | Path | None = None,
) -> NetworkVisualizer:
    """Create an editable HTML visualisation of an :class:`EngineNetwork`.

    ``layout`` can be the JSON file downloaded from a prior visualisation.
    """
    return NetworkVisualizer(network=network, title=title, layout=layout)


_HTML_TEMPLATE = r'''<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width,initial-scale=1">
  <title>__TITLE__</title>
  <script>
  __X6_SOURCE__
  </script>
  <style>
    :root { color-scheme: light; }
    * { box-sizing: border-box; }
    body { margin: 0; font-family: system-ui, -apple-system, "Segoe UI", sans-serif; color: #17212b; background: #ffffff; }
    .toolbar { display: flex; align-items: center; gap: 8px; flex-wrap: wrap; padding: 10px 14px; border-bottom: 1px solid #d9e0e7; background: #f8fafc; }
    .toolbar h1 { margin: 0 auto 0 0; font-size: 16px; font-weight: 650; }
    .toolbar button { padding: 6px 10px; border: 1px solid #aab8c5; border-radius: 5px; background: #ffffff; color: #17212b; cursor: pointer; font: inherit; font-size: 13px; }
    .toolbar button:hover { background: #edf3f8; }
    .hint { width: 100%; color: #5d6975; font-size: 12px; }
    #network-canvas { width: 100%; height: calc(100vh - 87px); min-height: 620px; background: #ffffff; }
    .error { margin: 24px; padding: 14px; border: 1px solid #d14343; border-radius: 6px; color: #7b1f1f; background: #fff5f5; }
  </style>
</head>
<body>
  <div class="toolbar">
    <h1>__TITLE__</h1>
    <div class="hint">
      Drag blocks to arrange the schematic. Scroll to pan; Ctrl/Cmd + scroll to zoom.
    </div>
  </div>
  <div id="network-canvas"></div>

  <script id="network-data" type="application/json">__PAYLOAD__</script>
  <script>
  (() => {
    const payload = JSON.parse(document.getElementById('network-data').textContent);
    const container = document.getElementById('network-canvas');

    if (!window.X6) {
      container.innerHTML = '<div class="error">Could not load AntV X6. Check your network connection, or pass a local X6 bundle through <code>x6_url</code>.</div>';
      return;
    }

    function makePortGroup(position, stroke) {
      return {
        position,
        attrs: {
          circle: {
            r: 4,
            magnet: false,
            stroke,
            strokeWidth: 1.25,
            fill: '#ffffff',
          },
        },
      };
    }

    const portGroups = {
      'station-left': makePortGroup('left', '#355267'),
      'station-right': makePortGroup('right', '#355267'),
      'station-top': makePortGroup('top', '#355267'),
      'station-bottom': makePortGroup('bottom', '#355267'),
      'signal-left': makePortGroup('left', '#7a5a18'),
      'signal-right': makePortGroup('right', '#7a5a18'),
      'signal-top': makePortGroup('top', '#7a5a18'),
      'signal-bottom': makePortGroup('bottom', '#7a5a18'),
    };

    const graph = new X6.Graph({
      container,
      width: Math.max(container.clientWidth, 900),
      height: Math.max(container.clientHeight, 620),
      grid: { visible: true, size: 10, type: 'mesh', args: [{ color: '#e8edf2', thickness: 1 }] },
      panning: true,
      mousewheel: { enabled: true, modifiers: ['ctrl', 'meta'] },
      connecting: { allowBlank: false, allowNode: false, allowEdge: false, allowPort: false },
    });

    const defaultLayout = structuredClone(payload.layout || { nodes: {}, edges: {} });
    let activeLayout = structuredClone(defaultLayout);

    function nodePosition(node) {
      const saved = activeLayout.nodes && activeLayout.nodes[node.id];
      if (saved && Number.isFinite(saved.x) && Number.isFinite(saved.y)) return saved;
      return { x: node.x, y: node.y };
    }

    for (const node of payload.nodes) {
      const position = nodePosition(node);
      const isSource = node.kind === 'source';
      const fillColor = isSource
        ? '#f7f9fb'
        : ((node.data && node.data.fill) || '#eef5fb');
      const strokeColor = isSource
        ? '#8092a1'
        : ((node.data && node.data.stroke) || '#2f6f96');

      graph.addNode({
        id: node.id,
        shape: 'rect',
        x: position.x,
        y: position.y,
        width: node.width,
        height: node.height,
        label: node.label,
        attrs: {
          body: {
            fill: fillColor,
            stroke: strokeColor,
            strokeWidth: 1.35,
            rx: 7,
            ry: 7,
          },
          label: {
            fill: '#17212b',
            fontSize: 12,
            fontWeight: 560,
            textWrap: { width: -16, height: -12, ellipsis: true },
          },
        },
        ports: { groups: portGroups, items: node.ports },
        data: node.data,
      });
    }

    for (const edge of payload.edges) {
      const isSignal = edge.kind === 'signal';
      const config = {
        id: edge.id,
        source: edge.source,
        target: edge.target,
        router: { name: 'manhattan', args: { padding: 16 } },
        connector: { name: 'rounded', args: { radius: 6 } },
        attrs: {
          line: {
            stroke: isSignal ? '#9a751e' : '#2f6f96',
            strokeWidth: isSignal ? 1.6 : 2.1,
            strokeDasharray: isSignal ? '6 4' : '',
            targetMarker: { name: 'classic', size: 7 },
          },
        },
        labels: [{ position: 0.5, attrs: { label: { text: edge.label, fill: isSignal ? '#79590e' : '#1f4e6b', fontSize: 10 } } }],
        data: { kind: edge.kind },
      };
      graph.addEdge(config);
    }

    function nodeCenter(node) {
      const bbox = node.getBBox();
      return {
        x: bbox.x + bbox.width / 2,
        y: bbox.y + bbox.height / 2,
      };
    }

    function sideToward(dx, dy) {
      return Math.abs(dx) > Math.abs(dy)
        ? (dx >= 0 ? 'right' : 'left')
        : (dy >= 0 ? 'bottom' : 'top');
    }

    function updatePortSides() {
      // A port faces the average position of the blocks connected to it.
      // This is a local routing heuristic: it makes the first and final
      // Manhattan segments as short and direct as possible without imposing
      // a global layout algorithm.
      const connections = new Map();

      function accumulate(nodeId, portId, kind, otherNodeId) {
        const node = graph.getCellById(nodeId);
        const other = graph.getCellById(otherNodeId);
        if (!node || !other || !portId) return;

        const own = nodeCenter(node);
        const peer = nodeCenter(other);
        const mapKey = `${nodeId}\u0000${portId}`;
        const item = connections.get(mapKey) || {
          node,
          portId,
          kind,
          dx: 0,
          dy: 0,
        };
        item.dx += peer.x - own.x;
        item.dy += peer.y - own.y;
        connections.set(mapKey, item);
      }

      for (const edge of payload.edges) {
        accumulate(edge.source.cell, edge.source.port, edge.kind, edge.target.cell);
        accumulate(edge.target.cell, edge.target.port, edge.kind, edge.source.cell);
      }

      for (const item of connections.values()) {
        const group = `${item.kind}-${sideToward(item.dx, item.dy)}`;
        const port = item.node.getPort(item.portId);
        if (port && port.group !== group) {
          item.node.setPortProp(item.portId, 'group', group);
        }
      }
    }

    let portUpdateFrame = null;
    function schedulePortSideUpdate() {
      if (portUpdateFrame !== null) return;
      portUpdateFrame = window.requestAnimationFrame(() => {
        portUpdateFrame = null;
        updatePortSides();
      });
    }

    graph.on('node:change:position', schedulePortSideUpdate);
    updatePortSides();

    const nodeUndoStack = [];
    const maxUndoEntries = 100;
    let dragStart = null;

    function isEditableElement(element) {
    return (
        element instanceof HTMLInputElement ||
        element instanceof HTMLTextAreaElement ||
        element instanceof HTMLSelectElement ||
        element?.isContentEditable
    );
    }

    graph.on('node:move', ({ node }) => {
    const position = node.position();
    dragStart = {
        id: node.id,
        x: position.x,
        y: position.y,
    };
    });

    graph.on('node:moved', ({ node }) => {
    if (!dragStart || dragStart.id !== node.id) return;

    const position = node.position();

    if (position.x !== dragStart.x || position.y !== dragStart.y) {
        nodeUndoStack.push(dragStart);

        if (nodeUndoStack.length > maxUndoEntries) {
        nodeUndoStack.shift();
        }
    }

    dragStart = null;
    });

    document.addEventListener('keydown', (event) => {
    const isUndoShortcut =
        (event.ctrlKey || event.metaKey) &&
        !event.shiftKey &&
        event.key.toLowerCase() === 'z';

    if (!isUndoShortcut || isEditableElement(event.target)) return;

    const previous = nodeUndoStack.pop();
    if (!previous) return;

    const node = graph.getCellById(previous.id);
    if (!node || !node.isNode()) return;

    event.preventDefault();
    node.position(previous.x, previous.y);
    });

    const resize = () => graph.resize(Math.max(container.clientWidth, 900), Math.max(container.clientHeight, 620));
    new ResizeObserver(resize).observe(container);
    resize();
  })();
  </script>
</body>
</html>
'''
