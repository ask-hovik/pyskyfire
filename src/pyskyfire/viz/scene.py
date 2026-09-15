"""Renderer-independent scenes. Coordinates are metres; rotations are radians."""
from __future__ import annotations

import base64
from dataclasses import dataclass, field
from pathlib import Path
from typing import Protocol

import numpy as np
from scipy.spatial.transform import Rotation

__all__ = ["Scene", "SceneNode", "MeshData", "FieldData", "Transform", "Anchor", "Component"]


class Transform:
    """Rigid local-to-parent transform, using column vectors.

    ``rotation`` is a scipy Rotation or extrinsic xyz Euler angles in radians.
    A supplied ``matrix`` is mutually exclusive with translation and rotation.
    """

    def __init__(self, translation=None, rotation=None, *, matrix=None):
        if matrix is not None:
            if translation is not None or rotation is not None:
                raise ValueError("matrix cannot be combined with translation or rotation")
            m = np.array(matrix, dtype=float, copy=True)
        else:
            m = np.eye(4)
            if translation is not None:
                t = np.asarray(translation, dtype=float)
                if t.shape != (3,):
                    raise ValueError("translation must have shape (3,)")
                m[:3, 3] = t
            if rotation is not None:
                r = rotation if isinstance(rotation, Rotation) else Rotation.from_euler("xyz", rotation)
                m[:3, :3] = r.as_matrix()
        if (m.shape != (4, 4) or not np.isfinite(m).all()
                or not np.allclose(m[3], [0, 0, 0, 1])
                or not np.allclose(m[:3, :3].T @ m[:3, :3], np.eye(3), atol=1e-9)
                or not np.isclose(np.linalg.det(m[:3, :3]), 1, atol=1e-9)):
            raise ValueError("Expected a finite rigid 4x4 transform (no scaling or reflection)")
        m.setflags(write=False)
        self.matrix = m

    def __matmul__(self, other):
        return Transform(matrix=self.matrix @ other.matrix)

    def inverse(self):
        return Transform(matrix=np.linalg.inv(self.matrix))

    def apply(self, points):
        return np.asarray(points) @ self.matrix[:3, :3].T + self.matrix[:3, 3]


@dataclass(frozen=True)
class Anchor:
    """Named coordinate frame. The +x axis is the outward port direction."""

    transform: Transform = field(default_factory=Transform)

    @classmethod
    def from_direction(cls, position, direction):
        x = np.asarray(direction, dtype=float)
        if x.shape != (3,) or not np.isfinite(x).all() or np.linalg.norm(x) < 1e-12:
            raise ValueError("Anchor direction must be a finite nonzero vector")
        x = x / np.linalg.norm(x)
        up = np.eye(3)[np.argmin(np.abs(x))]
        y = np.cross(up, x)
        y /= np.linalg.norm(y)
        m = np.eye(4)
        m[:3, :3] = np.column_stack((x, y, np.cross(x, y)))
        m[:3, 3] = position
        return cls(Transform(matrix=m))

    @property
    def position(self):
        return self.transform.matrix[:3, 3].copy()


@dataclass(frozen=True)
class MeshData:
    """Surface geometry in local Cartesian coordinates; indices describe triangles."""

    points: np.ndarray
    triangles: np.ndarray

    def __post_init__(self):
        p = np.array(self.points, dtype=float, copy=True)
        raw = np.asarray(self.triangles)
        if raw.dtype.kind not in "iu":
            raise ValueError("Triangle indices must be integers")
        t = np.array(raw, dtype=np.int64, copy=True)
        if p.ndim != 2 or p.shape[1] != 3 or not np.isfinite(p).all():
            raise ValueError("Mesh points must be finite with shape (N,3)")
        if t.ndim != 2 or t.shape[1] != 3 or (t.size and (t.min() < 0 or t.max() >= len(p))):
            raise ValueError("Triangles must have shape (M,3) and valid vertex indices")
        p.setflags(write=False)
        t.setflags(write=False)
        object.__setattr__(self, "points", p)
        object.__setattr__(self, "triangles", t)


@dataclass(frozen=True)
class FieldData:
    """Scalar mesh values in base SI units. NaN represents unavailable data.

    Units are descriptive metadata, never conversion inputs. Association is
    ``point`` or ``cell`` (one value per triangle).
    """

    values: np.ndarray
    units: str
    label: str
    association: str = "point"

    def __post_init__(self):
        values = np.array(self.values, dtype=float, copy=True)
        if values.ndim != 1 or np.isinf(values).any():
            raise ValueError("Field values must be a 1D array without infinities")
        if self.association not in ("point", "cell"):
            raise ValueError("Field association must be 'point' or 'cell'")
        values.setflags(write=False)
        object.__setattr__(self, "values", values)


@dataclass
class SceneNode:
    """A selectable part or group, independent of the renderer's actor layout."""

    name: str
    mesh: MeshData | None = None
    fields: dict[str, FieldData] = field(default_factory=dict)
    children: list[SceneNode] = field(default_factory=list)
    anchors: dict[str, Anchor] = field(default_factory=dict)
    transform: Transform = field(default_factory=Transform)
    color: str = "lightsteelblue"
    opacity: float = 1.0
    visible: bool = True
    metadata: dict = field(default_factory=dict)
    _parent: SceneNode | None = field(default=None, init=False, repr=False)

    @property
    def id(self):
        return f"{self._parent.id}/{self.name}" if self._parent else self.name

    @property
    def world_transform(self):
        return self._parent.world_transform @ self.transform if self._parent else self.transform

    def anchor(self, name):
        """Return a snapshot of a named anchor in scene coordinates."""
        return Anchor(self.world_transform @ self.anchors[name].transform)


class Component(Protocol):
    def build(self) -> SceneNode:
        """Return a fresh hierarchy, optionally sharing immutable mesh data."""
        ...


class Scene:
    """Assemble components for desktop or standalone HTML display.

    ``anchor`` names the added component's local anchor; ``attach_to`` is a
    world-space Anchor. Attachment mates their outward x axes in opposition.
    It is a one-time placement calculation, with no fluid-network constraints.
    Explicit placement and anchor placement are mutually exclusive.

    ``show_axes`` adds a corner orientation triad marking the scene axes
    (x red, y green, z blue). It is an orientation cue only, not a measuring
    grid: it carries no ticks and does not indicate the origin or scale.
    """

    def __init__(self, *, title="Pyskyfire 3D", color_by=None, clim=None, show_axes=False):
        self.title = title
        self.color_by = color_by
        self.clim = clim
        self.show_axes = show_axes
        self.nodes = []

    def add(self, component, *, name=None, transform=None, position=None,
            rotation=None, parent=None, anchor=None, attach_to=None):
        if transform is not None and (position is not None or rotation is not None):
            raise ValueError("transform and position/rotation overconstrain placement")
        attached = anchor is not None or attach_to is not None
        if attached and (anchor is None or attach_to is None):
            raise ValueError("Both anchor and attach_to are required")
        if attached and any(v is not None for v in (transform, position, rotation)):
            raise ValueError("Anchor and explicit placement overconstrain placement")
        if parent is not None and not any(node is parent for node in self.walk()):
            raise ValueError("parent must belong to this scene")
        node = component.build()
        if not isinstance(node, SceneNode):
            raise TypeError("Component.build() must return a SceneNode")
        if any(existing is node for existing in self.walk()):
            raise ValueError("Component.build() must return a fresh hierarchy")
        if name is not None:
            node.name = name
        siblings = self.nodes if parent is None else parent.children
        if any(n.name == node.name for n in siblings):
            raise ValueError(f"Duplicate scene name: {node.name}")
        if attached:
            if not isinstance(attach_to, Anchor):
                raise TypeError("attach_to must be an Anchor")
            mate = Transform(rotation=(0, 0, np.pi))
            world = attach_to.transform @ mate @ node.anchors[anchor].transform.inverse()
            node.transform = parent.world_transform.inverse() @ world if parent else world
        elif transform is not None or position is not None or rotation is not None:
            node.transform = transform if transform is not None else Transform(position, rotation)
        # Validate and adopt before inserting, so failed additions are atomic.
        seen = set()
        existing = {id(n) for n in self.walk()}

        def adopt(n, owner):
            if id(n) in seen or id(n) in existing:
                raise ValueError("Component hierarchy contains reused nodes or a cycle")
            seen.add(id(n))
            if not isinstance(n.name, str) or not n.name or "/" in n.name:
                raise ValueError("Node names must be nonempty and contain no '/'")
            if not isinstance(n.transform, Transform):
                raise TypeError("transform must be a Transform")
            if not np.isfinite(n.opacity) or not 0 <= n.opacity <= 1:
                raise ValueError("opacity must lie between zero and one")
            if len({c.name for c in n.children}) != len(n.children):
                raise ValueError("Sibling names must be unique")
            for f in n.fields.values():
                count = 0 if n.mesh is None else len(n.mesh.points if f.association == "point" else n.mesh.triangles)
                if len(f.values) != count or n.mesh is None:
                    raise ValueError("Field length must match its mesh association")
            n._parent = owner
            for child in n.children:
                adopt(child, n)

        adopt(node, parent)
        siblings.append(node)
        return node

    def walk(self):
        def visit(node):
            yield node
            for child in node.children:
                yield from visit(child)
        for node in self.nodes:
            yield from visit(node)

    def to_plotter(self, *, off_screen=False):
        """Create a PyVista plotter; the caller owns and closes it (also for GUI embedding)."""
        from .backends.pyvista import make_plotter
        return make_plotter(self, off_screen=off_screen)

    def show(self):
        plotter = self.to_plotter()
        try:
            return plotter.show()
        finally:
            plotter.close()
            plotter.deep_clean()

    def to_html(self):
        from .backends.threejs import to_html
        return to_html(self)

    def save_html(self, path):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.to_html(), encoding="utf-8")
        return path

    @property
    def data_url(self):
        """Interactive HTML for Report.add_iframe()."""
        return "data:text/html;base64," + base64.b64encode(self.to_html().encode()).decode("ascii")
