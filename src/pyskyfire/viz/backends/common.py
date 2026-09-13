"""Scalar display validation shared by desktop and HTML backends."""
import numpy as np


def is_visible(node):
    return node.visible and (node._parent is None or is_visible(node._parent))


def field_style(scene):
    if scene.color_by is None:
        return None
    fields = [n.fields[scene.color_by] for n in scene.walk()
              if scene.color_by in n.fields and is_visible(n)]
    if not fields:
        raise ValueError(f"No visible component provides field {scene.color_by!r}")
    if len({f.units for f in fields}) > 1:
        raise ValueError("A shared color field must have consistent SI units")
    finite = [f.values[np.isfinite(f.values)] for f in fields]
    finite = np.concatenate(finite)
    if scene.clim is not None:
        limits = np.asarray(scene.clim, dtype=float)
        if limits.shape != (2,) or not np.isfinite(limits).all() or limits[0] >= limits[1]:
            raise ValueError("clim must contain two finite increasing SI values")
    elif finite.size:
        low, high = float(finite.min()), float(finite.max())
        if low == high:
            pad = max(abs(low)*1e-6, 1e-12)
            low, high = low-pad, high+pad
        limits = np.array([low, high])
    else:
        limits = np.array([0., 1.])
    return {"name": scene.color_by, "label": fields[0].label,
            "units": fields[0].units, "limits": limits.tolist(),
            "has_values": bool(finite.size)}
