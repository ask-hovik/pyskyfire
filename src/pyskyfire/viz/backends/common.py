"""Scalar display validation shared by desktop and HTML backends."""
import numpy as np


def is_visible(node):
    return node.visible and (node._parent is None or is_visible(node._parent))


def _style_for(scene, name, clim):
    fields = [n.fields[name] for n in scene.walk() if name in n.fields and is_visible(n)]
    if not fields:
        raise ValueError(f"No visible component provides field {name!r}")
    if len({f.units for f in fields}) > 1:
        raise ValueError("A shared color field must have consistent SI units")
    finite = [f.values[np.isfinite(f.values)] for f in fields]
    finite = np.concatenate(finite)
    if clim is not None:
        limits = np.asarray(clim, dtype=float)
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
    return {"name": name, "label": fields[0].label,
            "units": fields[0].units, "limits": limits.tolist(),
            "has_values": bool(finite.size)}


def field_style(scene):
    if scene.color_by is None:
        return None
    return _style_for(scene, scene.color_by, scene.clim)


def field_styles(scene):
    """Color ranges for every field the HTML viewer's field selector can offer.

    ``clim`` pins the range of the initially selected field only; the rest are
    ranged over their own visible values. A field whose units disagree between
    components is omitted rather than offered with a meaningless shared scale,
    except when it is ``color_by``, where it still raises.
    """
    names = list(dict.fromkeys(name for node in scene.walk() if is_visible(node)
                               for name in node.fields))
    styles = {}
    for name in names:
        clim = scene.clim if name == scene.color_by else None
        try:
            styles[name] = _style_for(scene, name, clim)
        except ValueError:
            if name == scene.color_by:
                raise
    return styles
