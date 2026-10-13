"""Sense of rotation about a component's own axis.

An impeller spins one way about its shaft and the volute around it must wind
the same way, since a volute wound against the discharge swirl collects the
flow backwards. Both state the sense with the same word here, so that a design
can be checked by reading it rather than by looking at a render.

The sense is named as seen looking down the component's own axis from its
positive end, with the axis pointing at the viewer. Counterclockwise is then
the right-handed sense about that axis and carries sign ``+1``. Which axis
that is belongs to the component: ``+z`` for a volute, ``+x`` for an impeller.

Naming the sense about the part's own axis keeps it independent of how the
part is later placed. A rigid placement carries the sense along with the axis,
so two parts agree exactly when their placed axes point the same way and their
senses match -- and disagree when a placement reverses one of the axes, which
is the failure this vocabulary exists to make visible.
"""

from __future__ import annotations

__all__ = [
    "CLOCKWISE",
    "COUNTERCLOCKWISE",
    "parse_rotation",
    "rotation_sign",
]

CLOCKWISE = "clockwise"
COUNTERCLOCKWISE = "counterclockwise"

_ALIASES = {
    "cw": CLOCKWISE,
    "clockwise": CLOCKWISE,
    "ccw": COUNTERCLOCKWISE,
    "counterclockwise": COUNTERCLOCKWISE,
}


def parse_rotation(value) -> str:
    """Normalise a rotation sense to :data:`CLOCKWISE` or :data:`COUNTERCLOCKWISE`.

    Accepts ``'cw'``, ``'clockwise'``, ``'ccw'`` and ``'counterclockwise'``,
    in any capitalisation.
    """
    if not isinstance(value, str):
        raise TypeError(
            f"rotation must be one of {sorted(_ALIASES)}, not a "
            f"{type(value).__name__}"
        )
    try:
        return _ALIASES[value.strip().lower()]
    except KeyError:
        raise ValueError(
            f"rotation must be one of {sorted(_ALIASES)}; got {value!r}"
        ) from None


def rotation_sign(value) -> float:
    """``+1.0`` for a counterclockwise sense, ``-1.0`` for a clockwise one."""
    return 1.0 if parse_rotation(value) == COUNTERCLOCKWISE else -1.0
