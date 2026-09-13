"""Keyword-dispatched design inputs for cooling-channel geometry.

Most user-facing geometry inputs in :mod:`pyskyfire.regen` accept either a
constant or a function of position.  Historically every such function had to
take exactly one positional argument ``x``, even when the value did not depend
on it.  This module replaces that convention: a design-input function declares
the quantities it needs *by parameter name*, and only those are supplied.

    rib_width = lambda: 1.0e-3                       # constant
    rib_width = lambda x: 1.0e-3 + 1.0e-5 * x        # grows downstream
    rib_width = lambda pitch_mid: 0.1 * pitch_mid    # 10 % of the local pitch

Which names are legal depends on where the function is plugged in.  The
geometry build order forms a dependency chain, so a function may only ask for
quantities that are already known when its slot is evaluated.  Those
quantities are grouped into three tiers:

``A``
    Contour geometry.  Known before any channel exists.
``B``
    Channel layout and wall stack.  Known once the circuits are configured.
``C``
    The resolved channel section.  Known only after heights, centerlines and
    angular widths have been built.

See ``docs/explanations/design-input-functions.md`` for the user-facing
description of the tiers and an exhaustive glossary.
"""

from __future__ import annotations

import inspect
from typing import Any, Callable, Iterable, Mapping

import numpy as np

__all__ = [
    "DesignInputError",
    "DesignContext",
    "TIER_A",
    "TIER_B",
    "TIER_C",
    "VOCABULARY",
    "SLOT_KEYS",
    "make_evaluator",
    "evaluate",
    "evaluate_profile",
    "contour_providers",
    "layout_providers",
    "section_providers",
]


class DesignInputError(ValueError):
    """Raised when a design-input function cannot be dispatched or evaluated."""


# ---------------------------------------------------------------------------
# Vocabulary
# ---------------------------------------------------------------------------
# Every name a design-input function may declare, with the one-line description
# used both in error messages and in the generated documentation glossary.

TIER_A: dict[str, str] = {
    "x": "Axial coordinate measured along the engine axis [m].",
    "r": "Hot-wall (contour) radius at ``x`` [m].",
    "circumference": "Hot-wall circumference ``2*pi*r`` at ``x`` [m].",
    "dr_dx": "Contour slope ``dr/dx`` at ``x`` [-].",
    "normal_angle": "Angle between the outward wall normal and the plane "
    "normal to the engine axis [rad].",
    "x_frac": "Signed normalised axial coordinate: -1 at the chamber start, "
    "0 at the throat, +1 at the nozzle exit [-].",
    "area_ratio": "Local area ratio ``(r/r_throat)**2`` [-].",
    "x_throat": "Axial coordinate of the throat [m].",
    "r_throat": "Throat radius [m].",
    "x_start": "Axial coordinate of the contour start [m].",
    "x_end": "Axial coordinate of the contour end [m].",
    "r_chamber": "Contour radius at the chamber start [m].",
    "r_exit": "Contour radius at the nozzle exit [m].",
    "contour": "The :class:`~pyskyfire.regen.contour.Contour` object itself, "
    "as an escape hatch for quantities not listed here.",
}

TIER_B: dict[str, str] = {
    "n_channel_positions": "Number of circumferential channel positions "
    "(leaves) in this circuit [-].",
    "n_channels": "Total number of channels in this circuit, i.e. leaves "
    "times channels per leaf [-].",
    "n_channels_all_circuits": "Number of surface-occluding channels present "
    "at ``x``, summed over every circuit [-].",
    "t_wall": "Combined thickness of all main wall layers forming the "
    "channel walls, excluding liners [m].",
    "t_wall_total": "Combined thickness of the whole wall stack [m].",
    "t_hot_layer": "Combined thickness of layers explicitly marked as liners "
    "on the combustion side of the main walls [m].",
}

TIER_C: dict[str, str] = {
    "h": "Supplied channel height [m], using the cross section's height_input "
    "convention (full height by default, optionally rounded base height).",
    "gamma": "Local helix angle of the channel, measured from the meridional "
    "direction [rad].",
    "theta": "Angular pitch of one channel measured in the plane normal to "
    "the channel axis [rad]. Already reduced by ``cos(gamma)``.",
    "r_channel_base": "Radius of the cross-section base point: the channel "
    "floor for a squared section, the inner lobe centre for a rounded one [m].",
    "r_channel_top": "Radius of the channel roof, ``r_channel_base + h`` [m].",
    "pitch": "Channel-to-channel pitch at the hot wall, ``r * theta`` [m].",
    "pitch_base": "Channel-to-channel pitch at the channel floor, "
    "``r_channel_base * theta`` [m].",
    "pitch_mid": "Channel-to-channel pitch at mid channel height, "
    "``(r_channel_base + h/2) * theta`` [m].",
}

VOCABULARY: dict[str, str] = {**TIER_A, **TIER_B, **TIER_C}

_TIER_OF: dict[str, str] = {
    **{name: "A" for name in TIER_A},
    **{name: "B" for name in TIER_B},
    **{name: "C" for name in TIER_C},
}

# Wall thicknesses may not ask about the wall stack they define, nor about
# anything downstream of it.
_WALL_KEYS = frozenset({"t_wall", "t_wall_total", "t_hot_layer"})

_A = frozenset(TIER_A)
_B = frozenset(TIER_B)
_C = frozenset(TIER_C)

#: Which vocabulary names each design-input slot accepts.
SLOT_KEYS: dict[str, frozenset[str]] = {
    # Tier A only: the wall stack is not yet defined while it is being defined.
    "wall_thickness": _A,
    # Tier A + B: everything that is a pure function of the contour and the
    # circuit configuration.
    "channel_height": _A | _B,
    # Channel heights are built before any helix angle is evaluated, and a
    # height may not depend on the helix, so 'h' is available here without
    # creating a cycle. It is what an area-preserving helix law needs.
    "helix_angle": _A | _B | frozenset({"h"}),
    # Tier A + B, plus the two section quantities that are built before the
    # angular widths are.
    "channel_angular_width": _A | _B | frozenset({"h", "gamma"}),
    # Evaluated against a fully resolved section: everything is available.
    "rib_width": _A | _B | _C,
    "roughness": _A | _B | _C,
    "hot_gas_surface_area_multiplier": _A | _B | _C,
}


def _slot_keys(slot: str) -> frozenset[str]:
    try:
        return SLOT_KEYS[slot]
    except KeyError:  # pragma: no cover - guards internal typos
        raise DesignInputError(
            f"unknown design-input slot {slot!r}; "
            f"expected one of {sorted(SLOT_KEYS)}"
        ) from None


# ---------------------------------------------------------------------------
# Context
# ---------------------------------------------------------------------------


class DesignContext(Mapping):
    """Lazily evaluated bundle of design-input quantities at one station.

    Providers are either plain values or zero-argument callables.  A callable
    is invoked at most once, the first time its name is requested, so building
    a context for a station that only needs ``x`` costs nothing.
    """

    __slots__ = ("_providers", "_cache")

    def __init__(self, providers: Mapping[str, Any]):
        self._providers = providers
        self._cache: dict[str, Any] = {}

    def __getitem__(self, key: str) -> Any:
        try:
            return self._cache[key]
        except KeyError:
            pass
        try:
            provider = self._providers[key]
        except KeyError:
            raise DesignInputError(
                f"{key!r} is not available in this context"
            ) from None
        value = provider() if callable(provider) else provider
        self._cache[key] = value
        return value

    def __iter__(self):
        return iter(self._providers)

    def __len__(self) -> int:
        return len(self._providers)


# ---------------------------------------------------------------------------
# Provider factories
# ---------------------------------------------------------------------------


def contour_providers(contour, x: float) -> dict[str, Any]:
    """Return the Tier A providers for one station on ``contour``."""

    def _need_contour():
        if contour is None:
            raise DesignInputError(
                "this design-input function requests a contour-derived "
                "quantity, but no contour is attached to the object that "
                "owns it"
            )
        return contour

    def _x_frac():
        c = _need_contour()
        x_t = float(c.x_t)
        if x < x_t:
            span = x_t - float(c.xs[0])
            return (x - x_t) / span if span > 0.0 else 0.0
        span = float(c.xs[-1]) - x_t
        return (x - x_t) / span if span > 0.0 else 0.0

    def _area_ratio():
        c = _need_contour()
        r_t = float(c.r_t)
        return (float(c.r(x)) / r_t) ** 2 if r_t > 0.0 else float("inf")

    return {
        "x": float(x),
        "r": lambda: float(_need_contour().r(x)),
        "circumference": lambda: 2.0 * np.pi * float(_need_contour().r(x)),
        "dr_dx": lambda: float(_need_contour().dr_dx(x)),
        "normal_angle": lambda: float(_need_contour().normal_angle(x)),
        "x_frac": _x_frac,
        "area_ratio": _area_ratio,
        "x_throat": lambda: float(_need_contour().x_t),
        "r_throat": lambda: float(_need_contour().r_t),
        "x_start": lambda: float(_need_contour().xs[0]),
        "x_end": lambda: float(_need_contour().xs[-1]),
        "r_chamber": lambda: float(_need_contour().r_c),
        "r_exit": lambda: float(_need_contour().r_e),
        "contour": _need_contour,
    }


def layout_providers(
    *,
    n_channel_positions=None,
    n_channels=None,
    n_channels_all_circuits=None,
    t_wall=None,
    t_wall_total=None,
    t_hot_layer=None,
) -> dict[str, Any]:
    """Return the Tier B providers, skipping any that are unavailable."""

    candidates = {
        "n_channel_positions": n_channel_positions,
        "n_channels": n_channels,
        "n_channels_all_circuits": n_channels_all_circuits,
        "t_wall": t_wall,
        "t_wall_total": t_wall_total,
        "t_hot_layer": t_hot_layer,
    }
    return {key: value for key, value in candidates.items() if value is not None}


def section_providers(
    *,
    r=None,
    r_base=None,
    h=None,
    gamma=None,
    theta=None,
    t_wall=None,
) -> dict[str, Any]:
    """Return the Tier C providers derived from a resolved channel section.

    ``r`` is the hot-wall (contour) radius and ``r_base`` the resolved
    cross-section base point. They differ by the wall stack projected onto the
    contour normal, so ``r_base`` is not simply ``r + t_wall``; it defaults to
    ``r`` only for a section that has not been placed yet.

    Quantities whose inputs are still unknown are simply omitted, which is how
    a slot evaluated part-way through the geometry build ends up seeing only
    the part of tier C that already exists.
    """

    providers: dict[str, Any] = {}
    if gamma is not None:
        providers["gamma"] = float(gamma)
    if h is not None:
        providers["h"] = float(h)
    if h is None or r is None or t_wall is None:
        return providers

    r = float(r)
    h = float(h)
    r_base = r if r_base is None else float(r_base)
    providers["r_channel_base"] = r_base
    providers["r_channel_top"] = r_base + h
    if theta is None:
        return providers

    theta = float(theta)
    providers["theta"] = theta
    providers["pitch"] = r * theta
    providers["pitch_base"] = r_base * theta
    providers["pitch_mid"] = (r_base + 0.5 * h) * theta
    return providers


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------


def _describe_rejection(key: str, slot: str, allowed: Iterable[str]) -> str:
    allowed = sorted(allowed)
    if key in VOCABULARY:
        tier = _TIER_OF[key]
        if key in _WALL_KEYS and slot == "wall_thickness":
            detail = (
                f"{key!r} is a tier {tier} quantity describing the wall stack, "
                "which is what this function defines; requesting it would be "
                "circular"
            )
        else:
            detail = (
                f"{key!r} is a tier {tier} quantity and is not yet known when "
                f"{slot!r} is evaluated"
            )
    else:
        detail = f"{key!r} is not a design-input quantity"
    return f"{detail}. Available here: {', '.join(allowed)}"


def _signature_keys(function: Callable, slot: str, label: str):
    """Return the ordered parameter names ``function`` wants, or ``None``.

    ``None`` means "pass every available name", which is what a function with
    a ``**kwargs`` catch-all asks for.
    """

    allowed = _slot_keys(slot)
    try:
        signature = inspect.signature(function)
    except (TypeError, ValueError):  # builtins without introspectable signatures
        raise DesignInputError(
            f"{label}: cannot inspect the signature of {function!r}; "
            "wrap it in a plain Python function that names the quantities it "
            f"needs (available here: {', '.join(sorted(allowed))})"
        ) from None

    keys: list[str] = []
    for name, parameter in signature.parameters.items():
        kind = parameter.kind
        if kind is inspect.Parameter.VAR_KEYWORD:
            return None
        if kind is inspect.Parameter.VAR_POSITIONAL:
            raise DesignInputError(
                f"{label}: *{name} is not supported. Design-input functions "
                "are called with keyword arguments; declare the quantities "
                f"you need by name (available here: "
                f"{', '.join(sorted(allowed))})"
            )
        if kind is inspect.Parameter.POSITIONAL_ONLY:
            raise DesignInputError(
                f"{label}: positional-only parameter {name!r} is not "
                "supported. Design-input functions are called with keyword "
                "arguments."
            )
        if name not in allowed:
            raise DesignInputError(
                f"{label}: {_describe_rejection(name, slot, allowed)}"
            )
        keys.append(name)
    return tuple(keys)


def make_evaluator(value, *, slot: str, label: str) -> Callable[[DesignContext], float]:
    """Return ``f(context) -> float`` for a constant or design-input function.

    The callable's signature is inspected once, here, rather than at every
    station.

    The returned evaluator carries two attributes that let hot callers skip
    context construction: ``required_keys`` is the tuple of names it will read
    (``None`` when the function has a ``**kwargs`` catch-all), and
    ``constant_value`` is the float it always returns, or ``None``.
    """

    if not callable(value):
        try:
            constant = float(value)
        except (TypeError, ValueError):
            raise DesignInputError(
                f"{label}: expected a float or a callable, got "
                f"{type(value).__name__}"
            ) from None

        def evaluate_constant(context: DesignContext) -> float:
            return constant

        return _tag(evaluate_constant, (), constant)

    keys = _signature_keys(value, slot, label)
    allowed = _slot_keys(slot)

    if keys is None:
        def evaluate_all(context: DesignContext) -> float:
            arguments = {
                key: context[key] for key in allowed if key in context
            }
            return _as_float(value(**arguments), label)

        return _tag(evaluate_all, None, None)

    if not keys:
        def evaluate_nullary(context: DesignContext) -> float:
            return _as_float(value(), label)

        return _tag(evaluate_nullary, (), None)

    def evaluate_named(context: DesignContext) -> float:
        arguments = {key: context[key] for key in keys}
        return _as_float(value(**arguments), label)

    return _tag(evaluate_named, keys, None)


def _tag(evaluator, keys, constant):
    evaluator.required_keys = keys
    evaluator.constant_value = constant
    return evaluator


def _as_float(result, label: str) -> float:
    array = np.asarray(result, dtype=float)
    if array.ndim != 0:
        raise DesignInputError(
            f"{label}: design-input functions are evaluated one station at a "
            f"time and must return a scalar, got an array of shape "
            f"{array.shape}"
        )
    return float(array)


def evaluate(value, context: DesignContext, *, slot: str, label: str) -> float:
    """Evaluate a single design input at a single station."""

    return make_evaluator(value, slot=slot, label=label)(context)


def evaluate_profile(
    value,
    contexts: Iterable[DesignContext],
    *,
    slot: str,
    label: str,
) -> np.ndarray:
    """Evaluate a design input at every station in ``contexts``."""

    evaluator = make_evaluator(value, slot=slot, label=label)
    return np.array([evaluator(context) for context in contexts], dtype=float)
