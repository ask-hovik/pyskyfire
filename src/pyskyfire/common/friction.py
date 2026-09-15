"""Duct friction physics shared by cooling channels and piping.

This module is deliberately geometry-agnostic. Every entry point takes plain
scalars -- a Reynolds number, a hydraulic diameter, an absolute roughness, a
signed radius of curvature -- so the same physics serves a regen cooling
channel marching on an axial coordinate and a pipe run marching on arc length.
Anything that knows about contours, centerlines or stations belongs elsewhere.

Two interchangeable backends are provided:

``"native"``
    The correlations pyskyfire has always used, inherited from the RL10
    modelling work of Binder et al.  Laminar ``64/Re``, Petukhov (smooth) or
    Colebrook-White (rough) in the turbulent branch, linearly blended across
    ``2300 < Re < 3500``, with an Ito-type curvature multiplier.

``"fluids"``
    The same quantities taken from the ``fluids`` library.  The straight-duct
    factor blends the laminar and Clamond branches smoothly from ``Re=2040``;
    curved ducts use the same transition policy starting at Schmidt's
    curvature-dependent critical Reynolds number. The curvature multiplier is
    the ratio of the blended Schmidt curved factor to the blended straight
    factor, with the same roughness in both.

The native backend remains the default for validation while the ``fluids``
backend is evaluated. Comparing them measures model sensitivity; it does not
establish a bound on modelling uncertainty.

Notes
-----
Both curvature multipliers are helical-coil correlations. They model curvature
as an enhanced local wall friction over a continuously curved duct, which fits
a cooling channel bending gently over many diameters. They are *not* a model
for a discrete pipe bend, where the loss is set by secondary-flow development
and its decay downstream. Use a fitting loss coefficient for those.
"""

from __future__ import annotations

import warnings
from typing import cast

import numpy as np
from scipy.optimize import brentq

__all__ = [
    "FRICTION_MODELS",
    "FrictionCorrelationWarning",
    "ReDh_laminar",
    "ReDh_turbulent",
    "darcy_friction_factor",
    "curvature_multiplier",
    "friction_pressure_gradient",
    "resolve_friction_model",
    "transition_weight",
]


#: Upper Reynolds number of the laminar branch.
ReDh_laminar = 2300

#: Lower Reynolds number of the fully turbulent branch. Between this and
#: :data:`ReDh_laminar` the native model interpolates linearly.
ReDh_turbulent = 3500

#: Backend names accepted by every dispatching entry point in this module.
FRICTION_MODELS = ("native", "fluids")

#: Ratio between the fully turbulent and critical Reynolds numbers used by
#: the shared ``fluids`` transition policy. This preserves the relative width
#: of pyskyfire's established 2300--3500 transition interval.
_TRANSITION_END_RATIO = ReDh_turbulent / ReDh_laminar

#: Critical Reynolds number used for a straight duct by the ``fluids`` backend.
_FLUIDS_STRAIGHT_RE_CRITICAL = 2040.0


class FrictionCorrelationWarning(UserWarning):
    """A friction correlation is being extrapolated beyond a documented limit."""


def transition_weight(Re, Re_critical):
    """Return the shared laminar-to-turbulent blending weight.

    The weight is zero through ``Re_critical`` and one from
    ``Re_critical * (3500/2300)`` onward. Between those limits it follows the
    cubic smoothstep ``3*x**2 - 2*x**3``, giving continuous values and first
    derivatives where the pure laminar and turbulent branches are joined.

    The fixed critical value ``2040`` is used for straight ducts. Curved ducts
    supply Schmidt's geometry-dependent critical Reynolds number.

    Parameters
    ----------
    Re : float
        Reynolds number.
    Re_critical : float
        Reynolds number at the start of transition.

    Returns
    -------
    float
        Turbulent-branch weight in the closed interval ``[0, 1]``.
    """
    if not np.isfinite(Re) or not np.isfinite(Re_critical):
        raise ValueError("transition Reynolds numbers must be finite")
    if Re_critical <= 0.0:
        raise ValueError("Re_critical must be positive")
    if Re <= Re_critical:
        return 0.0
    Re_turbulent = Re_critical * _TRANSITION_END_RATIO
    if Re >= Re_turbulent:
        return 1.0
    fraction = (Re - Re_critical) / (Re_turbulent - Re_critical)
    return float(fraction**2 * (3.0 - 2.0 * fraction))


def resolve_friction_model(model):
    """Return a validated backend name.

    Parameters
    ----------
    model : str
        Backend name, case-insensitive. One of :data:`FRICTION_MODELS`.

    Returns
    -------
    str
        The normalized (lower-case) backend name.

    Raises
    ------
    ValueError
        If ``model`` is not a recognized backend.
    """
    if not isinstance(model, str):
        raise ValueError(
            f"friction model must be a string, got {type(model).__name__}"
        )
    resolved = model.lower()
    if resolved not in FRICTION_MODELS:
        raise ValueError(
            f"unknown friction model {model!r}; expected one of "
            f"{', '.join(FRICTION_MODELS)}"
        )
    return resolved


# ---------------------------------------------------------------------------
# Native backend
# ---------------------------------------------------------------------------


def _native_laminar(ReDh):
    """Return the laminar Darcy friction factor ``64/Re``."""
    return 64.0 / ReDh


def _native_turbulent(ReDh, Dh, roughness):
    """Return the turbulent Darcy friction factor.

    Uses Petukhov's smooth-wall correlation when ``roughness`` is zero or
    ``None``, and otherwise solves Colebrook-White iteratively.
    """
    if roughness is None or roughness == 0.0:
        # Petukhov (1970), smooth wall.
        return (0.79 * np.log(ReDh) - 1.64) ** (-2)

    relative_roughness = roughness / Dh
    if not np.isfinite(ReDh) or ReDh <= 0.0:
        raise ValueError("native Colebrook solve requires finite positive ReDh")
    if not np.isfinite(relative_roughness) or not 0.0 < relative_roughness < 3.71:
        raise ValueError("native Colebrook solve requires 0 < roughness/Dh < 3.71")

    # Solve for y = 1/sqrt(f). The residual is strictly increasing, negative
    # at zero and positive at the fully rough upper bound. Bracketing avoids
    # negative trial friction factors and an unchecked initial-guess solve.
    def residual(y):
        return y + 2.0 * np.log10(relative_roughness / 3.71 + 2.51 * y / ReDh)

    upper = -2.0 * np.log10(relative_roughness / 3.71)
    inverse_sqrt_f = cast(float, brentq(residual, 0.0, upper, xtol=1e-12))
    return 1.0 / inverse_sqrt_f**2


def _native_darcy(ReDh, Dh, roughness):
    """Composite laminar-turbulent Darcy factor with a linear transition."""
    if ReDh < ReDh_laminar:
        return _native_laminar(ReDh)

    if ReDh < ReDh_turbulent:
        f_laminar = _native_laminar(ReDh)
        f_turbulent = _native_turbulent(ReDh, Dh, roughness)
        return np.interp(
            ReDh,
            [ReDh_laminar, ReDh_turbulent],
            [f_laminar, f_turbulent],
        )

    return _native_turbulent(ReDh, Dh, roughness)


def _curvature_group(Re, length_scale, R_curv):
    """Return ``Re * (length_scale / |R|)^2`` or zero if effectively straight."""
    values = (Re, length_scale, R_curv)
    if not all(np.isfinite(value) for value in values):
        return 0.0
    if Re <= 0.0 or length_scale <= 0.0 or R_curv == 0.0:
        return 0.0
    return float(Re * (length_scale / abs(R_curv)) ** 2)


def _native_curvature(Re, Dh, R_curv, roughness=0.0):
    """RTE/Ito curvature multiplier for the Darcy friction factor.

    The correlation uses the hydraulic radius ``r_h = D_h/4`` and applies only
    for ``Re * (r_h/R)^2 > 6``. It is independent of bend direction.
    """
    group = _curvature_group(Re, 0.25 * Dh, R_curv)
    if group <= 6.0:
        return 1.0
    return group**0.05


# ---------------------------------------------------------------------------
# fluids backend
# ---------------------------------------------------------------------------


def _fluids_darcy(ReDh, Dh, roughness):
    """Return a smoothly transitioned ``fluids`` straight-duct Darcy factor.

    Uses the exact laminar factor through ``Re=2040``, Clamond's
    near-machine-precision Colebrook solver above the shared transition, and
    :func:`transition_weight` to blend the two branches between those limits.
    Calling the branch functions directly avoids the abrupt regime switch in
    the generic ``fluids.friction_factor`` dispatcher.
    """
    import fluids

    relative_roughness = 0.0 if not roughness else roughness / Dh
    weight = transition_weight(ReDh, _FLUIDS_STRAIGHT_RE_CRITICAL)
    if weight == 0.0:
        return float(fluids.friction.friction_laminar(ReDh))
    turbulent = float(fluids.friction.Clamond(ReDh, relative_roughness))
    if weight == 1.0:
        return turbulent
    laminar = float(fluids.friction.friction_laminar(ReDh))
    return (1.0 - weight) * laminar + weight * turbulent


def _fluids_curvature(Re, Dh, R_curv, roughness=0.0):
    """Curvature multiplier from the ``fluids`` helical-coil correlations.

    Returns ``blended_curved / blended_straight``, so it is directly comparable
    with :func:`_native_curvature`. Both factors use the same roughness and the
    same transition policy. The curved transition starts at Schmidt's critical
    Reynolds number rather than the fixed straight-duct value. Schmidt's
    laminar model itself is independent of roughness. The sign of ``R_curv``
    is discarded: helical-coil correlations do not distinguish concave from
    convex bending.
    """
    import fluids

    if not np.isfinite(R_curv) or R_curv == 0.0 or not np.isfinite(Re) or Re <= 0.0:
        return 1.0

    # ``Dc`` is the coil diameter measured centre-to-centre of the tube, which
    # is twice the centerline radius of curvature.
    coil_diameter = 2.0 * abs(float(R_curv))
    # These are the limits documented by fluids for the selected Schmidt
    # correlations, not general limits on curved-pipe flow. Keep messages
    # constant so Python's default filter does not warn at every march station.
    if Dh / coil_diameter >= 0.14:
        warnings.warn(
            "Schmidt transition criterion recommends Dh/(2*abs(R_curv)) < 0.14; "
            "high-curvature result is extrapolated.",
            FrictionCorrelationWarning,
            stacklevel=3,
        )
    if Re <= 100.0:
        warnings.warn(
            "Schmidt laminar curved friction is documented for Re > 100; "
            "low-Reynolds-number result is extrapolated.",
            FrictionCorrelationWarning,
            stacklevel=3,
        )
    if Re > 150_000.0:
        warnings.warn(
            "Schmidt turbulent curved friction is documented up to Re = 150000; "
            "high-Reynolds-number result is extrapolated.",
            FrictionCorrelationWarning,
            stacklevel=3,
        )
    roughness = 0.0 if roughness is None else roughness
    straight = _fluids_darcy(Re, Dh, roughness)
    Re_critical = float(
        fluids.friction.helical_transition_Re_Schmidt(Dh, coil_diameter)
    )
    weight = transition_weight(Re, Re_critical)
    if weight == 0.0:
        curved = float(
            fluids.friction.helical_laminar_fd_Schmidt(
                Re, Dh, coil_diameter,
            )
        )
    else:
        curved_turbulent = float(
            fluids.friction.helical_turbulent_fd_Schmidt(
                Re, Dh, coil_diameter, roughness,
            )
        )
        if weight == 1.0:
            curved = curved_turbulent
        else:
            curved_laminar = float(
                fluids.friction.helical_laminar_fd_Schmidt(
                    Re, Dh, coil_diameter,
                )
            )
            curved = (
                (1.0 - weight) * curved_laminar
                + weight * curved_turbulent
            )
    if straight <= 0.0:
        return 1.0
    return curved / straight


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------

_DARCY_BACKENDS = {"native": _native_darcy, "fluids": _fluids_darcy}
_CURVATURE_BACKENDS = {"native": _native_curvature, "fluids": _fluids_curvature}


def darcy_friction_factor(ReDh, Dh, roughness=0.0, *, model="native"):
    """Return the Darcy friction factor for a straight duct.

    Parameters
    ----------
    ReDh : float
        Reynolds number on the hydraulic diameter.
    Dh : float
        Hydraulic diameter [m].
    roughness : float, optional
        Absolute wall roughness [m]. Zero or ``None`` means hydraulically
        smooth. Default 0.
    model : {'native', 'fluids'}, optional
        Backend to evaluate. Default ``'native'``.

    Returns
    -------
    float
        Darcy friction factor (dimensionless).
    """
    return _DARCY_BACKENDS[resolve_friction_model(model)](ReDh, Dh, roughness)


def curvature_multiplier(ReDh, Dh, R_curv, *, roughness=0.0, model="native"):
    """Return the friction multiplier due to streamwise curvature.

    Parameters
    ----------
    ReDh : float
        Reynolds number on the hydraulic diameter.
    Dh : float
        Hydraulic diameter [m].
    R_curv : float
        Signed radius of curvature of the duct centerline [m]. Use ``np.inf``
        for a straight duct. The native backend ignores the sign; so does the
        ``fluids`` backend.
    roughness : float, optional
        Absolute wall roughness [m]. Zero or ``None`` means smooth. Used in
        both factors of the ``fluids`` curved/straight ratio; ignored by the
        native curvature correlation. Use the same value as for the straight
        factor being multiplied.
    model : {'native', 'fluids'}, optional
        Backend to evaluate. Default ``'native'``.

    Returns
    -------
    float
        Multiplier on the straight-duct Darcy factor. Exactly ``1.0`` for a
        straight duct; native also returns one below its activation threshold.

    Warns
    -----
    FrictionCorrelationWarning
        For the ``fluids`` curved model, when ``ReDh <= 100``,
        ``ReDh > 150000``, or ``Dh/(2*abs(R_curv)) >= 0.14``. The last limit
        is the recommendation for the Schmidt transition criterion; upstream
        does not specify a separate curvature range for the friction fits.
        Extrapolated results are returned without clipping or model switching.

    Notes
    -----
    Limits are from the ``fluids.friction`` documentation for
    ``helical_laminar_fd_Schmidt``, ``helical_turbulent_fd_Schmidt`` and
    ``helical_transition_Re_Schmidt``. These warnings concern the coil models,
    not straight-pipe Colebrook numerical accuracy or two-phase validity.
    """
    return _CURVATURE_BACKENDS[resolve_friction_model(model)](
        ReDh, Dh, R_curv, roughness,
    )


def friction_pressure_gradient(
    ReDh,
    Dh,
    rho,
    velocity,
    *,
    roughness=0.0,
    R_curv=np.inf,
    curvature_correction=True,
    model="native",
):
    """Return the irreversible pressure gradient ``dp/ds`` [Pa m⁻¹].

    This is the friction term of the one-dimensional momentum equation only.
    The reversible acceleration term ``-rho u du/ds`` depends on the state at
    both ends of a segment rather than on a local gradient, so the caller
    applies it while marching.

    Parameters
    ----------
    ReDh : float
        Reynolds number on the hydraulic diameter.
    Dh : float
        Hydraulic diameter [m].
    rho : float
        Bulk density [kg m⁻³].
    velocity : float
        Bulk velocity [m s⁻¹].
    roughness : float, optional
        Absolute wall roughness [m]. Default 0.
    R_curv : float, optional
        Signed radius of curvature [m]. Default ``np.inf`` (straight).
    curvature_correction : bool, optional
        Apply :func:`curvature_multiplier`. Default True.
    model : {'native', 'fluids'}, optional
        Backend to evaluate. Default ``'native'``.

    Returns
    -------
    float
        Pressure gradient with respect to arc length [Pa m⁻¹]. Negative for
        flow in the positive arc-length direction.
    """
    resolved = resolve_friction_model(model)
    friction = darcy_friction_factor(ReDh, Dh, roughness, model=resolved)
    if curvature_correction:
        friction *= curvature_multiplier(
            ReDh, Dh, R_curv, roughness=roughness, model=resolved,
        )
    return -friction / Dh * rho * velocity**2 / 2.0
