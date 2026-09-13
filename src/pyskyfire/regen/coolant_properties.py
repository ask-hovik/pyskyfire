"""Coolant film properties and two-phase modelling limits for the regenerative march.

Why the march is on enthalpy, not temperature
---------------------------------------------

The obvious way to advance the coolant along a channel is
:math:`\\Delta T = Q / (\\dot m\\, c_p)`, and that is what
:meth:`~pyskyfire.regen.coupled_solver.CoupledHeatExchangerPhysics.coolant_temperature_rate`
does. It is only valid while the coolant is single-phase. Once the bulk
reaches saturation, added heat goes into latent heat at constant temperature,
and temperature stops being a state variable at all: inside the dome ``T`` and
``p`` are not independent, so ``(T, p)`` no longer identifies a state.

That failure is silent rather than loud. CoolProp answers a ``(T, p)`` query
just inside the dome by picking a single-phase branch, so a temperature march
that steps past :math:`T_\\text{sat}` lands on *superheated vapour*. For
ethanol at 20 bar, 0.2 K either side of saturation:

===================  ===========================
``T_sat - 0.1 K``    :math:`\\rho = 596.9` kg m⁻³
``T_sat + 0.1 K``    :math:`\\rho = 31.4` kg m⁻³
===================  ===========================

A 19x density jump, and with it a 19x jump in ``rho u^2`` and therefore in the
pressure gradient. No exception is raised; the march simply continues on the
wrong branch and returns a plausible-looking profile.

Specific enthalpy has none of these problems. It is monotonic in the heat
added, it stays a valid state variable through the dome, and ``(h, p)``
identifies the state everywhere -- subcooled, saturated, superheated, or
supercritical. So the march integrates

.. math:: h_{i+1} = h_i + Q_\\text{cold} / \\dot m_\\text{channel}

and :meth:`~pyskyfire.thermo.state.FluidStateResolver.from_hp` resolves ``(T, quality, rho)`` afterwards.

What this module does *not* model
---------------------------------

**Two-phase heat transfer.** In the dome the coolant-side coefficient is still
the single-phase Colburn correlation, evaluated with *saturated liquid*
properties. That is a placeholder in both directions: it misses the flow-boiling
enhancement (which would raise it, by 2-10x) and it misses dryout (which would
collapse it). Treat the quality profile as telling you *where* the coolant
boils, not how hot the wall gets there. Nucleate boiling and CHF are not
modelled at all.

**Subcooled boiling.** Where the bulk is liquid but the wall is above
:math:`T_\\text{sat}`, real channels boil at the wall. That is not modelled
either. What *is* handled is the property-lookup artefact it used to cause: the
Colburn film temperature :math:`(T_\\text{cool} + T_\\text{cw})/2` would cross
saturation and CoolProp would return *vapour* properties for what is physically
superheated liquid, putting a spurious cliff in ``h_cold``. Film properties are
now clamped to saturated liquid instead (see
:func:`film_properties`). Still an approximation, but no longer a
discontinuity, and wrong in the conservative direction.

**Two-phase pressure drop.** Density in the dome is the homogeneous (no-slip)
value :math:`1/(x/\\rho_g + (1-x)/\\rho_f)`. This captures the large
acceleration effect but not the two-phase friction multiplier, so channel
:math:`\\Delta p` is under-predicted once quality is appreciable.

Where the enthalpy march is unavailable
---------------------------------------

:func:`~pyskyfire.thermo.state.probe_fluid` decides per fluid, once, at solver start:

* **Mixtures** (``HEOS::A[..]&B[..]``) fall back to the temperature march. Not
  because the physics is different but because it is computationally
  impossible: a mixture ``(h, p)`` flash costs ~374 ms against ~57 us for a
  pure fluid, a factor of 6500, which would put a 100-node run into the hours.
  Mixtures also have no ``pcrit`` through ``Props1SI`` and no surface tension,
  so most of the two-phase machinery has nothing to stand on. Revisit when
  tabulated backends land.
* **Incompressible backends** (``INCOMP::``) have no saturation data.
* **Supercritical operation** (:math:`p > p_\\text{crit}`) keeps the enthalpy
  march -- it is still the better-conditioned variable -- but there is no dome,
  so quality is undefined and reported as ``nan``. This is the normal state of
  affairs for hydrogen (:math:`p_\\text{crit}` = 13.0 bar) and for kerosene-like
  fluids (n-dodecane, 18.2 bar) at any realistic channel pressure. Heat-transfer
  deterioration near the pseudo-critical line is not modelled.
"""

from __future__ import annotations

__all__ = ["film_properties"]


# Temperature step used to back away from CoolProp's saturation guard band [K].
_SATURATION_NUDGE = 1.0e-3


def film_properties(resolver, T_film: float, p: float, getters) -> tuple:
    """Coolant film properties for the Colburn correlation.

    Parameters
    ----------
    resolver : pyskyfire.thermo.state.FluidStateResolver
        Fluid state resolver providing cached saturation properties.
    T_film : float
        Film temperature :math:`(T_\\text{cool} + T_\\text{cw})/2` [K].
    p : float
        Coolant pressure [Pa].
    getters : sequence of callable
        Property getters taking ``(T, p)``, in the order
        ``(conductivity, cp, viscosity)``.

    Returns
    -------
    tuple
        ``(k, cp, mu)`` at the film condition.

    Notes
    -----
    The film temperature routinely exceeds :math:`T_\\text{sat}` while the
    bulk is still subcooled liquid, because the wall is hot. A raw ``(T, p)``
    lookup there returns *vapour* properties, which is physically wrong --
    the fluid at the wall is superheated liquid on the verge of nucleating,
    not vapour -- and puts a cliff in ``h_cold`` that the local wall solve
    can fail to cross. Clamping to saturated liquid removes the cliff. It is
    still not a boiling model: the real coefficient in that region is higher
    than this returns.
    """
    k_get, cp_get, mu_get = getters
    sat = resolver.saturation(p)

    if sat is not None and T_film >= sat.T_sat:
        return sat.k_f, sat.cp_f, sat.mu_f

    return (
        _off_saturation(k_get, T_film, p),
        _off_saturation(cp_get, T_film, p),
        _off_saturation(mu_get, T_film, p),
    )


def _off_saturation(getter, T, p):
    """Evaluate a property, stepping off CoolProp's saturation guard band.

    CoolProp refuses a ``(T, p)`` pair within 1e-4 % of the saturation line
    because both branches are valid there and it will not choose. Retry a
    millikelvin to the liquid side. This is a crash guard, not a model.
    """
    try:
        return getter(T, p)
    except ValueError:
        return getter(T - _SATURATION_NUDGE, p)
