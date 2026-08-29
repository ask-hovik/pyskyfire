"""Efficiency corrections from ideal (ODE) to delivered engine performance.

CEA, and therefore :class:`~pyskyfire.skycea.aerothermodynamics.Aerothermodynamics`,
solves a one-dimensional equilibrium expansion. Real nozzle flow loses
performance to two-dimensional divergence, the wall boundary layer and
finite-rate chemistry, and real injectors do not release the full ideal
combustion energy. The standard engineering treatment is to keep the ideal
solution and multiply it by efficiencies obtained from higher-fidelity codes:

.. math::

    (I_{sp})_{actual} = (I_{sp})_{ideal}\\, \\eta_{cf}\\, \\eta_{c^*}

``eta_cf`` is a *thrust-coefficient* efficiency and carries everything
downstream of the throat. ``eta_cstar`` is a *characteristic-velocity*
efficiency and carries everything upstream of it. This is the lumped form of
the JANNAF simplified methodology: formal JANNAF resolves the nozzle side into
separate kinetic (ODE->ODK), divergence (ODK->TDK) and boundary-layer terms,
whereas a single ``eta_cf`` taken from a TDK/ODE thrust-coefficient ratio rolls
all three into one number.

The reference implementation of the constants here is Binder, M., Tomsik, T.,
and Veres, J. P., "RL10A-3-3A Rocket Engine Modeling Project", NASA
TM-107318, 1997 (NTRS 19970010379), which states the relation in appendix E.2
and tabulates :math:`\\eta_{cf}` in table E3.

Notes
-----
The nozzle discharge coefficient is deliberately *not* part of the specific
impulse chain. It rescales the mass flow that a given throat area passes, not
the impulse each unit of propellant delivers, and in the reference above it is
applied separately (sec. 4.2.3.3). :meth:`NozzleLoss.mass_flow` is provided for
that purpose; :meth:`NozzleLoss.apply` leaves mass flow untouched.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Union

import numpy as np

__all__ = [
    "EfficiencyTable",
    "NozzleLoss",
    "DeliveredPerformance",
    "RL10A_3_3A_ETA_CF",
    "RL10A_3_3A_LOSS",
]

# Matches pyskyfire.skycea.aerothermodynamics.G0. Duplicated rather than
# imported so this module stays free of the CEA dependency; test_nozzle_loss
# asserts the two stay equal.
G0 = 9.81
SEA_LEVEL_PRESSURE = 1.01325e5
PSI_TO_PA = 6894.757

Efficiency = Union[float, Callable[[float, float], float], "EfficiencyTable"]


class EfficiencyTable:
    """Efficiency tabulated against mixture ratio and chamber pressure.

    Parameters
    ----------
    mixture_ratio : array_like
        Mixture-ratio axis, strictly increasing [-].
    chamber_pressure : array_like
        Chamber-pressure axis, strictly increasing [Pa].
    values : array_like
        Efficiency values, shape ``(len(mixture_ratio), len(chamber_pressure))``.
    log_pressure : bool, optional
        Interpolate against ``log(p_c)`` rather than ``p_c``. Default ``True``,
        which is almost always what you want: published tables span decades of
        chamber pressure on a handful of columns, and efficiency tracks the log
        far more closely than the linear value.

    Notes
    -----
    Queries outside the tabulated ranges are **clamped** to the edge values
    rather than extrapolated. Efficiency tables of this kind are fits over a
    stated envelope and go badly wrong outside it -- the reference report notes
    that the vendor tables it compared against can return negative specific
    impulse when pushed past their intended range.
    """

    def __init__(self, mixture_ratio, chamber_pressure, values, *,
                 log_pressure: bool = True):
        self.mixture_ratio = np.asarray(mixture_ratio, dtype=float)
        self.chamber_pressure = np.asarray(chamber_pressure, dtype=float)
        self.values = np.asarray(values, dtype=float)
        self.log_pressure = bool(log_pressure)

        if self.mixture_ratio.ndim != 1 or self.mixture_ratio.size < 2:
            raise ValueError("mixture_ratio must be 1-D with at least two entries")
        if self.chamber_pressure.ndim != 1 or self.chamber_pressure.size < 2:
            raise ValueError("chamber_pressure must be 1-D with at least two entries")
        if np.any(np.diff(self.mixture_ratio) <= 0.0):
            raise ValueError("mixture_ratio must be strictly increasing")
        if np.any(np.diff(self.chamber_pressure) <= 0.0):
            raise ValueError("chamber_pressure must be strictly increasing")
        if self.values.shape != (self.mixture_ratio.size, self.chamber_pressure.size):
            raise ValueError(
                f"values must have shape {(self.mixture_ratio.size, self.chamber_pressure.size)}, "
                f"got {self.values.shape}"
            )
        if not np.all(np.isfinite(self.values)):
            raise ValueError("values must all be finite")
        if np.any(self.values <= 0.0):
            raise ValueError("values must all be positive")
        if self.log_pressure and self.chamber_pressure[0] <= 0.0:
            raise ValueError("log_pressure requires a positive chamber_pressure axis")

        self._pressure_axis = (
            np.log(self.chamber_pressure) if self.log_pressure else self.chamber_pressure
        )

    def __call__(self, p_c, MR):
        """Interpolate the table at ``p_c`` [Pa] and mixture ratio ``MR`` [-].

        Scalars return a float; array inputs broadcast against each other and
        return an array of the broadcast shape.
        """
        pressure = np.asarray(p_c, dtype=float)
        ratio = np.asarray(MR, dtype=float)
        shape = np.broadcast_shapes(pressure.shape, ratio.shape)
        pressure = np.broadcast_to(pressure, shape).ravel()
        ratio = np.broadcast_to(ratio, shape).ravel()

        if self.log_pressure:
            if np.any(pressure <= 0.0):
                raise ValueError("p_c must be positive when log_pressure is set")
            pressure = np.log(pressure)

        # np.interp clamps to the endpoints, which is the behaviour we want.
        rows = np.stack(
            [np.interp(pressure, self._pressure_axis, row) for row in self.values]
        )

        upper = np.clip(
            np.searchsorted(self.mixture_ratio, ratio, side="right"),
            1,
            self.mixture_ratio.size - 1,
        )
        lower = upper - 1
        span = self.mixture_ratio[upper] - self.mixture_ratio[lower]
        weight = np.clip((ratio - self.mixture_ratio[lower]) / span, 0.0, 1.0)

        columns = np.arange(ratio.size)
        low = rows[lower, columns]
        high = rows[upper, columns]
        result = (low + weight * (high - low)).reshape(shape)

        return float(result) if result.ndim == 0 else result


def _evaluate(efficiency: Efficiency, p_c: float, MR: float) -> float:
    """Resolve a constant, callable or table efficiency to a number."""
    if callable(efficiency):
        value = float(efficiency(p_c, MR))
    else:
        value = float(efficiency)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(
            f"efficiency must evaluate to a finite positive number, got {value!r}"
        )
    return value


@dataclass(frozen=True)
class DeliveredPerformance:
    """Performance after efficiencies, at unchanged propellant flow.

    Attributes
    ----------
    eta_cf, eta_cstar, eta_Isp : float
        The thrust-coefficient, characteristic-velocity and combined specific
        impulse efficiencies actually applied [-].
    c_star : float
        Delivered characteristic velocity [m s⁻¹].
    CF_vac, CF_amb, CF_SL : float
        Delivered thrust coefficients in vacuum, at the design ambient
        pressure, and at sea level [-].
    Isp_vac, Isp_amb, Isp_SL : float
        Delivered specific impulses [s].
    F_vac, F_amb : float
        Delivered thrust at the unchanged mass flow [N].
    mdot : float
        Propellant mass flow, echoed unchanged [kg s⁻¹].
    """

    eta_cf: float
    eta_cstar: float
    eta_Isp: float
    c_star: float
    CF_vac: float
    CF_amb: float
    CF_SL: float
    Isp_vac: float
    Isp_amb: float
    Isp_SL: float
    F_vac: float
    F_amb: float
    mdot: float


@dataclass(frozen=True)
class NozzleLoss:
    """Multiplicative efficiency model taking ideal performance to delivered.

    Parameters
    ----------
    eta_cf : float, callable or EfficiencyTable
        Thrust-coefficient efficiency [-]. A callable is invoked as
        ``eta_cf(p_c, MR)`` with ``p_c`` in Pa.
    eta_cstar : float, callable or EfficiencyTable, optional
        Characteristic-velocity efficiency [-]. Default 1.0.
    Cd : float, optional
        Nozzle discharge coefficient [-]. Used only by :meth:`mass_flow`; it is
        not part of the specific-impulse chain. Default 1.0.
    name : str, optional
        Free-text label, for reports.

    Examples
    --------
    >>> loss = NozzleLoss(eta_cf=0.9479, eta_cstar=0.9892)
    >>> round(loss.specific_impulse(469.11, p_c=3.32e6, MR=5.25), 1)
    439.9
    """

    eta_cf: Efficiency
    eta_cstar: Efficiency = 1.0
    Cd: float = 1.0
    name: str | None = None

    def __post_init__(self):
        # Constants can be checked now; the dummy conditions are ignored for
        # them. Callables and tables can only be checked when they are called.
        if not callable(self.eta_cf):
            _evaluate(self.eta_cf, 0.0, 0.0)
        if not callable(self.eta_cstar):
            _evaluate(self.eta_cstar, 0.0, 0.0)
        if not np.isfinite(self.Cd) or self.Cd <= 0.0:
            raise ValueError(f"Cd must be finite and positive, got {self.Cd!r}")

    # ------------------------------------------------------------------
    # efficiencies
    # ------------------------------------------------------------------
    def thrust_coefficient_efficiency(self, p_c: float, MR: float) -> float:
        """Thrust-coefficient efficiency at ``p_c`` [Pa] and ``MR`` [-]."""
        return _evaluate(self.eta_cf, p_c, MR)

    def c_star_efficiency(self, p_c: float, MR: float) -> float:
        """Characteristic-velocity efficiency at ``p_c`` [Pa] and ``MR`` [-]."""
        return _evaluate(self.eta_cstar, p_c, MR)

    def overall_efficiency(self, p_c: float, MR: float) -> float:
        """Combined specific-impulse efficiency ``eta_cf * eta_cstar`` [-]."""
        return (
            self.thrust_coefficient_efficiency(p_c, MR)
            * self.c_star_efficiency(p_c, MR)
        )

    # ------------------------------------------------------------------
    # scalar corrections
    # ------------------------------------------------------------------
    def specific_impulse(self, Isp_ideal: float, p_c: float, MR: float) -> float:
        """Delivered specific impulse [s] from the ideal ODE value [s]."""
        return float(Isp_ideal) * self.overall_efficiency(p_c, MR)

    def characteristic_velocity(self, c_star_ideal: float, p_c: float,
                                MR: float) -> float:
        """Delivered characteristic velocity [m s⁻¹] from the ideal value."""
        return float(c_star_ideal) * self.c_star_efficiency(p_c, MR)

    def thrust_coefficient(self, CF_ideal: float, p_c: float, MR: float) -> float:
        """Delivered *vacuum* thrust coefficient [-] from the ideal value.

        Ambient operation is handled by :meth:`apply`, which subtracts the
        pressure term after scaling. Do not apply ``eta_cf`` to an
        ambient-corrected coefficient: the ``p_amb/p_c * eps`` term is exact
        geometry and pressure, and carries no nozzle loss.
        """
        return float(CF_ideal) * self.thrust_coefficient_efficiency(p_c, MR)

    def mass_flow(self, p_c: float, A_t: float, c_star_ideal: float,
                  MR: float) -> float:
        """Mass flow [kg s⁻¹] through a geometric throat area ``A_t`` [m²].

        Applies both the discharge coefficient and the characteristic-velocity
        efficiency:

        .. math::

            \\dot{m} = \\frac{C_d\\, p_c\\, A_t}{\\eta_{c^*}\\, c^*_{ideal}}
        """
        return (
            self.Cd * float(p_c) * float(A_t)
            / self.characteristic_velocity(c_star_ideal, p_c, MR)
        )

    # ------------------------------------------------------------------
    # whole-engine
    # ------------------------------------------------------------------
    def apply(self, aero, p_amb: float | None = None) -> DeliveredPerformance:
        """Degrade an :class:`Aerothermodynamics` design point.

        Parameters
        ----------
        aero : object
            Anything exposing the ideal design point: ``p_c``, ``MR``,
            ``c_star``, ``CF_vac``, ``eps``, ``mdot`` and ``p_amb``. Specific
            impulses are rebuilt from ``CF`` and ``c_star`` rather than read
            off ``aero``, so that the ambient pressure term stays unscaled.
        p_amb : float, optional
            Ambient pressure [Pa] for the ``*_amb`` outputs. Defaults to the
            value stored on ``aero``.

        Returns
        -------
        DeliveredPerformance

        Notes
        -----
        Mass flow is passed through unchanged -- the engine burns what it
        burns, and these efficiencies say how much impulse that flow delivers.
        Use :meth:`mass_flow` if you instead want the flow a throat area
        passes.
        """
        p_c = float(aero.p_c)
        MR = float(aero.MR)
        eps = float(aero.eps)
        mdot = float(aero.mdot)
        ambient = float(aero.p_amb if p_amb is None else p_amb)

        eta_cf = self.thrust_coefficient_efficiency(p_c, MR)
        eta_cstar = self.c_star_efficiency(p_c, MR)

        c_star = float(aero.c_star) * eta_cstar
        CF_vac = float(aero.CF_vac) * eta_cf
        CF_amb = CF_vac - ambient / p_c * eps
        CF_SL = CF_vac - SEA_LEVEL_PRESSURE / p_c * eps

        return DeliveredPerformance(
            eta_cf=eta_cf,
            eta_cstar=eta_cstar,
            eta_Isp=eta_cf * eta_cstar,
            c_star=c_star,
            CF_vac=CF_vac,
            CF_amb=CF_amb,
            CF_SL=CF_SL,
            Isp_vac=CF_vac * c_star / G0,
            Isp_amb=CF_amb * c_star / G0,
            Isp_SL=CF_SL * c_star / G0,
            F_vac=mdot * CF_vac * c_star,
            F_amb=mdot * CF_amb * c_star,
            mdot=mdot,
        )


# ----------------------------------------------------------------------
# RL10A-3-3A, Binder et al. NASA TM-107318 (NTRS 19970010379)
# ----------------------------------------------------------------------

#: Thrust-coefficient efficiency from table E3, p. 140. Generated by taking the
#: ratio of TDK to ODE thrust coefficient, so it carries two-dimensional
#: divergence, boundary-layer and finite-rate chemistry losses together.
RL10A_3_3A_ETA_CF = EfficiencyTable(
    mixture_ratio=[1.0, 4.0, 6.0, 8.0, 10.0, 20.0, 50.0],
    chamber_pressure=[5.0 * PSI_TO_PA, 50.0 * PSI_TO_PA,
                      150.0 * PSI_TO_PA, 500.0 * PSI_TO_PA],
    values=[
        [0.9381, 0.9534, 0.9567, 0.9593],   # O/F = 1
        [0.8860, 0.9269, 0.9419, 0.9517],   # O/F = 4
        [0.8508, 0.9024, 0.9282, 0.9469],   # O/F = 6
        [0.8463, 0.8950, 0.9193, 0.9414],   # O/F = 8
        [0.8507, 0.8978, 0.9214, 0.9419],   # O/F = 10
        [0.8853, 0.9226, 0.9374, 0.9486],   # O/F = 20
        [0.9568, 0.9537, 0.9568, 0.9596],   # O/F = 50
    ],
)

#: Full RL10A-3-3A loss set. ``eta_cstar`` is the design-point value from table
#: 2.5.1; the report shows it varying with mixture ratio and chamber pressure
#: (fig. 4.2.8) but only tabulates the single number, so it is held constant
#: here. ``Cd`` is the constant 0.975 adopted for the engine system model in
#: sec. 4.2.3.3, in preference to the 0.982 implied by the ratio of P&W's
#: effective flow area to the physical throat.
RL10A_3_3A_LOSS = NozzleLoss(
    eta_cf=RL10A_3_3A_ETA_CF,
    eta_cstar=0.9892,
    Cd=0.975,
    name="RL10A-3-3A (Binder et al., NASA TM-107318)",
)
