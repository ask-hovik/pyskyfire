"""Orifices: short restrictions described by a flow area and a discharge coefficient.

An :class:`Orifice` is a hole that discharges into a large volume: an injector
orifice into the chamber, a valve seat into the line behind it, a bleed into a
cavity. Its jet's kinetic energy is not recovered, so the whole drop from the
inlet to the discharge pressure is lost. An :class:`OrificePlate` is a thin
round hole sitting in a pipe: the flow approaches it at the pipe velocity, and
part of the drop recovers downstream as the jet re-expands to the pipe bore.
The plate follows ISO 5167; see its notes.

An orifice's ideal flow expands isentropically from the inlet stagnation state
to the throat pressure through the fluid's
:class:`~pyskyfire.thermo.eos.EquationOfState`:

.. math::

    \\dot m = C_d A \\rho_t \\sqrt{2 (h_1 - h_t)}

For a constant-density fluid this is the familiar
:math:`C_d A \\sqrt{2 \\rho \\Delta p}`. For a gas it is the isentropic nozzle
flow, and it chokes where the ideal mass flux peaks.
For a liquid expanded below its vapour pressure, a real equation of state
gives the homogeneous-equilibrium two-phase flow, whose flux peaks close to
the vapour pressure: the flow then passes about
:math:`C_d A \\sqrt{2 \\rho (p_1 - p_{sat})}` whatever the discharge pressure,
which is how a cavitating orifice behaves. Metastable (delayed) flashing is
not modelled; the homogeneous model is known to under-predict flow near
saturation.

The discharge coefficient is always supplied, except on a plate, where the
Reader-Harris/Gallagher correlation of ISO 5167-2 is the default.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import warnings

import CoolProp.CoolProp as CP
import fluids
from scipy.optimize import brentq

from pyskyfire.thermo.eos import CoolPropEOS, EquationOfState, as_equation_of_state
from pyskyfire.thermo.state import FluidState

__all__ = [
    "Orifice",
    "OrificePlate",
    "OrificeResult",
    "OrificePlateResult",
    "OrificeCapacityError",
    "OrificeCorrelationWarning",
    "required_area",
]

_GOLDEN = (math.sqrt(5.0) - 1.0) / 2.0
#: Lowest throat pressure searched, as a fraction of the inlet pressure.
_FLOOR_FRACTION = 1.0e-3
#: Pressure taps the Reader-Harris/Gallagher correlation is defined for.
_TAPS = ("corner", "flange", "D", "D/2")


class OrificeCapacityError(ValueError):
    """The inlet state cannot drive the requested mass flow through the orifice.

    Either the flow chokes below the requested mass flow, or passing it would
    need a discharge pressure below the searched floor. ``max_mass_flow`` is
    the largest flow found [kg/s].
    """

    def __init__(self, message: str, max_mass_flow: float):
        super().__init__(message)
        self.max_mass_flow = max_mass_flow


class OrificeCorrelationWarning(RuntimeWarning):
    """A discharge-coefficient correlation is used outside its documented range."""


@dataclass(frozen=True)
class OrificeResult:
    """States and diagnostics of the flow through an orifice.

    ``throat_state`` is the ideal isentropic state at the vena contracta and
    ``velocity`` the ideal jet velocity there [m/s]. The real jet is slower by
    the velocity coefficient, which the discharge coefficient includes but
    does not separate out. ``outlet_state`` is at the discharge pressure with
    the inlet enthalpy: the jet's kinetic energy has been dissipated.
    ``choked`` is True when the discharge pressure lies below the pressure at
    which the ideal mass flux peaks, so lowering it further adds no flow.
    """

    inlet_state: FluidState
    throat_state: FluidState
    outlet_state: FluidState
    mass_flow: float
    discharge_coefficient: float
    velocity: float
    choked: bool

    @property
    def pressure_loss(self) -> float:
        """Irreversible pressure loss from inlet to outlet [Pa]."""
        return self.inlet_state.p - self.outlet_state.p


@dataclass(frozen=True)
class OrificePlateResult(OrificeResult):
    """An :class:`OrificeResult` with the plate's metering quantities.

    ``throat_state`` is at the downstream tap, and
    :attr:`differential_pressure` is the tap-to-tap drop a flow meter reads.
    ``outlet_state`` is where the flow has recovered to the pipe bore, so
    :attr:`pressure_loss` is the permanent loss, which is smaller.
    """

    beta: float
    reynolds: float

    @property
    def differential_pressure(self) -> float:
        """Drop between the upstream and downstream taps [Pa]."""
        return self.inlet_state.p - self.throat_state.p


# --------------------------------------------------------------------------
# Flow model shared by orifices and plates
# --------------------------------------------------------------------------

def _inlet(eos: EquationOfState, state: FluidState) -> FluidState:
    """The inlet as this EOS knows it.

    A state carrying an enthalpy is taken as it is, so two-phase inlets can be
    given; the enthalpy must then be on this EOS's datum. A state without one
    is resolved from its pressure and temperature.
    """
    if not isinstance(state, FluidState):
        raise TypeError("inlet_state must be a thermo.FluidState")
    if not (math.isfinite(state.p) and state.p > 0.0):
        raise ValueError("inlet pressure must be finite and positive")
    if math.isfinite(state.h):
        if math.isfinite(state.rho):
            return state
        return eos.state_ph(state.p, state.h)
    if not (math.isfinite(state.T) and state.T > 0.0):
        raise ValueError("inlet state needs an enthalpy or a finite positive temperature")
    return eos.state_pT(state.p, state.T)


def _ideal_flux(eos, inlet, p_throat, area_ratio=0.0):
    """Ideal mass flux [kg m⁻² s⁻¹], throat state and velocity at ``p_throat``.

    ``area_ratio`` is throat over upstream area, for the approach velocity.
    """
    throat = eos.isentropic(inlet, p_throat)
    dh = inlet.h - throat.h
    if not dh > 0.0:
        return 0.0, throat, 0.0
    carry = (throat.rho * area_ratio / inlet.rho) ** 2
    velocity = math.sqrt(2.0 * dh / (1.0 - carry))
    return throat.rho * velocity, throat, velocity


def _peak_flux(eos, inlet, p_low, p_high):
    """Golden-section search for the throat pressure of peak ideal flux.

    The ideal flux rises as the throat pressure falls, peaks, and then falls
    again, so it is unimodal in ``[p_low, p_high]``.
    """
    a, b = p_low, p_high
    c = b - _GOLDEN * (b - a)
    d = a + _GOLDEN * (b - a)
    fc = _ideal_flux(eos, inlet, c)[0]
    fd = _ideal_flux(eos, inlet, d)[0]
    while b - a > 1.0e-7 * inlet.p:
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - _GOLDEN * (b - a)
            fc = _ideal_flux(eos, inlet, c)[0]
        else:
            a, c, fc = c, d, fd
            d = a + _GOLDEN * (b - a)
            fd = _ideal_flux(eos, inlet, d)[0]
    p_star = 0.5 * (a + b)
    return p_star, _ideal_flux(eos, inlet, p_star)[0]


def _throat_at_discharge(eos, inlet, p_out):
    """Throat pressure for a discharge pressure, and whether the flow is choked."""
    step = 1.0e-3 * (inlet.p - p_out)
    at_out = _ideal_flux(eos, inlet, p_out)[0]
    above = _ideal_flux(eos, inlet, p_out + step)[0]
    if at_out >= above:
        # Still rising as the pressure falls to p_out: the peak is at or below it.
        return p_out, False
    p_star, _ = _peak_flux(eos, inlet, p_out, inlet.p)
    return p_star, True


def _throat_for_flux(eos, inlet, target, what):
    """Throat pressure at which the ideal flux equals ``target``, unchoked."""
    p0 = inlet.p
    floor = _FLOOR_FRACTION * p0
    dp = max(target**2 / (2.0 * inlet.rho), 1.0e-9 * p0)
    history = [(p0, 0.0)]
    best = 0.0
    while True:
        # Never more than halve the throat pressure in one step, so a peak
        # is seen as a falling flux before the search reaches pressures the
        # equation of state cannot expand to.
        p = max(p0 - dp, 0.5 * history[-1][0], floor)
        try:
            flux = _ideal_flux(eos, inlet, p)[0]
        except ValueError as exc:
            raise OrificeCapacityError(
                f"{what} cannot be expanded below {p:g} Pa to pass the flow ({exc})", best
            ) from exc
        best = max(best, flux)
        if flux >= target:
            p_hi = history[-1][0]
            return brentq(lambda x: _ideal_flux(eos, inlet, x)[0] - target,
                          p, p_hi, xtol=1.0e-10 * p0, rtol=1.0e-12)
        if flux < history[-1][1]:
            # Past the peak: it lies between this point and the one two back.
            upper = history[-2][0] if len(history) > 1 else p0
            p_star, peak = _peak_flux(eos, inlet, p, upper)
            if peak < target:
                raise OrificeCapacityError(
                    f"{what} chokes at a throat pressure of {p_star:g} Pa and cannot pass "
                    f"the requested flow; the ideal flux peaks at {peak:g} kg/m²s "
                    f"against {target:g} kg/m²s required", peak)
            return brentq(lambda x: _ideal_flux(eos, inlet, x)[0] - target,
                          p_star, p0, xtol=1.0e-10 * p0, rtol=1.0e-12)
        if p <= floor:
            raise OrificeCapacityError(
                f"{what} would need a discharge pressure below {floor:g} Pa "
                f"to pass the requested flow", best)
        history.append((p, flux))
        dp = 2.0 * (p0 - p)


def _positive(value, label):
    value = float(value)
    if not (math.isfinite(value) and value > 0.0):
        raise ValueError(f"{label} must be finite and positive")
    return value


def _discharge_coefficient(cd):
    cd = float(cd)
    if not (math.isfinite(cd) and 0.0 < cd <= 1.0):
        raise ValueError("cd must lie in (0, 1]")
    return cd


def _mass_flow(mass_flow):
    mass_flow = float(mass_flow)
    if not (math.isfinite(mass_flow) and mass_flow > 0.0):
        raise ValueError("mass_flow must be finite and positive")
    return mass_flow


def _discharge_state(eos, p, h, fallback_T):
    state = eos.state_ph(p, h)
    if not math.isfinite(state.T):
        # An EOS without cp cannot recover T from h; an adiabatic,
        # constant-density throttle keeps it.
        state = FluidState(T=fallback_T, p=state.p, h=state.h, quality=state.quality,
                           rho=state.rho, phase=state.phase, T_sat=state.T_sat)
    return state


# --------------------------------------------------------------------------
# Orifice
# --------------------------------------------------------------------------

class Orifice:
    """A short restriction discharging into a large volume.

    Parameters
    ----------
    area : float
        Geometric flow area [m²].
    cd : float
        Discharge coefficient, in (0, 1]. It relates the real flow to the
        ideal isentropic flow through ``area`` and is supplied, not computed:
        about 0.6 for a sharp-edged hole, 0.8 and above for a rounded or long
        one, and much lower for a swirl element whose air core blocks the exit.
    hydraulic_diameter : float, optional
        Hydraulic diameter of the flow area [m]. The constructors set it.
    length : float, optional
        Length of the passage [m]. It is geometry only here; it does not
        change ``cd``.

    Notes
    -----
    Use :meth:`circular`, :meth:`annular` or :meth:`slot` for those shapes.
    The generic constructor takes any area, so an effective area quoted as
    ``Cd·A`` can be given as ``Orifice(area=CdA, cd=1.0)``.
    """

    def __init__(self, area, *, cd, hydraulic_diameter=None, length=0.0):
        self.area = _positive(area, "area")
        self.cd = _discharge_coefficient(cd)
        self.hydraulic_diameter = (None if hydraulic_diameter is None
                                   else _positive(hydraulic_diameter, "hydraulic_diameter"))
        length = float(length)
        if not (math.isfinite(length) and length >= 0.0):
            raise ValueError("length must be finite and non-negative")
        self.length = length

    def __repr__(self):
        return f"{type(self).__name__}(area={self.area:g}, cd={self.cd:g})"

    @classmethod
    def circular(cls, diameter, *, cd, length=0.0):
        """A round hole of ``diameter`` [m]."""
        diameter = _positive(diameter, "diameter")
        return cls(math.pi * diameter**2 / 4.0, cd=cd, hydraulic_diameter=diameter, length=length)

    @classmethod
    def annular(cls, inner_diameter, outer_diameter, *, cd, length=0.0):
        """The annulus between two concentric diameters [m], as around a coaxial post."""
        inner = _positive(inner_diameter, "inner_diameter")
        outer = _positive(outer_diameter, "outer_diameter")
        if not outer > inner:
            raise ValueError("outer_diameter must exceed inner_diameter")
        return cls(math.pi * (outer**2 - inner**2) / 4.0, cd=cd,
                   hydraulic_diameter=outer - inner, length=length)

    @classmethod
    def slot(cls, width, height, *, cd, length=0.0):
        """A rectangular slot of ``width`` by ``height`` [m]."""
        width = _positive(width, "width")
        height = _positive(height, "height")
        return cls(width * height, cd=cd, hydraulic_diameter=2.0 * width * height / (width + height),
                   length=length)

    @classmethod
    def sized(cls, *, mass_flow, inlet_state, dp, fluid, cd, length=0.0):
        """The round hole that passes ``mass_flow`` with a pressure drop ``dp``.

        See :func:`required_area`.
        """
        area = required_area(mass_flow=mass_flow, inlet_state=inlet_state, dp=dp,
                             fluid=fluid, cd=cd)
        return cls.circular(math.sqrt(4.0 * area / math.pi), cd=cd, length=length)

    @property
    def effective_area(self) -> float:
        """Discharge coefficient times area, ``Cd·A`` [m²]."""
        return self.cd * self.area

    @property
    def length_to_diameter(self) -> float:
        """Passage length over hydraulic diameter, ``L/D``. NaN without a diameter."""
        if self.hydraulic_diameter is None:
            return math.nan
        return self.length / self.hydraulic_diameter

    def _result(self, eos, inlet, p_throat, p_out, choked):
        flux, throat, velocity = _ideal_flux(eos, inlet, p_throat)
        return OrificeResult(
            inlet_state=inlet,
            throat_state=throat,
            outlet_state=_discharge_state(eos, p_out, inlet.h, inlet.T),
            mass_flow=self.cd * self.area * flux,
            discharge_coefficient=self.cd,
            velocity=velocity,
            choked=choked,
        )

    def solve(self, *, inlet_state: FluidState, mass_flow: float, fluid) -> OrificeResult:
        """The discharge pressure that drives ``mass_flow`` through the orifice.

        Parameters
        ----------
        inlet_state : FluidState
            Upstream stagnation state. Its enthalpy is used when finite, so
            it must then be on ``fluid``'s datum; otherwise the state is
            resolved from its pressure and temperature.
        mass_flow : float
            Mass flow [kg/s].
        fluid : str or Fluid or EquationOfState
            The fluid, as :func:`~pyskyfire.thermo.eos.as_equation_of_state`
            reads it. :class:`~pyskyfire.thermo.eos.Incompressible` gives the
            constant-density orifice equation.

        Raises
        ------
        OrificeCapacityError
            When the flow would choke below ``mass_flow``.
        """
        eos = as_equation_of_state(fluid)
        inlet = _inlet(eos, inlet_state)
        target = _mass_flow(mass_flow) / (self.cd * self.area)
        p_throat = _throat_for_flux(eos, inlet, target, "the orifice")
        return self._result(eos, inlet, p_throat, p_throat, False)

    def mass_flow(self, *, inlet_state: FluidState, p_out: float, fluid) -> OrificeResult:
        """The flow the orifice passes from ``inlet_state`` to a discharge pressure.

        Parameters are as for :meth:`solve`, with ``p_out`` the discharge
        pressure [Pa] instead of a mass flow. Below the choking pressure the
        flow no longer depends on ``p_out``, and the result says so.
        """
        eos = as_equation_of_state(fluid)
        inlet = _inlet(eos, inlet_state)
        p_out = _positive(p_out, "p_out")
        if not p_out < inlet.p:
            raise ValueError("p_out must be below the inlet pressure")
        p_throat, choked = _throat_at_discharge(eos, inlet, p_out)
        return self._result(eos, inlet, p_throat, p_out, choked)

    def solve_inlet_pressure(self, *, mass_flow: float, p_out: float, T_in: float,
                             fluid) -> OrificeResult:
        """The inlet pressure that drives ``mass_flow`` to a discharge at ``p_out``.

        The inlet is at temperature ``T_in`` [K], and its state is resolved
        from that and each trial pressure. This is the question a feed
        plenum asks: what it must hold to deliver a flow into a chamber at
        known pressure. Parameters are otherwise as for :meth:`solve`.
        """
        eos = as_equation_of_state(fluid)
        mass_flow = _mass_flow(mass_flow)
        p_out = _positive(p_out, "p_out")
        T_in = _positive(T_in, "T_in")

        def flow(p_in):
            inlet = eos.state_pT(p_in, T_in)
            p_throat, _ = _throat_at_discharge(eos, inlet, p_out)
            return self.cd * self.area * _ideal_flux(eos, inlet, p_throat)[0]

        rho = eos.state_pT(p_out, T_in).rho
        dp = (mass_flow / (self.cd * self.area))**2 / (2.0 * rho)
        for _ in range(60):
            if flow(p_out + dp) >= mass_flow:
                break
            dp *= 2.0
        else:
            raise OrificeCapacityError("no inlet pressure found that passes the flow", math.nan)
        p_in = brentq(lambda p: (flow(p) if p > p_out else 0.0) - mass_flow,
                      p_out, p_out + dp, xtol=1.0e-10 * p_out, rtol=1.0e-12)
        inlet = eos.state_pT(p_in, T_in)
        p_throat, choked = _throat_at_discharge(eos, inlet, p_out)
        return self._result(eos, inlet, p_throat, p_out, choked)


def required_area(*, mass_flow, inlet_state, dp, fluid, cd) -> float:
    """Geometric area that passes ``mass_flow`` with a pressure drop ``dp`` [m²].

    The inverse of :meth:`Orifice.solve`, for any hole shape: the orifice
    discharges at ``inlet_state.p - dp``, with discharge coefficient ``cd``.

    Raises
    ------
    ValueError
        When ``dp`` is larger than the drop at which the flow chokes. Such an
        orifice would pass the flow at a smaller drop, so the drop could not
        be the one asked for.
    """
    eos = as_equation_of_state(fluid)
    inlet = _inlet(eos, inlet_state)
    mass_flow = _mass_flow(mass_flow)
    cd = _discharge_coefficient(cd)
    dp = _positive(dp, "dp")
    p_out = inlet.p - dp
    if not p_out > 0.0:
        raise ValueError("dp must be below the inlet pressure")
    p_throat, choked = _throat_at_discharge(eos, inlet, p_out)
    if choked:
        raise ValueError(
            f"a drop of {dp:g} Pa chokes the flow; the largest unchoked drop is "
            f"{inlet.p - p_throat:g} Pa")
    return mass_flow / (cd * _ideal_flux(eos, inlet, p_out)[0])


# --------------------------------------------------------------------------
# Orifice plate
# --------------------------------------------------------------------------

class OrificePlate:
    """A thin round orifice plate in a straight pipe.

    Parameters
    ----------
    pipe_diameter : float
        Internal diameter of the pipe upstream and downstream [m].
    orifice_diameter : float
        Diameter of the bore through the plate [m].
    cd : float, optional
        Discharge coefficient. By default it is the Reader-Harris/Gallagher
        correlation of ISO 5167-2 (through ``fluids``), evaluated at each
        solve's Reynolds number.
    taps : {'flange', 'corner', 'D', 'D/2'}, optional
        Pressure-tap arrangement the correlation, and so
        :attr:`OrificePlateResult.differential_pressure`, refer to.
    thickness : float, optional
        Plate thickness [m]. Geometry only.

    Notes
    -----
    The flow follows ISO 5167,
    :math:`\\dot m = C \\varepsilon A \\sqrt{2 \\rho_1 \\Delta p / (1 - \\beta^4)}`,
    with :meth:`expansibility` for ε. That ε is fitted to thin plates, whose
    jet keeps contracting past the bore, and it is not the isentropic nozzle
    expansion :class:`Orifice` uses: for a gas at a 2.5 % tap drop the
    nozzle expansion would overstate the drop by about 1.6 %.

    The permanent loss is ISO 5167-1's recovery relation,

    .. math::

        \\Delta \\bar\\omega = \\frac{\\sqrt{1 - \\beta^4 (1 - C^2)} - C \\beta^2}
        {\\sqrt{1 - \\beta^4 (1 - C^2)} + C \\beta^2} \\Delta p.
    """

    def __init__(self, pipe_diameter, orifice_diameter, *, cd=None, taps="flange", thickness=0.0):
        self.pipe_diameter = _positive(pipe_diameter, "pipe_diameter")
        self.orifice_diameter = _positive(orifice_diameter, "orifice_diameter")
        if not self.orifice_diameter < self.pipe_diameter:
            raise ValueError("orifice_diameter must be smaller than pipe_diameter")
        self.cd = None if cd is None else _discharge_coefficient(cd)
        if taps not in _TAPS:
            raise ValueError(f"taps must be one of {_TAPS}")
        self.taps = taps
        thickness = float(thickness)
        if not (math.isfinite(thickness) and thickness >= 0.0):
            raise ValueError("thickness must be finite and non-negative")
        self.thickness = thickness

    def __repr__(self):
        cd = "Reader-Harris/Gallagher" if self.cd is None else f"{self.cd:g}"
        return (f"OrificePlate(pipe_diameter={self.pipe_diameter:g}, "
                f"orifice_diameter={self.orifice_diameter:g}, cd={cd}, taps={self.taps!r})")

    @property
    def area(self) -> float:
        """Bore area [m²]."""
        return math.pi * self.orifice_diameter**2 / 4.0

    @property
    def beta(self) -> float:
        """Diameter ratio, bore over pipe."""
        return self.orifice_diameter / self.pipe_diameter

    def discharge_coefficient(self, *, reynolds: float, density: float = 1.0) -> float:
        """The plate's discharge coefficient at a pipe Reynolds number.

        A supplied ``cd`` is returned as it is. Otherwise this is
        Reader-Harris/Gallagher, which depends only on the Reynolds number
        and geometry; ``density`` does not change the result.
        """
        if self.cd is not None:
            return self.cd
        reynolds = _positive(reynolds, "reynolds")
        self._check_range(reynolds)
        # fluids takes the Reynolds number through rho, mu and m; any
        # combination that gives the same Re_D gives the same C.
        mu = 1.0
        m = reynolds * math.pi * self.pipe_diameter * mu / 4.0
        return float(fluids.flow_meter.C_Reader_Harris_Gallagher(
            D=self.pipe_diameter, Do=self.orifice_diameter, rho=density, mu=mu, m=m,
            taps=self.taps))

    def _check_range(self, reynolds):
        beta = self.beta
        D_mm = self.pipe_diameter * 1e3
        problems = []
        if not 0.1 <= beta <= 0.75:
            problems.append(f"beta={beta:.3g} outside 0.1-0.75")
        if not 50.0 <= D_mm <= 1000.0:
            problems.append(f"pipe diameter {D_mm:.3g} mm outside 50-1000 mm")
        if self.orifice_diameter < 12.5e-3:
            problems.append(f"bore {self.orifice_diameter*1e3:.3g} mm below 12.5 mm")
        if self.taps == "flange":
            re_min = max(5000.0, 170.0 * beta**2 * D_mm)
        else:
            re_min = 5000.0 if beta <= 0.56 else 16000.0 * beta**2
        if reynolds < re_min:
            problems.append(f"Re_D={reynolds:.3g} below {re_min:.3g}")
        if problems:
            warnings.warn(
                "Reader-Harris/Gallagher is extrapolated: " + "; ".join(problems),
                OrificeCorrelationWarning,
                stacklevel=4,
            )

    def expansibility(self, eos: EquationOfState, inlet: FluidState, p_tap: float) -> float:
        """ISO 5167-2's expansibility factor, ε, between the taps.

        The isentropic exponent is the secant ``ln(p1/p2)/ln(ρ1/ρ2)`` along
        the isentrope between the tap pressures, so it applies to real
        fluids, and a constant density gives ε = 1.
        """
        tap = eos.isentropic(inlet, p_tap)
        density_ratio = math.log(inlet.rho / tap.rho)
        if not density_ratio > 0.0:
            return 1.0
        kappa = math.log(inlet.p / p_tap) / density_ratio
        b4 = self.beta**4
        return 1.0 - (0.351 + 0.256 * b4 + 0.93 * b4**2) * (1.0 - (p_tap / inlet.p)**(1.0 / kappa))

    def permanent_loss_fraction(self, cd: float) -> float:
        """Permanent loss over tap-to-tap drop, from ISO 5167-1."""
        b2 = self.beta**2
        root = math.sqrt(1.0 - b2**2 * (1.0 - cd**2))
        return (root - cd * b2) / (root + cd * b2)

    def solve(self, *, inlet_state: FluidState, mass_flow: float, fluid,
              viscosity: float | None = None) -> OrificePlateResult:
        """The tap-to-tap drop and permanent loss at ``mass_flow``.

        Parameters
        ----------
        inlet_state : FluidState
            State at the upstream tap. See :meth:`Orifice.solve`.
        mass_flow : float
            Mass flow [kg/s].
        fluid : str or Fluid or EquationOfState
            The fluid; see :meth:`Orifice.solve`.
        viscosity : float, optional
            Dynamic viscosity at the inlet [Pa s], for the Reynolds number.
            Looked up from CoolProp when ``fluid`` is a CoolProp fluid; it
            must be given otherwise. Not used when ``cd`` was supplied, except
            to report the Reynolds number.
        """
        eos = as_equation_of_state(fluid)
        inlet = _inlet(eos, inlet_state)
        mass_flow = _mass_flow(mass_flow)
        mu = self._viscosity(eos, inlet, viscosity)
        reynolds = 4.0 * mass_flow / (math.pi * self.pipe_diameter * mu) if mu else math.nan
        cd = self.discharge_coefficient(reynolds=reynolds, density=inlet.rho)
        b4 = self.beta**4
        liquid_dp = (mass_flow / (cd * self.area))**2 * (1.0 - b4) / (2.0 * inlet.rho)
        dp = liquid_dp
        for _ in range(100):
            if not dp < inlet.p:
                raise OrificeCapacityError(
                    "the orifice plate would need a drop above its inlet pressure",
                    math.nan)
            new = liquid_dp / self.expansibility(eos, inlet, inlet.p - dp)**2
            if abs(new - dp) <= 1.0e-12 * inlet.p:
                dp = new
                break
            dp = new
        else:
            raise RuntimeError("orifice-plate expansibility iteration did not converge")
        p_tap = inlet.p - dp
        if p_tap / inlet.p < 0.75:
            warnings.warn(
                f"tap pressure ratio {p_tap / inlet.p:.3g} is below ISO 5167's 0.75",
                OrificeCorrelationWarning,
                stacklevel=2,
            )
        _, throat, velocity = _ideal_flux(eos, inlet, p_tap, self.beta**2)
        loss = self.permanent_loss_fraction(cd) * dp
        return OrificePlateResult(
            inlet_state=inlet,
            throat_state=throat,
            outlet_state=_discharge_state(eos, inlet.p - loss, inlet.h, inlet.T),
            mass_flow=mass_flow,
            discharge_coefficient=cd,
            velocity=velocity,
            choked=False,
            beta=self.beta,
            reynolds=reynolds,
        )

    def _viscosity(self, eos, inlet, viscosity):
        if viscosity is not None:
            return _positive(viscosity, "viscosity")
        if isinstance(eos, CoolPropEOS):
            return float(CP.PropsSI("V", "Hmass", inlet.h, "P", inlet.p, eos.fluid))
        if self.cd is None:
            raise ValueError("viscosity is required for Reader-Harris/Gallagher "
                             "when the fluid is not a CoolProp fluid")
        return None
