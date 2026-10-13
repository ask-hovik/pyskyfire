"""Meanline design of an axial turbine, from its duty to a :class:`TurbineDefinition`.

The design runs at mean radius, which is the same at every station, so blade
speed is too. Each stage is a stator and a rotor, and for each one you choose:

- a flow coefficient ``phi = c_x/u``, which sets the axial velocity through
  that stage;
- either its reaction, the share of the stage's static enthalpy drop taken
  in the rotor, or the swirl it leaves with. Given one, the other follows.

The work each stage does is its share of ``efficiency * Δh_s``, where the
isentropic head and its basis (total or static outlet) are the duty's. Euler's
equation then fixes how much the rotor turns the flow.

The losses matter to the design only through density. The entropy rise that
the efficiency implies is shared among the rows, evenly unless ``loss_split``
says otherwise. Each station's static state comes from its enthalpy and the
entropy accumulated so far, and continuity at that density gives the blade
height. Efficiency is taken as given, however it was arrived at.

Row settings then turn the triangles into blades, each kind of row with its
own section profile. Blade height and aspect ratio give the chord. Pitch-to-chord or a Zweifel coefficient gives the blade
count. Metal angles are the flow angles plus incidence and deviation. Across
the span each row follows a free vortex.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from pyskyfire.common.rotation import COUNTERCLOCKWISE
from pyskyfire.thermo.state import FluidState

from .definition import (BladeRow, FreeVortexTwist, ThicknessDistribution, TurbineDefinition,
                         metal_angles)
from .duty import TurbineDuty

__all__ = ["Stage", "SectionProfile", "RowSettings", "StationFlow", "StageFlow", "AxialTurbine"]


def _positive(value, label):
    if value is None or not (np.isfinite(value) and value > 0):
        raise ValueError(f"{label} must be a positive finite number, got {value!r}")
    return float(value)


@dataclass(frozen=True)
class Stage:
    """The velocity-triangle choices for one stage.

    Give ``reaction`` or ``exit_swirl``, not both.

    Parameters
    ----------
    flow_coefficient : float
        ``c_x/u`` through the stage.
    reaction : float, optional
        Rotor share of the stage's static enthalpy drop.
    exit_swirl : float, optional
        Absolute swirl leaving the rotor, positive in the rotation direction
        [m s⁻¹]. Zero leaves the stage axially.
    """

    flow_coefficient: float
    reaction: float | None = None
    exit_swirl: float | None = None

    def __post_init__(self):
        _positive(self.flow_coefficient, "flow_coefficient")
        if (self.reaction is None) == (self.exit_swirl is None):
            raise ValueError("A stage takes a reaction or an exit swirl, one of the two")
        for name in ("reaction", "exit_swirl"):
            value = getattr(self, name)
            if value is not None and not np.isfinite(value):
                raise ValueError(f"{name} must be finite")


@dataclass(frozen=True)
class SectionProfile:
    """The section shape of one kind of row, in fractions of its chord.

    Parameters
    ----------
    max_thickness : float
        Section thickness over chord.
    leading_edge_diameter : float
        Leading-edge thickness over chord.
    trailing_edge_thickness : float
        Trailing-edge thickness over chord.
    trailing_edge_wedge_deg : (float, float)
        Suction- and pressure-side trailing-edge wedge angles.
    distribution : ThicknessDistribution
        How the thickness is laid along the camber line.
    stagger_fraction : float
        Where the chord line sits between the two metal angles: 0 along the
        inlet metal angle, 1 along the outlet. The default, 0.5, is their
        mean, which is also where a symmetric impulse blade has it. A row that
        turns the flow a long way from axial needs a larger fraction for the
        built section's throat to leave the flow near its outlet metal angle;
        the measured gauging angle shows how close it comes.
    stagger_deg : float, optional
        Stagger at mean radius, used instead of ``stagger_fraction``.

    :meth:`impulse` is a starting point for impulse rotor blades: thick at
    mid-chord and thin at both edges, so the leading edge does not narrow the
    steep inlet the flow comes in through.
    """

    max_thickness: float = 0.2
    leading_edge_diameter: float = 0.06
    trailing_edge_thickness: float = 0.03
    trailing_edge_wedge_deg: tuple = (5.0, 2.0)
    distribution: ThicknessDistribution = ThicknessDistribution()
    stagger_fraction: float = 0.5
    stagger_deg: float | None = None

    def __post_init__(self):
        for name in ("max_thickness", "leading_edge_diameter", "trailing_edge_thickness"):
            _positive(getattr(self, name), name)
        if not 0 <= self.stagger_fraction <= 1:
            raise ValueError("stagger_fraction must lie in [0, 1]")
        if not isinstance(self.distribution, ThicknessDistribution):
            raise TypeError("distribution must be a ThicknessDistribution")

    @classmethod
    def impulse(cls, **changes):
        """An impulse blade: symmetric thickness and equally thin edges."""
        base = dict(max_thickness=0.15, leading_edge_diameter=0.03, trailing_edge_thickness=0.03,
                    trailing_edge_wedge_deg=(4.0, 4.0), distribution=ThicknessDistribution.symmetric())
        base.update(changes)
        return cls(**base)

    def stagger(self, inlet_metal, outlet_metal):
        """Stagger at mean radius for a section between these metal angles [rad]."""
        if self.stagger_deg is not None:
            return math.radians(self.stagger_deg)
        return inlet_metal + self.stagger_fraction*(outlet_metal - inlet_metal)


@dataclass(frozen=True)
class RowSettings:
    """How the triangles become blades.

    Parameters
    ----------
    aspect_ratio : float
        Mean blade height over true chord.
    pitch_chord : float, optional
        Blade spacing over true chord, at mean radius.
    zweifel : float, optional
        Zweifel's tangential loading coefficient, used instead of
        ``pitch_chord`` to set the spacing from the turning. With neither
        given, 0.8.
    incidence_deg, deviation_deg : float
        See :func:`~pyskyfire.turbine.definition.metal_angles`.
    stator, rotor : SectionProfile
        Section shape of every stator row and every rotor row.
    axial_gap : float
        Gap between rows over the upstream row's axial chord.
    inlet_length, outlet_length : float
        Straight annulus ahead of the first row and behind the last, over that
        row's axial chord.
    tip_clearance : float
        Rotor tip gap [m]. Zero is a shrouded rotor.
    hub_clearance : float
        Stator hub gap [m], for cantilevered stators over a rotating drum.
    """

    aspect_ratio: float = 1.5
    pitch_chord: float | None = None
    zweifel: float | None = None
    incidence_deg: float = 0.0
    deviation_deg: float = 0.0
    stator: SectionProfile = SectionProfile()
    rotor: SectionProfile = SectionProfile()
    axial_gap: float = 0.3
    inlet_length: float = 0.5
    outlet_length: float = 0.5
    tip_clearance: float = 0.0
    hub_clearance: float = 0.0

    def __post_init__(self):
        if self.pitch_chord is not None and self.zweifel is not None:
            raise ValueError("RowSettings takes pitch_chord or zweifel, not both")
        if self.pitch_chord is None and self.zweifel is None:
            object.__setattr__(self, "zweifel", 0.8)
        _positive(self.aspect_ratio, "aspect_ratio")
        for name in ("pitch_chord", "zweifel"):
            if getattr(self, name) is not None:
                _positive(getattr(self, name), name)
        for name in ("stator", "rotor"):
            if not isinstance(getattr(self, name), SectionProfile):
                raise TypeError(f"{name} must be a SectionProfile")
        for name in ("axial_gap", "inlet_length", "outlet_length", "tip_clearance", "hub_clearance"):
            if not (np.isfinite(getattr(self, name)) and getattr(self, name) >= 0):
                raise ValueError(f"{name} must be nonnegative")


@dataclass(frozen=True)
class StationFlow:
    """The flow at mean radius between two rows (or at either end).

    Angles are radians from axial, positive in the rotation direction.
    Relative quantities use the blade speed ``u``, which is the same at every
    station; they only mean something next to a rotor.
    """

    name: str
    c_x: float
    c_u: float
    u: float
    h0: float
    static: FluidState
    p0: float
    speed_of_sound: float

    @property
    def c(self) -> float:
        return math.hypot(self.c_x, self.c_u)

    @property
    def alpha(self) -> float:
        """Absolute flow angle."""
        return math.atan2(self.c_u, self.c_x)

    @property
    def w_u(self) -> float:
        return self.c_u - self.u

    @property
    def w(self) -> float:
        return math.hypot(self.c_x, self.w_u)

    @property
    def beta(self) -> float:
        """Relative flow angle."""
        return math.atan2(self.w_u, self.c_x)

    @property
    def mach(self) -> float:
        return self.c/self.speed_of_sound

    @property
    def mach_relative(self) -> float:
        return self.w/self.speed_of_sound


@dataclass(frozen=True)
class StageFlow:
    """One stage's triangles: stator inlet (1), rotor inlet (2), rotor exit (3)."""

    inlet: StationFlow
    middle: StationFlow
    exit: StationFlow

    @property
    def work(self) -> float:
        """Euler work ``u*(c_u2 - c_u3)`` [J kg⁻¹]."""
        return self.middle.u*(self.middle.c_u - self.exit.c_u)

    @property
    def loading(self) -> float:
        """``Δh0/u²``."""
        return self.work/self.middle.u**2

    @property
    def flow_coefficient(self) -> float:
        return self.exit.c_x/self.exit.u

    @property
    def reaction(self) -> float:
        """Rotor share of the stage's static enthalpy drop."""
        h1, h2, h3 = (s.static.h for s in (self.inlet, self.middle, self.exit))
        return (h2 - h3)/(h1 - h3)


class AxialTurbine:
    """Design an axial turbine at mean radius. See the module docstring.

    Parameters
    ----------
    duty : TurbineDuty
        What the turbine is asked to do.
    efficiency : float
        Actual work over ``duty.isentropic_head``, on the duty's own basis.
    stages : sequence of Stage
        Stator-rotor pairs, in flow order.
    mean_diameter : float, optional
        Blade mean diameter [m].
    blade_jet_ratio : float, optional
        ``u/c0``, used instead of ``mean_diameter``.
    work_split : sequence of float, optional
        Each stage's share of the work. Default: equal.
    inlet_swirl : float
        Absolute swirl entering the first stator [m s⁻¹].
    loss_split : sequence of float, optional
        Each row's share of the entropy rise. Default: equal.
    rows : RowSettings
        How the triangles become blades.
    rotation : str
        Which way the rotor turns about +x.

    Attributes
    ----------
    stations : tuple of StationFlow
        The flow at every row edge, inlet to exit.
    stages : tuple of StageFlow
    definition : TurbineDefinition
    """

    def __init__(self, duty, efficiency, stages, *, mean_diameter=None, blade_jet_ratio=None,
                 work_split=None, inlet_swirl=0.0, loss_split=None, rows=None,
                 rotation=COUNTERCLOCKWISE):
        if not isinstance(duty, TurbineDuty):
            raise TypeError("duty must be a TurbineDuty")
        if not 0 < efficiency <= 1:
            raise ValueError("efficiency must lie in (0, 1]")
        stages = tuple(stages)
        if not stages or not all(isinstance(s, Stage) for s in stages):
            raise ValueError("stages must be a non-empty sequence of Stage")
        if (mean_diameter is None) == (blade_jet_ratio is None):
            raise ValueError("Give mean_diameter or blade_jet_ratio, one of the two")
        n_rows = 2*len(stages)
        work_split = self._shares(work_split, len(stages), "work_split")
        loss_split = self._shares(loss_split, n_rows, "loss_split")
        if not np.isfinite(inlet_swirl):
            raise ValueError("inlet_swirl must be finite")
        self.duty = duty
        self.efficiency = float(efficiency)
        self.stage_settings = stages
        self.work_split = work_split
        self.loss_split = loss_split
        self.inlet_swirl = float(inlet_swirl)
        self.settings = rows if rows is not None else RowSettings()
        if not isinstance(self.settings, RowSettings):
            raise TypeError("rows must be RowSettings")
        self.rotation = rotation

        if mean_diameter is None:
            mean_diameter = 2*_positive(blade_jet_ratio, "blade_jet_ratio") \
                * duty.spouting_velocity/duty.omega
        self.mean_radius = _positive(mean_diameter, "mean_diameter")/2
        self.blade_speed = duty.omega*self.mean_radius
        self.work = self.efficiency*duty.isentropic_head

        self._solve_thermodynamics()
        self.definition = self._build_definition()

    @staticmethod
    def _shares(values, n, label):
        if values is None:
            return (1.0/n,)*n
        values = tuple(float(v) for v in values)
        if len(values) != n or any(not (v > 0) for v in values):
            raise ValueError(f"{label} needs {n} positive shares")
        total = sum(values)
        return tuple(v/total for v in values)

    # -- meanline -----------------------------------------------------------
    def _swirls(self, stage, c_in, c_x, work):
        """Rotor inlet and exit swirl for one stage, from its reaction or exit swirl."""
        u = self.blade_speed
        turning = work/u
        if stage.exit_swirl is not None:
            c_u3 = stage.exit_swirl
            return c_u3 + turning, c_u3
        # R*(h1 - h3) = h2 - h3 with c_x2 = c_x3 and Euler's c_u2 = c_u3 + work/u,
        # which is a quadratic in c_u3. The root wanted is the one nearest the
        # repeating-stage answer, where c1 = c3.
        R = stage.reaction
        a = R/2
        b = turning
        c = R*(work - 0.5*(c_in**2 - c_x**2)) - work + 0.5*turning**2
        guess = u*(1 - R) - 0.5*turning
        if abs(a) < 1e-12:
            c_u3 = -c/b
        else:
            disc = b*b - 4*a*c
            if disc < 0:
                raise ValueError(f"No stage has reaction {R:g} with this work and flow "
                                 "coefficient: the rotor would have to add enthalpy")
            roots = ((-b + math.sqrt(disc))/(2*a), (-b - math.sqrt(disc))/(2*a))
            c_u3 = min(roots, key=lambda r: abs(r - guess))
        return c_u3 + turning, c_u3

    def _solve_thermodynamics(self):
        duty, eos, u = self.duty, self.duty.eos, self.blade_speed
        s_in = eos.entropy(duty.inlet)
        h01 = duty.inlet.h

        # Swirl and axial velocity at every station, then the enthalpy each carries.
        velocities, h0s = [], []
        c_x_in = self.stage_settings[0].flow_coefficient*u
        c_u = self.inlet_swirl
        velocities.append((c_x_in, c_u))
        h0s.append(h01)
        h0 = h01
        for stage, share in zip(self.stage_settings, self.work_split):
            c_x = stage.flow_coefficient*u
            work = share*self.work
            c_in = math.hypot(*velocities[-1])
            c_u2, c_u3 = self._swirls(stage, c_in, c_x, work)
            velocities += [(c_x, c_u2), (c_x, c_u3)]
            h0s += [h0, h0 - work]
            h0 -= work

        # The entropy the efficiency implies, from the outlet it lands at.
        c_exit = math.hypot(*velocities[-1])
        if duty.outlet_basis == "total":
            s_out = eos.entropy(eos.state_ph(duty.outlet.p, h0))
        else:
            s_out = eos.entropy(eos.state_ph(duty.outlet.p, h0 - 0.5*c_exit**2))
        rise = s_out - s_in
        if rise < -1e-9*abs(s_in or 1):
            if duty.outlet_basis == "static":
                raise ValueError(
                    f"A total-to-static efficiency of {self.efficiency:g} leaves "
                    f"{(1 - self.efficiency)*duty.isentropic_head/1e3:.3g} kJ/kg of loss, less than the "
                    f"{0.5*c_exit**2/1e3:.3g} kJ/kg the exit velocity carries away")
            raise ValueError("The efficiency and duty give an entropy fall")
        entropy = s_in + rise*np.r_[0.0, np.cumsum(self.loss_split)]

        names = ["inlet"]
        for k in range(1, len(self.stage_settings) + 1):
            names += [f"S{k} exit", f"R{k} exit"]
        stations = []
        for name, (c_x, c_u), h0_i, s_i in zip(names, velocities, h0s, entropy):
            static = eos.state_hs(h0_i - 0.5*(c_x**2 + c_u**2), s_i)
            total = eos.state_hs(h0_i, s_i)
            stations.append(StationFlow(name, c_x, c_u, u, h0_i, static, total.p,
                                        eos.speed_of_sound(static)))
        self.stations = tuple(stations)
        self.stages = tuple(StageFlow(*self.stations[2*k:2*k + 3])
                            for k in range(len(self.stage_settings)))

    @property
    def exit_swirl(self) -> float:
        """Absolute swirl leaving the last rotor [m s⁻¹]."""
        return self.stations[-1].c_u

    @property
    def exit_flow_angle(self) -> float:
        return self.stations[-1].alpha

    def station_height(self, station) -> float:
        """Annulus height continuity asks for at a station [m]."""
        area = self.duty.mdot/(station.static.rho*station.c_x)
        return area/(2*math.pi*self.mean_radius)

    # -- hardware ---------------------------------------------------------------
    def _blade_count(self, flow_in, flow_out, chord, axial_chord):
        s = self.settings
        if s.pitch_chord is not None:
            pitch = s.pitch_chord*chord
        else:
            # Zweifel: Z = 2 (s/b) cos²a2 |tan a2 - tan a1|, for axial chord b.
            turning = abs(math.tan(flow_out) - math.tan(flow_in))
            pitch = s.zweifel*axial_chord/(2*math.cos(flow_out)**2*max(turning, 1e-6))
        return max(2, round(2*math.pi*self.mean_radius/pitch))

    def _build_definition(self):
        s, r_m = self.settings, self.mean_radius
        heights = [self.station_height(st) for st in self.stations]
        incidence, deviation = math.radians(s.incidence_deg), math.radians(s.deviation_deg)
        rows, x = [], 0.0
        previous_axial_chord = None
        for k, stage in enumerate(self.stages):
            for j, (kind, a, b, h_in, h_out) in enumerate((
                    ("stator", stage.inlet, stage.middle, heights[2*k], heights[2*k + 1]),
                    ("rotor", stage.middle, stage.exit, heights[2*k + 1], heights[2*k + 2]))):
                rotor = kind == "rotor"
                flow_in, flow_out = (a.beta, b.beta) if rotor else (a.alpha, b.alpha)
                inlet_metal, outlet_metal = metal_angles(flow_in, flow_out, incidence, deviation)
                profile = s.rotor if rotor else s.stator
                stagger = profile.stagger(inlet_metal, outlet_metal)
                chord = 0.5*(h_in + h_out)/s.aspect_ratio
                axial_chord = chord*math.cos(stagger)
                if previous_axial_chord is not None:
                    x += s.axial_gap*previous_axial_chord
                twist = FreeVortexTwist(r_m, a.c_x, b.c_x, a.c_u, b.c_u,
                                        self.blade_speed if rotor else 0.0, incidence, deviation)
                rows.append(BladeRow(
                    name=f"{'R' if rotor else 'S'}{k + 1}", kind=kind,
                    blade_count=self._blade_count(flow_in, flow_out, chord, axial_chord),
                    x_le=x, axial_chord=axial_chord,
                    r_hub_in=r_m - h_in/2, r_tip_in=r_m + h_in/2,
                    r_hub_out=r_m - h_out/2, r_tip_out=r_m + h_out/2,
                    inlet_metal_angle=inlet_metal, outlet_metal_angle=outlet_metal,
                    stagger=stagger,
                    leading_edge_diameter=profile.leading_edge_diameter*chord,
                    trailing_edge_thickness=profile.trailing_edge_thickness*chord,
                    max_thickness=profile.max_thickness*chord,
                    trailing_edge_wedge_angles=tuple(math.radians(w)
                                                     for w in profile.trailing_edge_wedge_deg),
                    thickness_distribution=profile.distribution,
                    tip_clearance=s.tip_clearance if rotor else 0.0,
                    hub_clearance=0.0 if rotor else s.hub_clearance,
                    twist=twist))
                x += axial_chord
                previous_axial_chord = axial_chord
        return TurbineDefinition(
            tuple(rows), rotation=self.rotation,
            inlet_length=s.inlet_length*rows[0].axial_chord,
            outlet_length=s.outlet_length*rows[-1].axial_chord)
