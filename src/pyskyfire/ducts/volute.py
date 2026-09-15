"""Volute geometry generation and one-dimensional loss analysis.

A :class:`Volute` resolves its own geometry on construction: it sweeps the
polar angle from the tongue, and at each station solves the area law for the
cross-section size that satisfies it. :meth:`Volute.solve` then runs Aungier's
meanline loss analysis over the resolved geometry.

Stations follow Aungier's numbering:

``5``
    Volute inlet, where the swirling flow is collected.
``6``
    Full-collection plane, one revolution from the tongue, where the whole mass
    flow has been gathered and turns into the exit cone.
``7``
    Exit cone discharge, at the round pipe the volute delivers to.

The analysis is incompressible. Aungier iterates on the full-collection-plane
velocity because gas density responds to the losses; for the liquids pyskyfire
pumps that feedback is weak, so density is re-evaluated at the mean pressure a
few times rather than solved as a thermodynamic loop. Aungier himself notes
that density variation through a volute is nearly negligible even for gas.

References
----------
Aungier, R. H., *Centrifugal Compressors: A Strategy for Aerodynamic Design and
Analysis*, ASME Press, 2000. Sections 5.7 (Eqs. 5-105 to 5-110) and 11.2
(Eqs. 11-1 to 11-5).

Weber, C. R. and Koronowski, M. E., "Meanline Performance Prediction of Volutes
in Centrifugal Compressors," ASME 86-GT-216, 1986.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import warnings

import fluids
import numpy as np
from scipy.optimize import brentq

from pyskyfire.thermo.state import FluidState

from .area_law import AreaLaw, ConstantAngularMomentum, VoluteInlet
from .exit_cone import ExitCone, resolve_exit_cone
from pyskyfire.common.friction import darcy_friction_factor
from .pipe import _property, _state
from .volute_section import VoluteSection

__all__ = [
    "Volute",
    "VoluteGeometry",
    "VoluteResult",
    "VoluteConsistencyWarning",
]

#: Relative mass-flow mismatch above which the inlet velocities are reported as
#: inconsistent with the mass flow being solved.
_CONSISTENCY_TOLERANCE = 0.01


class VoluteConsistencyWarning(RuntimeWarning):
    """Volute inputs are mutually inconsistent or outside a model's range."""


@dataclass(frozen=True)
class VoluteGeometry:
    """Resolved passage geometry around the volute.

    Attributes
    ----------
    theta : np.ndarray
        Polar angle from the tongue [rad, shape (N,)], ending at ``2*pi``.
    size : np.ndarray
        Section size parameter at each station [m, shape (N,)].
    aspect : np.ndarray
        Section aspect ratio at each station [shape (N,)]. Equal to the
        section family's design aspect ratio wherever ``size`` alone could
        satisfy the area law, and flattened above it near the tongue where it
        could not.
    area : np.ndarray
        Passage area at each station [m^2, shape (N,)].
    centroid_radius : np.ndarray
        Radius of the section area centroid [m, shape (N,)].
    clipped : np.ndarray
        True where the area law asked for less area than the family can
        produce even at its smallest size and flattest aspect, so that
        smallest section was substituted [bool, shape (N,)].
    """

    theta: np.ndarray
    size: np.ndarray
    aspect: np.ndarray
    area: np.ndarray
    centroid_radius: np.ndarray
    clipped: np.ndarray

    @property
    def tongue_clip_angle(self) -> float:
        """Azimuth up to which the area law could not be satisfied [rad].

        Zero when the law is satisfied everywhere. A finite value means the
        passage near the tongue is larger than conservation of angular momentum
        wants, because the cross-section family cannot shrink any further --
        at its smallest size *and* its flattest allowed aspect ratio.
        """
        if not self.clipped.any():
            return 0.0
        return float(self.theta[self.clipped].max())

    @property
    def flattened(self):
        """True where the aspect ratio was pushed past its design value.

        These are the stations near the tongue that a fixed-aspect section
        would have had to clip, and that flattening kept on the area law
        instead. Flattening only ever raises the aspect ratio, so the smallest
        one around the volute is the design value it was raised from.
        """
        return self.aspect > self.aspect.min() + 1.0e-12


@dataclass(frozen=True)
class VoluteResult:
    """Loss breakdown and station conditions for one volute operating point.

    All loss coefficients are referred to the inlet velocity head, so a total
    pressure loss follows as ``omega * rho * C5**2 / 2``.

    Attributes
    ----------
    geometry : VoluteGeometry
        The geometry the analysis ran on.
    SP : float
        Realised sizing parameter, ``r5*C_u5 / (r6*C6)`` (Aungier Eq. 5-106).
        Unity is the minimum-loss value.
    omega_meridional, omega_tangential, omega_friction, omega_exit_cone : float
        Loss coefficients from Aungier Eqs. 5-105, 5-107 (with 5-108 applied
        for a collector), 5-109 and 5-110.
    omega_exit_friction, omega_exit_bend : float
        Wall friction along the exit-cone duct, and the loss of turning it.
        Both are zero without an attached
        :class:`~pyskyfire.ducts.ExitCone`, since there is then no duct to
        charge them against. See :meth:`loss_breakdown` on how they sit
        alongside Aungier's own exit-cone term.
    omega_total : float
        Their sum.
    pressure_loss : float
        Total pressure loss across the volute and its exit cone [Pa].
    C_5, C_6, C_7 : float
        Velocity magnitudes at the three stations [m/s].
    A_6, r_6, d_H : float
        Full-collection-plane area [m^2], centroid radius [m], and Aungier's
        equivalent hydraulic diameter ``sqrt(4*A6/pi)`` [m].
    path_length : float
        Mean flow path length ``pi*(r5 + r6)/2`` used by the friction term [m].
    reynolds, friction_factor : float
        Reynolds number on ``d_H`` and the Fanning friction factor ``c_f``.
    mass_flow, density : float
        The operating point the analysis ran at [kg/s] and the density it
        resolved [kg/m^3].
    states : tuple of FluidState
        Inlet and discharge states.
    """

    geometry: VoluteGeometry
    SP: float
    omega_meridional: float
    omega_tangential: float
    omega_friction: float
    omega_exit_cone: float
    omega_exit_friction: float
    omega_exit_bend: float
    omega_total: float
    pressure_loss: float
    C_5: float
    C_6: float
    C_7: float
    A_6: float
    r_6: float
    d_H: float
    path_length: float
    reynolds: float
    friction_factor: float
    mass_flow: float
    density: float
    states: tuple

    @property
    def total_pressure_loss(self) -> float:
        """Alias matching :class:`~pyskyfire.ducts.pipe.PipeResult`."""
        return self.pressure_loss

    def collected_velocity(self):
        """Through-flow velocity at each azimuthal station [m/s].

        Only the fraction ``theta / (2*pi)`` of the mass flow has been gathered
        by azimuth ``theta``, so this is what the passage area at that station
        actually has to carry. It ends at ``C_6``.

        Under conservation of angular momentum it is the *product* with the
        centroid radius that stays constant, not the velocity: the flow slows
        as the passage carries it outward. A constant-velocity law would be the
        one holding this uniform.
        """
        theta = self.geometry.theta
        return (theta / (2.0 * math.pi)) * self.mass_flow / (
            self.density * self.geometry.area
        )

    def collection_loss(self) -> float:
        """Total pressure lost between the inlet and the collection plane [Pa].

        Everything except the exit cone's three terms, which are spent in the
        duct downstream of it.
        """
        collected = (self.omega_meridional + self.omega_tangential
                     + self.omega_friction)
        return collected * self.density * self.C_5**2 / 2.0

    def stagnation_pressure(self):
        """Total pressure at each azimuthal station [Pa].

        Aungier's analysis is lumped: it charges one coefficient per mechanism
        between the inlet and the collection plane, and says nothing about
        where round the volute each is spent. Spreading them linearly with
        azimuth is an assumption made so the passage can be coloured by
        pressure, not a result of the model. It is a defensible one -- the dump
        terms are paid as flow is collected, and collected flow grows linearly
        with azimuth, while the friction term accumulates along a path length
        that also grows with azimuth -- but the endpoints are what the model
        actually pins down.
        """
        fraction = self.geometry.theta / (2.0 * math.pi)
        return self.states[0].p - fraction * self.collection_loss()

    def static_pressure(self):
        """Static pressure at each azimuthal station [Pa].

        The local velocity head taken off the stagnation pressure above, so it
        inherits the same caveat. It starts at the inlet stagnation pressure,
        where nothing has been collected and the passage is still.
        """
        velocity = self.collected_velocity()
        return self.stagnation_pressure() - self.density * velocity**2 / 2.0

    def exit_cone_pressures(self, stations=9):
        """Stagnation and static pressure along the exit cone [Pa].

        Unlike the volute passage, the duct's losses really are spent between
        two stations the model pins down, so the stagnation pressure ramps
        linearly between them rather than apportioning a lumped total. The
        velocity is eased from ``C_6`` to ``C_7`` on the same smoothstep the
        swept duct's section morphs along, so the static pressure follows the
        rendered shape; only its endpoints carry model weight.

        Returns ``(stagnation, static)``, each with ``stations`` values, the
        first at the collection plane and the last at the flange.
        """
        stations = int(stations)
        if stations < 2:
            raise ValueError("stations must be at least two")
        fraction = np.linspace(0.0, 1.0, stations)
        blend = fraction**2 * (3.0 - 2.0*fraction)
        start = self.states[0].p - self.collection_loss()
        duct_loss = (self.omega_exit_cone + self.omega_exit_friction
                     + self.omega_exit_bend) * self.density * self.C_5**2 / 2.0
        stagnation = start - fraction * duct_loss
        velocity = (1.0 - blend) * self.C_6 + blend * self.C_7
        return stagnation, stagnation - self.density * velocity**2 / 2.0

    def loss_breakdown(self) -> dict:
        """Loss coefficients by mechanism, for reporting.

        ``exit_cone`` is Aungier's Eq. 5-110, which charges the velocity change
        between the collection plane and the flange and knows nothing of the
        duct's length or how sharply it turns. ``exit_friction`` and
        ``exit_bend`` add exactly those two, from the resolved duct geometry.

        The split is not perfectly clean. Aungier notes Eq. 5-110 is
        deliberately pessimistic -- "substantially greater" than Weber and
        Koronowski's -- precisely because the flow entering the exit cone is
        distorted and swirling, so some allowance for the turn is already
        buried in it. Summing all three is therefore conservative rather than
        additive-from-first-principles. They are reported separately so that
        can be seen rather than hidden.
        """
        return {
            "meridional": self.omega_meridional,
            "tangential": self.omega_tangential,
            "friction": self.omega_friction,
            "exit_cone": self.omega_exit_cone,
            "exit_friction": self.omega_exit_friction,
            "exit_bend": self.omega_exit_bend,
        }


class Volute:
    """A volute or collector resolved from a section family and an area law.

    Parameters
    ----------
    section : VoluteSection
        Cross-section family. Its fixed inlet width is the ``b5`` of the design
        equation.
    inlet : VoluteInlet
        Inlet radius and velocity components.
    area_law : AreaLaw, optional
        How the passage grows with azimuth. Defaults to
        :class:`~pyskyfire.ducts.area_law.ConstantAngularMomentum` at ``SP=1``,
        the minimum-loss shaped volute.
    exit_cone : ExitCone, optional
        The duct carrying the collected flow to its pipe flange. When given,
        its final diameter is what :meth:`solve` uses for station 7, so the
        flange size is stated once. Without one the volute stops at the
        full-collection plane and ``solve`` needs an explicit ``exit_diameter``.
    tilt : float, optional
        Inclination of the inlet plane away from radial, toward the machine
        axis [rad]. Zero collects the flow radially, with the slot spanning the
        axial direction; a pump whose discharge points a little downward, or a
        manifold following a nozzle's slope, has a nonzero tilt. Constant
        around the circumference.
    n_stations : int, optional
        Number of azimuthal stations, including both the tongue and the
        full-collection plane.

    Notes
    -----
    Geometry is resolved on construction and is independent of the operating
    point, so one volute can be solved at several mass flows.

    Tilt rotates the section plane rigidly, so it changes where the passage
    sits in space but leaves every scalar in the loss analysis unchanged. It
    matters for packaging and for :meth:`sweep`, not for pressure loss.
    """

    def __init__(self, section, inlet, area_law=None, *, exit_cone=None,
                 tilt=0.0, n_stations=181):
        if not isinstance(section, VoluteSection):
            raise TypeError("section must be a VoluteSection")
        if not isinstance(inlet, VoluteInlet):
            raise TypeError("inlet must be a VoluteInlet")
        if area_law is None:
            area_law = ConstantAngularMomentum()
        if not isinstance(area_law, AreaLaw):
            raise TypeError("area_law must be an AreaLaw")
        n_stations = int(n_stations)
        if n_stations < 3:
            raise ValueError("n_stations must be at least three")
        tilt = float(tilt)
        if not np.isfinite(tilt):
            raise ValueError("tilt must be finite")

        self.section = section
        self.inlet = inlet
        self.area_law = area_law
        self.tilt = tilt
        self.geometry = self._resolve(n_stations)
        if exit_cone is not None and not isinstance(exit_cone, ExitCone):
            raise TypeError("exit_cone must be an ExitCone")
        self.exit_cone = exit_cone
        self.exit_geometry = None if exit_cone is None else self._resolve_exit(exit_cone)

    def outlet_center(self):
        """Centre of the volute outlet bore, in machine coordinates [m].

        This is the point an exit cone leaves from, and the axis its flange is
        carried along.
        """
        local = self.section.fan_origin(self.geometry.size[-1:],
                                        self.geometry.aspect[-1:])[0]
        theta = float(self.geometry.theta[-1])
        radial = np.array([math.cos(theta), math.sin(theta), 0.0])
        axial = np.array([0.0, 0.0, 1.0])
        u_hat = math.cos(self.tilt)*radial + math.sin(self.tilt)*axial
        v_hat = -math.sin(self.tilt)*radial + math.cos(self.tilt)*axial
        return self.inlet.radius*radial + local[0]*u_hat + local[1]*v_hat

    def _resolve_exit(self, exit_cone):
        size = self.geometry.size[-1:]
        aspect = self.geometry.aspect[-1:]
        angles = np.linspace(0.0, 2.0*math.pi, 361, endpoint=False)
        return resolve_exit_cone(
            exit_cone,
            start=self.outlet_center(),
            outlet_diameter=2.0*float(size[0]),
            max_section_radius=float(
                self.section.radial_profile(size, angles, aspect).max()),
        )

    def sweep_exit_cone(self, n_points=64, n_stations=17):
        """The exit-cone duct swept into machine coordinates.

        The duct starts as the volute's own outlet section, sampled exactly as
        :meth:`sweep` samples it -- the same call, at the same size, with the
        same ``n_points`` -- so its first ring is bit-identical to the
        volute's own last one rather than merely tracing the same curve at
        different points along it. That identity is what makes the two sweeps
        stitch into one seamless duct wall.

        Each of those vertices then keeps its own angle about the outlet
        centre fixed for the whole duct, while its radius blends by a
        smoothstep toward the flange radius, so the wall leaves the volute and
        meets the flange without a crease. The frame it blends within starts as
        the volute's own tilted section-plane basis, not the untilted
        machine-frame default, so a tilted volute's cutwater edge continues to
        rotate smoothly into the flange instead of snapping to meet it.

        Returns
        -------
        np.ndarray
            Closed boundary rings [m, shape (n_stations, n_points, 3)].
        """
        if self.exit_geometry is None:
            raise ValueError("this volute has no exit cone")
        n_points = int(n_points)
        if n_points < 3:
            raise ValueError("n_points must be at least three")
        cone = self.exit_geometry
        size = self.geometry.size[-1:]
        aspect = self.geometry.aspect[-1:]

        # The outlet ring, in the section's own local (u, v) plane -- the same
        # coordinates fan_origin is expressed in -- gives each vertex's fixed
        # angle and starting radius about the outlet centre.
        ring_uv = self.section.polygon(size, n_points=n_points, aspect=aspect)[0][:-1]
        offset = ring_uv - self.section.fan_origin(size, aspect)[0]
        psi = np.arctan2(offset[:, 1], offset[:, 0])
        start_radius = np.linalg.norm(offset, axis=-1)

        theta = float(self.geometry.theta[-1])
        radial = np.array([math.cos(theta), math.sin(theta), 0.0])
        axial = np.array([0.0, 0.0, 1.0])
        u_hat = math.cos(self.tilt)*radial + math.sin(self.tilt)*axial
        v_hat = -math.sin(self.tilt)*radial + math.cos(self.tilt)*axial

        points, u, v = cone.frames(n_stations, start_frame=(u_hat, v_hat))
        fraction = np.linspace(0.0, 1.0, len(points))
        blend = fraction**2 * (3.0 - 2.0*fraction)
        radius = ((1.0 - blend)[:, None]*start_radius[None, :]
                  + blend[:, None]*(cone.final_diameter/2.0))
        spoke = (np.cos(psi)[None, :, None]*u[:, None, :]
                 + np.sin(psi)[None, :, None]*v[:, None, :])
        return points[:, None, :] + radius[..., None]*spoke

    def sweep(self, n_points=64):
        """Section boundaries swept into machine coordinates.

        The machine axis is ``z`` and the tongue sits at azimuth zero. Each
        station's section is placed in the plane through the axis at its
        azimuth, rotated about the local tangential direction by
        :attr:`tilt`, with the slot mouth centred on the inlet radius.

        Parameters
        ----------
        n_points : int, optional
            Boundary samples per piece of the section outline.

        Returns
        -------
        np.ndarray
            Closed boundary rings [m, shape (N, M, 3)], one per azimuthal
            station.
        """
        ring = self.section.polygon(self.geometry.size, n_points=n_points,
                                    aspect=self.geometry.aspect)
        theta = self.geometry.theta
        radial = np.column_stack((np.cos(theta), np.sin(theta), np.zeros_like(theta)))
        axial = np.zeros_like(radial)
        axial[:, 2] = 1.0
        u_hat = math.cos(self.tilt) * radial + math.sin(self.tilt) * axial
        v_hat = -math.sin(self.tilt) * radial + math.cos(self.tilt) * axial
        origin = self.inlet.radius * radial
        u, v = ring[..., 0], ring[..., 1]
        return (origin[:, None, :]
                + u[..., None] * u_hat[:, None, :]
                + v[..., None] * v_hat[:, None, :])

    # -- geometry ---------------------------------------------------------

    def _residual(self, size, theta, aspect=None):
        """Area-law miss for one trial section."""
        area, centroid, _ = self.section.metrics(size, aspect)
        radius = self.inlet.radius + centroid[:, 0]
        return self.area_law.residual(
            theta, float(area[0]), float(radius[0]),
            self.section.inlet_width, self.inlet,
        )

    def _solve_size(self, theta, minimum):
        """The smallest section satisfying the area law at ``theta``.

        Returns ``(size, aspect, clipped)``, working down three regimes as the
        tongue is approached and the law asks for less and less area.

        Away from the tongue the design aspect ratio can hold the law, and only
        the size is solved for. Closer in, the law wants less area than even
        the smallest section at that aspect: the bore can never be narrower
        than the slot feeding it, so the size has bottomed out. Flattening is
        the way past that -- holding the size at its floor and stretching the
        section across the slot shrinks its area without narrowing the bore --
        so the aspect ratio is solved for instead, up to the family's limit.
        Only when the flattest allowed section is still too big does the law
        get clipped, and it is clipped to that flattest section rather than to
        the design-aspect one, which is as close to the law as the family can
        get.
        """
        design = self.section.aspect_ratio
        ceiling = self.section.max_aspect_ratio

        if self._residual(minimum, theta, design) < 0.0:
            high = max(minimum * 2.0, minimum + self.section.inlet_width)
            for _ in range(60):
                if self._residual(high, theta, design) > 0.0:
                    break
                high *= 2.0
            else:
                raise RuntimeError(
                    f"could not bracket a section size at theta={theta:g} rad; "
                    "the area law does not appear to be satisfiable by this "
                    "section"
                )
            size = brentq(self._residual, minimum, high, args=(theta, design),
                          xtol=1e-14, rtol=1e-12)
            return size, design, False

        if ceiling > design and self._residual(minimum, theta, ceiling) < 0.0:
            aspect = brentq(
                lambda value: self._residual(minimum, theta, value),
                design, ceiling, xtol=1e-12, rtol=1e-12,
            )
            return minimum, aspect, False

        return minimum, ceiling, True

    def _resolve(self, n_stations):
        theta = np.linspace(0.0, 2.0 * math.pi, n_stations)
        minimum = self.section.minimum_size
        sizes = np.empty(n_stations)
        aspects = np.empty(n_stations)
        clipped = np.zeros(n_stations, dtype=bool)
        for index, angle in enumerate(theta):
            sizes[index], aspects[index], clipped[index] = self._solve_size(
                float(angle), minimum)
        area, centroid, _ = self.section.metrics(sizes, aspects)
        return VoluteGeometry(
            theta=theta,
            size=sizes,
            aspect=aspects,
            area=area,
            centroid_radius=self.inlet.radius + centroid[:, 0],
            clipped=clipped,
        )

    # -- performance ------------------------------------------------------

    def _exit_cone_losses(self, rho, mu, C_5, C_6, C_7, roughness, friction_model):
        """Wall friction and bend loss of the exit-cone duct.

        Both are referred to the inlet velocity head, like every other
        coefficient here, and both are zero without an attached exit cone.

        The duct's area varies along its length, from the volute's outlet
        section to the flange, so friction is charged at the mean of the two
        velocity heads over the mean of the two hydraulic diameters rather
        than integrated -- a duct one or two diameters long does not earn a
        marching solution. The bend uses the same ``fluids`` correlation the
        pipe solver uses for a rounded bend, at the duct's own centreline
        radius and turn.
        """
        cone = self.exit_geometry
        if cone is None:
            return 0.0, 0.0

        entry = math.sqrt(4.0 * self.geometry.area[-1] / math.pi)
        diameter = 0.5 * (entry + cone.final_diameter)
        velocity = 0.5 * (C_6 + C_7)
        reynolds = rho * velocity * diameter / mu
        darcy = darcy_friction_factor(reynolds, diameter, roughness,
                                      model=friction_model)
        friction = darcy * (cone.arc_length / diameter) * (velocity / C_5) ** 2

        if cone.is_straight:
            return friction, 0.0
        ratio = cone.bend_radius / diameter
        if ratio < 0.5:
            warnings.warn(
                f"exit-cone bend radius over diameter is {ratio:.2f}, below the "
                "0.5 where the fluids bend correlations stop being fitted; the "
                "bend loss is extrapolated",
                VoluteConsistencyWarning,
                stacklevel=3,
            )
        bend = fluids.fittings.bend_rounded(
            Di=diameter, angle=math.degrees(cone.turn_angle),
            rc=cone.bend_radius, Re=reynolds, roughness=roughness,
            method="Miller",
        ) * (velocity / C_5) ** 2
        return friction, bend

    def solve(
        self,
        *,
        inlet_state: FluidState,
        mass_flow: float,
        transport,
        exit_diameter: float | None = None,
        roughness: float = 0.0,
        friction_model: str = "fluids",
    ) -> VoluteResult:
        """Run Aungier's meanline loss analysis at one operating point.

        Parameters
        ----------
        inlet_state : FluidState
            Conditions at the volute inlet.
        mass_flow : float
            Mass flow through the volute [kg/s]. Checked against the inlet
            velocity components and the inlet slot area, which are independent
            inputs; a mismatch is warned about rather than silently absorbed.
        transport : object
            Property source exposing ``get_rho(T, p)`` and ``get_mu(T, p)``.
        exit_diameter : float, optional
            Diameter of the round pipe the exit cone discharges to [m], which
            sets station 7. Defaults to the attached exit cone's final
            diameter, and is required when there is no exit cone.
        roughness : float, optional
            Absolute wall roughness [m].
        friction_model : str, optional
            Backend passed to
            :func:`~pyskyfire.common.friction.darcy_friction_factor`.
        """
        if not isinstance(inlet_state, FluidState):
            raise TypeError("inlet_state must be a thermo.FluidState")
        mass_flow = float(mass_flow)
        if exit_diameter is None:
            if self.exit_geometry is None:
                raise ValueError(
                    "exit_diameter is required when the volute has no exit cone")
            exit_diameter = self.exit_geometry.final_diameter
        exit_diameter = float(exit_diameter)
        roughness = float(roughness)
        if not np.isfinite(mass_flow) or mass_flow <= 0.0:
            raise ValueError("mass_flow must be finite and positive")
        if not np.isfinite(exit_diameter) or exit_diameter <= 0.0:
            raise ValueError("exit_diameter must be finite and positive")
        if not np.isfinite(roughness) or roughness < 0.0:
            raise ValueError("roughness must be finite and nonnegative")
        T_in, p_in = float(inlet_state.T), float(inlet_state.p)
        if not np.isfinite(T_in) or T_in <= 0.0 or not np.isfinite(p_in) or p_in <= 0.0:
            raise ValueError("inlet state must have finite positive T and p")

        inlet = self.inlet
        geometry = self.geometry
        A_6 = float(geometry.area[-1])
        r_6 = float(geometry.centroid_radius[-1])
        A_7 = math.pi * exit_diameter**2 / 4.0
        d_H = math.sqrt(4.0 * A_6 / math.pi)
        path_length = math.pi * (inlet.radius + r_6) / 2.0
        C_5 = inlet.C

        # The volute is adiabatic here, so only density is re-evaluated as the
        # pressure falls. For a liquid this converges on the first pass.
        p_out = p_in
        for _ in range(8):
            rho = _property(transport, "rho", T_in, 0.5 * (p_in + p_out))
            mu = _property(transport, "mu", T_in, 0.5 * (p_in + p_out))

            C_6 = mass_flow / (rho * A_6)
            C_7 = mass_flow / (rho * A_7)
            SP = inlet.swirl / (r_6 * C_6)

            # Eq. 5-105: the meridional head becomes swirl inside the passage
            # and is never recovered.
            omega_m = (inlet.C_m / C_5) ** 2

            # Eq. 5-107, in two branches about the optimum SP = 1.
            swirl_head = inlet.radius * inlet.C_u**2 / (r_6 * C_5**2)
            if SP >= 1.0:
                omega_u = 0.5 * swirl_head * (1.0 - 1.0 / SP**2)
            else:
                omega_u = swirl_head * (1.0 - 1.0 / SP) ** 2
            if self.area_law.collector:
                # Eq. 5-108: a uniform passage also tries to diffuse the
                # entering flow to rest where collection begins.
                omega_u = (1.0 + omega_u) / 2.0

            # Eq. 5-109, with c_f the Fanning factor.
            reynolds = rho * C_6 * d_H / mu
            c_f = darcy_friction_factor(
                reynolds, d_H, roughness, model=friction_model
            ) / 4.0
            omega_sf = 4.0 * c_f * (C_6 / C_5) ** 2 * path_length / d_H

            # Eq. 5-110.
            omega_ec = ((C_6 - C_7) / C_5) ** 2

            # The duct that carries the flow to the flange: its own wall
            # friction and, where it turns, its bend loss. Aungier's Eq. 5-110
            # above charges only the velocity change and is blind to both.
            omega_ef, omega_eb = self._exit_cone_losses(
                rho, mu, C_5, C_6, C_7, roughness, friction_model)

            omega_total = omega_m + omega_u + omega_sf + omega_ec + omega_ef + omega_eb
            loss = omega_total * rho * C_5**2 / 2.0
            next_pressure = p_in - loss
            if next_pressure <= 0.0:
                raise RuntimeError(
                    f"volute pressure fell to {next_pressure:g} Pa"
                )
            if abs(next_pressure - p_out) < 1.0e-6:
                p_out = next_pressure
                break
            p_out = next_pressure

        slot_flow = rho * inlet.C_m * 2.0 * math.pi * inlet.radius * self.section.inlet_width
        if abs(slot_flow - mass_flow) > _CONSISTENCY_TOLERANCE * mass_flow:
            warnings.warn(
                f"inlet velocities imply {slot_flow:.4g} kg/s through the inlet "
                f"slot but {mass_flow:.4g} kg/s was requested; C_m, the inlet "
                "radius and the slot width are independent inputs and should "
                "agree with the mass flow",
                VoluteConsistencyWarning,
                stacklevel=2,
            )

        return VoluteResult(
            geometry=geometry,
            SP=SP,
            omega_meridional=omega_m,
            omega_tangential=omega_u,
            omega_friction=omega_sf,
            omega_exit_cone=omega_ec,
            omega_exit_friction=omega_ef,
            omega_exit_bend=omega_eb,
            omega_total=omega_total,
            pressure_loss=loss,
            C_5=C_5,
            C_6=C_6,
            C_7=C_7,
            A_6=A_6,
            r_6=r_6,
            d_H=d_H,
            path_length=path_length,
            reynolds=reynolds,
            friction_factor=c_f,
            mass_flow=mass_flow,
            density=rho,
            states=(_state(T_in, p_in, rho), _state(T_in, p_out, rho)),
        )
