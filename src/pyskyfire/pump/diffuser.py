"""Vaned radial diffuser between an impeller and its volute.

A vaneless volute collects the impeller discharge at nearly the speed it left
the wheel, and it carries it round the whole circumference that way, so its
friction is paid at the impeller's own dynamic head. Vanes take the swirl out
sooner. They are fixed, so they do no work: they push on the flow, remove its
angular momentum and turn it outward, so it reaches the volute radius in a
short path, and the channels between them widen, so it slows by their area
ratio rather than only by the radius ratio. The volute behind then collects
slow, low-swirl flow.

The design follows Gülich, chapter 7.9. The inlet is placed and sized by his
rules: the gap to the impeller from the pressure-pulsation limits (Table
10.2), the inlet width from the impeller outlet width, the vane angle from the
approach flow (T7.7.4 to T7.7.6) and the throat from conservation of angular
momentum times a factor (T7.7.10), or from a chosen deceleration (T7.7.11).
Each channel then opens from its throat at the rate that gives it the area
ratio recovering the most pressure for its length (T1.7.13), so its recovery
is Gülich's optimum (T1.7.14), corrected for the channel's curvature
(Eq. 1.46), and not a value read off a chart at some other ratio.

The vanes are thin, of the thickness Gülich gives for a diffuser vane (7.9.2),
with rounded ends. Each is a camber line: a log spiral from the leading edge,
at the angle that makes the throat to the next vane's leading edge the design
throat, then turning outward just fast enough that the channel between two
vanes -- about the pitch times the sine of their angle, less a vane
thickness -- widens at that rate. The channels therefore curve. A straight
channel would leave the throat far steeper than the flow arriving at it, and
the vane between two of them would be a wedge filling most of the ring.

Gülich's quality criteria for the finished diffuser are reported on
:attr:`VanedDiffuser.checks`, never enforced. Only a geometry that cannot be
drawn raises.

References
----------
Gülich, J. F., *Centrifugal Pumps*, 4th ed., Springer, 2020. Tables 1.7, 3.8
(2), 7.7 and 10.2; chapter 7.9.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import atan, cos, degrees, exp, log10, pi, radians, sin, sqrt, tan

import numpy as np

from pyskyfire.common.blade_surface import BladeSurface
from pyskyfire.common.rotation import COUNTERCLOCKWISE, parse_rotation, rotation_sign

__all__ = ["VanedDiffuser", "DiffuserOutlet", "DiffuserChecks", "DiffuserGeometry",
           "DiffuserResult", "minimum_gap_ratio", "vane_count_conflicts"]

G0 = 9.80665


def minimum_gap_ratio(*, nq, rho, H):
    """Smallest d3/d2 Gülich allows for low pressure pulsations (Table 10.2).

    ``H`` is the stage head [m] and ``rho`` its density [kg/m³]. The gap grows
    with the energy density rho*H, which is what drives the rotor-stator
    interaction; at the heads of a rocket pump it is well over ten percent.
    """
    if nq > 40:
        return 1.04+0.001*(nq-40)
    loading = rho*H/(1000*1000)-0.1
    return 1.015+0.08*max(loading, 0.)**0.8


def vane_count_conflicts(vane_count, blade_count):
    """Vane and blade count pairings Gülich says to avoid (7.9.1 item 5, Table 10.2).

    Each entry is ``(nu2, nu3, m)`` with ``m = nu3*vane_count - nu2*blade_count``
    of 0 or +-1 for ``nu`` up to 3: those put a harmonic of the blade passing
    on a harmonic of the vane passing, or one lobe off it, and drive pressure
    pulsations. An empty tuple is a pairing he would accept.
    """
    return tuple((v2, v3, v3*vane_count-v2*blade_count)
                 for v2 in (1, 2, 3) for v3 in (1, 2, 3)
                 if abs(v3*vane_count-v2*blade_count) <= 1)


def _positive(value, name):
    if not np.isscalar(value) or not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return float(value)


def _cf(velocity, length, nu, roughness):
    """Gülich's boundary-layer friction coefficient, Eq. (T3.8.3)."""
    reynolds = max(velocity*length/nu, 1e3)
    return 0.136/(-log10(0.2*roughness/length+12.5/reynolds))**2.15


def _offset(curve, distance):
    """Offset a plane curve to its left (toward the axis for an outward spiral)."""
    tangent = np.gradient(curve, axis=0)
    tangent /= np.linalg.norm(tangent, axis=1)[:, None]
    return curve+np.asarray(distance)[..., None]*np.column_stack((-tangent[:, 1], tangent[:, 0]))


def _rotate(points, angle):
    c, s = cos(angle), sin(angle)
    return np.asarray(points) @ np.array([[c, s], [-s, c]])


def _signed_distance(points, line):
    """Distance from each point to a polyline, positive to its left."""
    a, b = line[:-1], line[1:]
    seg = b-a
    out = np.empty(len(points))
    for i, p in enumerate(points):
        t = np.clip(np.einsum("ij,ij->i", p-a, seg)/np.einsum("ij,ij->i", seg, seg), 0, 1)
        near = a+t[:, None]*seg
        k = int(np.argmin(np.linalg.norm(p-near, axis=1)))
        side = seg[k, 0]*(p-a[k])[1]-seg[k, 1]*(p-a[k])[0]
        out[i] = np.copysign(np.linalg.norm(p-near[k]), side)
    return out


def _resample(curve, count):
    s = np.r_[0., np.cumsum(np.linalg.norm(np.diff(curve, axis=0), axis=1))]
    target = np.linspace(0, s[-1], count)
    return np.column_stack([np.interp(target, s, curve[:, k]) for k in range(2)])


@dataclass(frozen=True)
class DiffuserOutlet:
    """What the diffuser hands the volute, circumferentially mixed.

    ``cm`` is the mean radial velocity over the whole outlet ring and ``cu``
    the mass-averaged swirl of the channel jets, both in m/s. Together with
    ``radius`` and ``width`` they are a volute inlet.
    """
    radius: float
    width: float
    cm: float
    cu: float


@dataclass(frozen=True)
class DiffuserChecks:
    """The finished diffuser against Gülich's criteria (7.9.1 item 10). Reported only.

    ``gap_ratio`` is d3/d2 against ``minimum_gap_ratio`` (Table 10.2).
    ``inlet_width_ratio`` b3/b2 should be 1.05 to 1.3. ``throat_factor`` is
    f_a3, normally 1.1 to 1.3. ``deceleration`` is c3q/c2: 0.7 to 0.85,
    higher allowed below nq 12 to 15, and Table 10.2 asks for not much below
    0.8. ``length_ratio`` L/a3 should lie between 2.5 and 6 and
    ``aspect_ratio`` b3/a3 between 0.8 and 2. ``exit_energy`` is c4²/(2gH),
    0.02 to 0.04 for multistage pumps. ``opening_angle_deg`` is the equivalent
    conical half-angle (Eq. 1.45), which must stay below
    ``allowed_opening_angle_deg`` (Eq. 1.44) or the channel stalls.
    ``vane_count_conflicts`` lists each (nu2, nu3, m) with
    ``|m| = |nu3 z_Le - nu2 z_La|`` of 0 or 1 for nu up to 3, which Gülich says
    to avoid (7.9.1 item 5, Table 10.2). ``incidence_deg`` is the vane's
    leading spiral against the approach flow, within +-3 degrees by his rule;
    it follows from the throat. ``deflection_deg`` is how far the channel
    turns the flow from throat to exit; Gülich applies straight-diffuser data
    up to about 25 degrees. ``drawn_throat_ratio`` and ``drawn_exit_ratio``
    are the throat and exit widths measured on the drawn vanes over the design
    values, which the camber construction meets only approximately.
    """
    gap_ratio: float
    minimum_gap_ratio: float
    inlet_width_ratio: float
    throat_factor: float
    deceleration: float
    length_ratio: float
    aspect_ratio: float
    exit_energy: float
    opening_angle_deg: float
    allowed_opening_angle_deg: float
    vane_count_conflicts: tuple
    incidence_deg: float
    deflection_deg: float
    drawn_throat_ratio: float
    drawn_exit_ratio: float


@dataclass(frozen=True)
class DiffuserGeometry:
    """One vane and the walls it stands between, in the impeller's own frame.

    ``vane`` is a :class:`~pyskyfire.common.BladeSurface` in (x, r, theta)
    about the shaft (+x), in the wheel's rotation sense; ``suction`` is the
    face toward the impeller and ``pressure`` the one away from it. ``hub``
    and ``shroud`` are the (x, r) wall curves from the impeller outlet to the
    diffuser outlet, symmetric about the impeller's outlet mid-plane at x = 0,
    with the hub at +x as on the impeller.
    """
    vane: BladeSurface
    vane_count: int
    hub: np.ndarray
    shroud: np.ndarray
    rotation_sign: float


@dataclass(frozen=True)
class DiffuserResult:
    """Losses of one diffuser operating point, Gülich Table 3.8 (2).

    Coefficients are referred to u2²/2, as Gülich writes them. ``vaneless``
    is the wall friction between impeller and vanes (T3.8.23), ``inlet`` the
    friction in the triangular inlet region (T3.8.18), ``shock`` the loss of
    decelerating c2 onto the throat (the 0.3 term of T3.8.19) and ``channel``
    the diffusing channel, 1 - cp - 1/AR² times the throat head. ``mixing``
    is the channel jets filling the outlet ring: they leave through part of
    the circumference only, the vanes' trailing edges taking the rest, and
    mixing out keeps their swirl but loses their excess radial velocity, a
    Borda-Carnot loss on that component. This is the part of Gülich's overflow
    loss that belongs to the diffuser; the energy left in the mixed flow is not
    charged, because the volute collects it and charges its own loss.
    """
    zeta_vaneless: float
    zeta_inlet: float
    zeta_shock: float
    zeta_channel: float
    zeta_mixing: float
    zeta_total: float
    pressure_loss: float
    c2: float
    c3q: float
    c4: float
    density: float
    mass_flow: float
    states: tuple

    @property
    def total_pressure_loss(self):
        return self.pressure_loss

    def loss_breakdown(self):
        return {"vaneless": self.zeta_vaneless, "inlet": self.zeta_inlet,
                "shock": self.zeta_shock, "channel": self.zeta_channel,
                "mixing": self.zeta_mixing}


class VanedDiffuser:
    """A vaned radial diffuser designed from the impeller it follows.

    Build one with :meth:`from_impeller`. The direct constructor takes the
    impeller outlet as numbers: ``flow`` through the impeller [m³/s], its
    outlet diameter ``d2`` and width ``b2`` [m], tip speed ``u2`` and outlet
    swirl ``c2u`` [m/s], blade count, speed ``n`` [rpm], head ``H`` [m] (the
    head at the impeller outlet, used for the gap and exit-energy rules),
    specific speed ``nq``, density ``rho`` [kg/m³] and ``rotation``.

    Design inputs, with Gülich's ranges:

    ``vane_count``
        Required. Match it against the impeller blade count; see
        :attr:`DiffuserChecks.vane_count_conflicts`.
    ``gap_ratio``
        d3/d2. Defaults to the Table 10.2 minimum.
    ``inlet_width_ratio``
        b3/b2, 1.05 to 1.3. Default 1.15.
    ``outer_diameter_ratio``, ``length_ratio``
        d4/d2, or the channel length over its throat width, L/a3, which
        then sets d4. Give at most one. By default L/a3 is 4, the middle of
        Gülich's 2.5 to 6. His own d4 rule, (1.05 to 1.15) + 0.01 nq, assumes
        a gap of a few percent; the gap a high-head pump needs eats most of
        that, so the channel is sized directly instead.
    ``outlet_width_ratio``
        b4/b3: the walls flare linearly from inlet to outlet radius. The
        flare opens the channel too, and the equivalent angle counts both.
        Default 1, a plane diffuser.
    ``throat_factor``, ``deceleration``
        f_a3 of T7.7.10, 1.1 to 1.3, or instead c3q/c2 of T7.7.11, 0.7 to
        0.85. Give at most one. By default c3q/c2 is 0.8: Table 10.2 asks for
        not much below 0.8 against part-load pulsations, and behind the wide
        gap of a high-head pump f_a3 = 1.2 would decelerate below 0.7.
    ``vane_thickness_ratio``
        Vane thickness over d2. Gülich gives 0.01 to 0.015 at the inlet.
        Default 0.012, held along the vane, with rounded ends.
    ``roughness``
        Wall roughness [m], for the friction terms.
    """

    def __init__(self, *, flow, d2, b2, u2, c2u, impeller_blade_count, n, H, nq, rho,
                 vane_count, gap_ratio=None, inlet_width_ratio=1.15,
                 outer_diameter_ratio=None, length_ratio=None, outlet_width_ratio=1.0,
                 throat_factor=None, deceleration=None, vane_thickness_ratio=0.012,
                 roughness=0.0, rotation=COUNTERCLOCKWISE):
        for name, value in (("flow", flow), ("d2", d2), ("b2", b2), ("u2", u2), ("c2u", c2u),
                            ("n", n), ("H", H), ("nq", nq), ("rho", rho),
                            ("inlet_width_ratio", inlet_width_ratio),
                            ("outlet_width_ratio", outlet_width_ratio),
                            ("vane_thickness_ratio", vane_thickness_ratio)):
            _positive(value, name)
        if isinstance(vane_count, bool) or not isinstance(vane_count, (int, np.integer)) or vane_count < 3:
            raise ValueError("vane_count must be an integer >= 3")
        if throat_factor is not None and deceleration is not None:
            raise ValueError("throat_factor and deceleration each set the throat; give at most one")
        if outer_diameter_ratio is not None and length_ratio is not None:
            raise ValueError("outer_diameter_ratio and length_ratio each set the outer diameter; give at most one")
        if not roughness >= 0:
            raise ValueError("roughness must be nonnegative")
        self.flow, self.d2, self.b2, self.u2, self.c2u = map(float, (flow, d2, b2, u2, c2u))
        self.n, self.H, self.nq, self.rho = map(float, (n, H, nq, rho))
        self.impeller_blade_count = int(impeller_blade_count)
        self.vane_count = z = int(vane_count)
        self.roughness = float(roughness)
        self.rotation = parse_rotation(rotation)
        self.c2m = self.flow/(pi*self.d2*self.b2)

        least = minimum_gap_ratio(nq=self.nq, rho=self.rho, H=self.H)
        gap = least if gap_ratio is None else _positive(gap_ratio, "gap_ratio")
        if gap <= 1:
            raise ValueError("gap_ratio must exceed 1: the vanes stand outside the impeller")
        self.d3, self.b3 = gap*self.d2, inlet_width_ratio*self.b2
        self.b4 = outlet_width_ratio*self.b3
        r2, r3 = self.d2/2, self.d3/2
        if outer_diameter_ratio is not None:
            r4 = _positive(outer_diameter_ratio, "outer_diameter_ratio")*r2
            if r4 <= 1.05*r3:
                raise ValueError("outer_diameter_ratio leaves no room for a channel: make d4 at least 5% over d3")
        else:
            r4 = None
            stretch = 4.0 if length_ratio is None else _positive(length_ratio, "length_ratio")

        # Approach flow (T7.7.3 to T7.7.5). The leading edge is rounded to a
        # point, so it blocks nothing where the angle is taken.
        self.c3u = r2*self.c2u/r3
        self.c3m = self.flow/(pi*self.d3*self.b3)
        self.flow_angle = atan(self.c3m/self.c3u)
        self.c2 = float(np.hypot(self.c2u, self.c2m))

        # Throat (T7.7.8a, T7.7.10, T7.7.11).
        X = exp(self.flow/(self.b3*self.c2u*r2*z))
        free_vortex = r3*(X-1)
        if throat_factor is None:
            ratio = 0.8 if deceleration is None else _positive(deceleration, "deceleration")
            self.a3 = self.flow/(z*self.b3*ratio*self.c2)
        else:
            self.a3 = _positive(throat_factor, "throat_factor")*free_vortex
        self.throat_factor = self.a3/free_vortex
        self.c3q = self.flow/(z*self.a3*self.b3)

        pitch = 2*pi/z
        e = self.vane_thickness = vane_thickness_ratio*self.d2
        nose = np.array([r3*cos(pitch), r3*sin(pitch)])       # the next vane's leading edge

        # The leading spiral, at the angle whose throat to the next leading
        # edge is the design throat.
        def spiral(angle):
            theta = np.linspace(0, 1.5*pitch, 600)
            radius = r3*np.exp(theta*tan(angle))
            return np.column_stack((radius*np.cos(theta), radius*np.sin(theta)))

        def throat_width(angle):
            return float(np.min(np.linalg.norm(_offset(spiral(angle), e/2)-nose, axis=1)))

        low, high = radians(0.2), radians(60.)
        if not throat_width(low) < self.a3 < throat_width(high):
            raise ValueError("No leading spiral draws this throat past vanes this thick; "
                             "use fewer or thinner vanes, or a smaller throat")
        from scipy.optimize import brentq
        self.vane_angle = brentq(lambda a: throat_width(a)-self.a3, low, high, xtol=1e-12)
        lead = spiral(self.vane_angle)
        at = int(np.argmin(np.linalg.norm(_offset(lead, e/2)-nose, axis=1)))
        start_theta = float(np.arctan2(lead[at, 1], lead[at, 0]))
        start_radius = float(np.hypot(*lead[at]))

        # The channel: from the throat its width a grows linearly along its
        # length at the rate that reaches the optimum area ratio (T1.7.13),
        # and the vanes turn so that pitch*r*sin(alpha) - e is that width.
        self.equivalent_radius = sqrt(self.a3*self.b3/pi)
        fixed = r4

        def march(length, overshoot):
            ratio = 1.05+0.184*length/self.equivalent_radius
            exit_width = ratio*self.a3*self.b3/self.b4
            grow = (overshoot*exit_width-self.a3)/length
            steps = 2000
            ds = (length if fixed is None else 4*(fixed-r3)/radians(5)+length)/steps
            radius, width, s_along = [start_radius], [self.a3], [0.]
            angle = [max(self.vane_angle, np.arcsin(min((self.a3+e)/(start_radius*pitch), 1.)))]
            for _ in range(steps):
                sine = (width[-1]+e)/(radius[-1]*pitch)
                if sine >= 0.97:
                    raise ValueError("The channel would have to run radially to keep widening; "
                                     "use more vanes or a shorter channel")
                a_here = max(np.arcsin(sine), angle[-1])
                radius.append(radius[-1]+sin(a_here)*ds)
                width.append(width[-1]+grow*ds)
                s_along.append(s_along[-1]+ds)
                angle.append(a_here)
                if fixed is None and s_along[-1] >= length-1e-12:
                    break
                if fixed is not None and radius[-1] >= fixed:
                    # Land the last step on the outer diameter exactly.
                    f = (fixed-radius[-2])/(radius[-1]-radius[-2])
                    radius[-1] = fixed
                    width[-1] = width[-2]+f*(width[-1]-width[-2])
                    s_along[-1] = s_along[-2]+f*ds
                    break
            return ratio, exit_width, np.array(radius), np.array(angle), np.array(s_along)

        # The width law takes both walls at one radius, but the vane on the
        # outside of a channel is further out and steeper than that, so the
        # drawn exit comes out narrower than the law says. Measure it and
        # open the law until the drawn exit is the design exit.
        def draw(overshoot):
            if fixed is None:
                length = stretch*self.a3
                ratio, exit_width, radius, angle, along = march(length, overshoot)
            else:
                length = 4*self.a3
                for _ in range(60):
                    ratio, exit_width, radius, angle, along = march(length, overshoot)
                    if radius[-1] < fixed-1e-9:
                        raise ValueError("The channel cannot reach this outer diameter")
                    updated = float(along[-1])
                    if abs(updated-length) < 1e-10:
                        break
                    length = updated
            r4 = float(radius[-1])
            self.d4 = 2*r4
            self.channel_length, self.area_ratio, self.a4 = float(along[-1]), ratio, exit_width
            self.half_angle = atan((exit_width-self.a3)/(2*self.channel_length))
            self.exit_angle = float(angle[-1])

            # Camber: the leading spiral out to the throat, then the channel law.
            channel_theta = start_theta+np.r_[0., np.cumsum(np.diff(radius)/(0.5*(radius[1:]+radius[:-1])
                                                                             * np.tan(0.5*(angle[1:]+angle[:-1]))))]
            theta = np.r_[np.linspace(0, start_theta, 80)[:-1], channel_theta]
            rad = np.r_[r3*np.exp(np.linspace(0, start_theta, 80)[:-1]*tan(self.vane_angle)), radius]
            camber = _resample(np.column_stack((rad*np.cos(theta), rad*np.sin(theta))), 241)
            heading = theta[-1]+pi/2-self.exit_angle-(start_theta+pi/2-angle[0])
            self.deflection = abs(heading)
            self.pressure_recovery = (0.36*(self.channel_length/self.equivalent_radius)**0.26
                                      * (1-0.13*degrees(self.deflection)/30))
            self.c4 = self.flow/(z*self.a4*self.b4)
            self.outlet = DiffuserOutlet(r4, self.b4, self.flow/(2*pi*r4*self.b4),
                                         self.c4*cos(self.exit_angle))

            # Faces: the camber either side by half the thickness, rounded to a
            # point at both ends so the two faces share their end stations.
            arc = np.r_[0., np.cumsum(np.linalg.norm(np.diff(camber, axis=0), axis=1))]
            ramp = 1.5*e
            half = 0.5*e*np.sqrt(np.clip(np.minimum(arc, arc[-1]-arc)/ramp, 0, 1))
            inner_face, outer_face = _offset(camber, half), _offset(camber, -half)
            # Channel 0 lies between vane 0's inner face and vane 1's outer face.
            neighbour = _rotate(outer_face, pitch)
            facing = inner_face[np.arctan2(inner_face[:, 1], inner_face[:, 0]) > pitch+1e-6]
            gaps = _signed_distance(facing, neighbour)
            if len(facing) == 0 or np.any(gaps >= 0):
                raise ValueError("Adjacent vanes cross: use fewer or thinner vanes")
            drawn_throat = float(np.min(np.linalg.norm(inner_face-nose, axis=1)))
            drawn_exit = float(abs(_signed_distance(inner_face[-2:-1], neighbour)[0]))
            self._faces = (inner_face, outer_face)

            return drawn_throat, drawn_exit

        overshoot = 1.
        for _ in range(30):
            drawn_throat, drawn_exit = draw(overshoot)
            if abs(drawn_exit/self.a4-1) < 1e-4:
                break
            overshoot *= (self.a4-self.a3)/max(drawn_exit-self.a3, 1e-9)*0.7+0.3
        conflicts = vane_count_conflicts(z, self.impeller_blade_count)
        allowed = degrees(radians(16.5)*sqrt(self.equivalent_radius/self.channel_length))
        equivalent = degrees(atan(self.equivalent_radius/self.channel_length*(sqrt(self.area_ratio)-1)))
        self.checks = DiffuserChecks(
            gap, least, inlet_width_ratio, self.throat_factor, self.c3q/self.c2,
            self.channel_length/self.a3, self.b3/self.a3, self.c4**2/(2*G0*self.H),
            equivalent, allowed, conflicts, degrees(self.vane_angle-self.flow_angle),
            degrees(self.deflection), drawn_throat/self.a3, drawn_exit/self.a4)
        self.geometry = None

    @classmethod
    def from_impeller(cls, impeller, *, rho=None, **design):
        """Design the diffuser for an :class:`~pyskyfire.pump.Impeller`'s outlet.

        The swirl is the impeller's slip prediction on its mean streamline, and
        the head its outlet head. ``rho`` defaults to the impeller's.
        """
        rho = impeller.rho if rho is None else rho
        if rho is None:
            raise ValueError("The gap rule needs a density: give rho or build the impeller with one")
        d = impeller.dimensions
        return cls(flow=d.internal_flow, d2=d.d2, b2=d.b2, u2=d.u2,
                   c2u=impeller.velocity_triangles.outlet[1].cu,
                   impeller_blade_count=d.blade_count, n=impeller.n, H=impeller.H, nq=d.nq,
                   rho=rho, rotation=impeller.rotation, **design)

    def build_geometry(self, *, chord_points=65):
        """Resolve the vane surfaces and walls for display. Returns and stores them."""
        if isinstance(chord_points, bool) or not isinstance(chord_points, (int, np.integer)) or chord_points < 9:
            raise ValueError("chord_points must be an integer >= 9")
        sign = rotation_sign(self.rotation)
        r3, r4 = self.d3/2, self.d4/2

        def width(r):
            return self.b3+(self.b4-self.b3)*np.clip((r-r3)/(r4-r3), 0, 1)

        sides = []
        stations = np.linspace(0, len(self._faces[0])-1, chord_points)
        for face in self._faces:
            # Both faces are sampled at the same camber stations, so the end
            # caps join matching points.
            plane = np.column_stack([np.interp(stations, np.arange(len(face)), face[:, k]) for k in range(2)])
            r = np.hypot(plane[:, 0], plane[:, 1])
            theta = sign*np.arctan2(plane[:, 1], plane[:, 0])
            half = width(r)/2
            sides.append(np.stack([np.column_stack((s*half, r, theta)) for s in (-1, 1)]))
        radii = np.r_[self.d2/2, np.linspace(r3, r4, 33)]
        walls = [np.column_stack((s*width(radii)/2, radii)) for s in (1, -1)]
        self.geometry = DiffuserGeometry(BladeSurface(*sides), self.vane_count, *walls, sign)
        return self.geometry

    def solve(self, *, inlet_state, mass_flow, transport):
        """Total pressure lost from impeller outlet to channel exit, Table 3.8 (2).

        The swirl entering is the design swirl and the through-flow follows the
        mass flow, as the volute treats it. Density and viscosity are taken at
        the inlet; the flow is a liquid.
        """
        from pyskyfire.ducts.pipe import _property, _state
        T, p = float(inlet_state.T), float(inlet_state.p)
        rho = _property(transport, "rho", T, p)
        nu = _property(transport, "mu", T, p)/rho
        flow = float(mass_flow)/rho
        z, r2, r3 = self.vane_count, self.d2/2, self.d3/2
        c2m = flow/(pi*self.d2*self.b2)
        c2 = float(np.hypot(self.c2u, c2m))
        c3q = flow/(z*self.a3*self.b3)
        c4 = flow/(z*self.a4*self.b4)
        angle = atan((flow/(pi*self.d3*self.b3))/self.c3u)
        # T3.8.23 over the vaneless gap only, r2 to r3, at its own width.
        path = (r3-r2)/sin(angle)
        cf = _cf(c2, path, nu, self.roughness)
        vaneless = 2*cf*r2/(self.b3*sin(angle)*cos(angle)**2)*(self.c2u/self.u2)**2*(1-r2/r3)
        # T3.8.18, the triangular inlet region, lengths over d2.
        a3s, b3s, b2s = self.a3/self.d2, self.b3/self.d2, self.b2/self.d2
        phi = c2m/self.u2
        cf = _cf(c3q, r3*2*pi/z, nu, self.roughness)
        inlet = ((cf+0.0015)*(a3s+b3s)*pi**3*(phi*b2s)**2/(8*(z*a3s*b3s)**3)
                 * (1+c2/c3q)**3)
        head = (c3q/self.u2)**2
        shock = head*0.3*(c2/c3q-1)**2
        channel = head*(1-self.pressure_recovery-1/self.area_ratio**2)
        jet = c4*sin(self.exit_angle)
        mixed = flow/(2*pi*self.d4/2*self.b4)
        mixing = ((jet-mixed)/self.u2)**2
        total = vaneless+inlet+shock+channel+mixing
        loss = total*rho*self.u2**2/2
        if loss >= p:
            raise RuntimeError(f"diffuser pressure fell to {p-loss:g} Pa")
        return DiffuserResult(vaneless, inlet, shock, channel, mixing, total, loss, c2, c3q, c4,
                              rho, float(mass_flow), (_state(T, p, rho), _state(T, p-loss, rho)))

    def __repr__(self):
        return (f"VanedDiffuser(vanes={self.vane_count}, d3={self.d3:.4g}, d4={self.d4:.4g}, "
                f"b3={self.b3:.4g}, a3={self.a3:.4g}, a4={self.a4:.4g}, L={self.channel_length:.4g})")
