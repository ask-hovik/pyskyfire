"""Preliminary radial impeller sizing and physical inputs for pyTurbo Aero.

No surface generator or efficiency predictor is run during construction.
Metal angles are measured from the backward circumferential direction toward
the meridional direction. Positive wrap denotes backward sweep for +x rotation.

The outlet triangle carries slip, and the design is closed on it. The flow does
not leave along the blades, so the swirl a wheel actually delivers is less than
its blade angle says, and it is that swirl the head comes from. Exactly one of
the blade outlet angle and the outlet diameter is therefore solved, so that the
swirl the blades are predicted to leave (Gülich Eq. T3.2.7) is the swirl the
head needs at the supplied hydraulic efficiency (Eq. T3.2.8). Solving the angle
is the order of Gülich's own design procedure; solving the diameter is its
converse, for when the blades are fixed.

The hydraulic efficiency is the impeller's alone: the losses between its inlet
and its outlet, Gülich's eta_h,La (Eq. 8.13a). The head ``H`` is therefore the
total head at the impeller outlet, which is the pump's delivered head plus
whatever the volute or diffuser after it will lose. Those losses belong to the
collector and are charged there.

References
----------
Gülich, J. F., *Centrifugal Pumps*, 4th ed., Springer, 2020. Table 3.2,
Eq. (3.26), Eq. (7.1) and chapter 7.2.
"""
from dataclasses import dataclass
from inspect import Parameter, signature
from math import pi

import numpy as np
from scipy.optimize import brentq

from pyskyfire.common.rotation import COUNTERCLOCKWISE, parse_rotation
from .geometry import (BladeSection, BladeSettings, ImpellerBladeDefinition,
                       MeridionalContours, SplitterSettings)

G0 = 9.80665


def _positive(value, name, *, zero=False):
    if not np.isscalar(value) or not np.isfinite(value) or (value < 0 if zero else value <= 0):
        raise ValueError(f"{name} must be finite and {'nonnegative' if zero else 'positive'}")
    return float(value)


def design_value(value, *, label, **context):
    """Dispatch a scalar callback by its requested keyword names.

    Angle/prewhirl callbacks receive span_fraction, x, r. Thickness callbacks
    additionally receive chord_fraction (normalized m-prime along the blade).
    Callbacks are evaluated at scalar positions and must return finite scalars.
    """
    if callable(value):
        kwargs = {}
        for name, p in signature(value).parameters.items():
            if p.kind == Parameter.VAR_KEYWORD:
                kwargs.update(context)
            elif p.kind in (Parameter.POSITIONAL_ONLY, Parameter.VAR_POSITIONAL):
                raise ValueError(f"{label}: use named arguments, not positional-only arguments or *args")
            elif name in context:
                kwargs[name] = context[name]
            elif p.default is Parameter.empty:
                raise ValueError(f"{label}: unknown argument {name!r}; available: {', '.join(context)}")
        value = value(**kwargs)
    if not np.isscalar(value) or not np.isfinite(value):
        raise ValueError(f"{label} must return a finite scalar")
    return float(value)


def specific_speed(n, Q, H):
    """Gülich nq convention: rpm, m³/s, m."""
    return n*np.sqrt(Q)/H**0.75


def pressure_coefficient(nq, *, fT=1.0):
    """Best-efficiency head coefficient of radial pumps, Gülich Eq. (3.26).

    A mean curve through measured industrial pumps, so it proposes an outlet
    diameter rather than fixing one. ``fT`` says where in the scatter band the
    design sits: Gülich finds 1.0 describes multistage pumps and 1.1 low-head
    single-stage volute pumps (EN 733). A higher coefficient gives a smaller
    wheel and a flatter, less stable head curve.
    """
    return 1.21*fT*np.exp(-0.77*nq/100)


def outlet_diameter(psi, n, H):
    return 60/(pi*n)*np.sqrt(2*G0*H/psi)


def outlet_width_ratio(nq):
    """Mean relative outlet width b2/d2 of radial impellers, Gülich Eq. (7.1)."""
    v = nq/100
    return 0.017+0.262*v-0.080*v*v+0.0093*v**3


def inlet_diameter_from_phi(*, Q, n, phi1, inlet_hub_diameter):
    """Solve annular continuity for a fixed wetted hub and Cm1 = phi1*Ueye."""
    for name, v in (("Q", Q), ("n", n), ("phi1", phi1), ("inlet_hub_diameter", inlet_hub_diameter)):
        _positive(v, name)
    d = inlet_hub_diameter
    k = 240*Q/(phi1*pi*pi*n)
    return brentq(lambda eye: eye*(eye*eye-d*d)-k, d, d+k**(1/3), xtol=1e-14)


@dataclass(frozen=True)
class OutletSlip:
    """Gülich's slip prediction at the impeller outlet, Table 3.2.

    ``blockage`` is tau_2 (T3.2.3), ``limiting_ratio`` epsilon_lim (T3.2.4),
    ``diameter_correction`` k_w (T3.2.5) and ``slip_factor`` gamma (T3.2.6).
    ``cm`` is the unblocked meridional velocity and ``cu`` the absolute swirl
    the blades leave (T3.2.7), both in m/s.
    """
    blade_angle_deg: float
    blockage: float
    limiting_ratio: float
    diameter_correction: float
    slip_factor: float
    cm: float
    cu: float


def outlet_slip(*, beta2_deg, blade_count, d2, b2, u2, flow, blade_thickness, d1m, f1=0.98):
    """Swirl a radial wheel leaves with, including slip (Gülich Table 3.2).

    ``flow`` is the flow through the impeller, leakage included. ``d1m`` is the
    mean leading-edge diameter, which lowers the slip factor once it reaches
    far enough out that the blades are too short to guide the flow.
    ``blade_thickness`` is the full body thickness: Gülich correlated the
    factor with it as the blockage, not with a thinned trailing edge. ``f1``
    is 0.98 for radial impellers. The correlation reproduces tests with a
    standard deviation of about 4 percent.
    """
    angle = np.radians(beta2_deg)
    margin = 1-blade_thickness*blade_count/(pi*d2*np.sin(angle))
    if margin <= 0:
        raise ValueError(f"At a {beta2_deg:.3g} degree outlet angle, {blade_count} blades "
                         f"{blade_thickness*1e3:.3g} mm thick fill the whole outlet circumference")
    limit = float(np.exp(-8.16*np.sin(angle)/blade_count))
    ratio = d1m/d2
    kw = 1. if ratio <= limit else 1-((ratio-limit)/(1-limit))**3
    gamma = f1*(1-np.sqrt(np.sin(angle))/blade_count**0.7)*kw
    cm = flow/(pi*d2*b2)
    cu = u2*gamma-cm/(margin*np.tan(angle))
    return OutletSlip(float(beta2_deg), float(1/margin), limit, float(kw), float(gamma), float(cm), float(cu))


@dataclass(frozen=True)
class ImpellerDimensions:
    nq: float
    psi: float
    d1: float
    inlet_hub_diameter: float
    shaft_diameter: float
    d2: float
    b2: float
    u1: float
    u2: float
    internal_flow: float
    blade_count: int


@dataclass(frozen=True)
class VelocityTriangle:
    span_fraction: float
    radius: float
    u: float
    cm: float
    cu: float
    w: float
    beta_flow_deg: float
    beta_metal_deg: float


@dataclass(frozen=True)
class ImpellerVelocityTriangles:
    inlet: tuple[VelocityTriangle, ...]
    outlet: tuple[VelocityTriangle, ...]
    euler_head: float
    hydraulic_head: float
    # The outlet triangles carry slip: cu is Gülich's prediction and
    # beta_flow_deg the unblocked relative flow angle it implies (T3.2.15).


@dataclass(frozen=True)
class ImpellerDiagnostics:
    """How the closed design sits against Gülich's guidance. Nothing is enforced.

    ``solved`` names what the head was closed on, ``"beta2"`` or ``"d2"``.
    ``psi`` is the resolved head coefficient 2gH/u2² and ``psi_correlation``
    what Eq. (3.26) proposes at this specific speed and ``fT``. ``slip`` is the
    mean-streamline prediction. ``deviation_deg`` is beta2B minus the flow
    angle, which Gülich suggests keeping below 10 to 14 degrees, and
    ``deceleration_ratio`` is w2 over w1 at the shroud, whose de Haller limit
    he puts at 0.7.
    """
    solved: str
    psi: float
    psi_correlation: float
    slip: OutletSlip
    deviation_deg: float
    deceleration_ratio: float


def meridional_contours(d1, dh, d2, b2, *, axial_length=None):
    """Smooth quarter-ellipse endwalls with axial inlet and radial outlet.

    The inlet is x=0; the shroud outlet is x=L and hub outlet x=L+b2.
    Axial length is a design choice, with a preliminary dimension-based default.
    """
    L = axial_length if axial_length is not None else 0.75*(d1-dh)/2+b2
    _positive(L, "axial_length")
    t = np.linspace(0, pi/2, 129)
    hub = np.column_stack(((L+b2)*np.sin(t), dh/2+(d2-dh)/2*(1-np.cos(t))))
    shroud = np.column_stack((L*np.sin(t), d1/2+(d2-d1)/2*(1-np.cos(t))))
    return MeridionalContours(hub, shroud)


def _section_path(passage, settings, span):
    """Normalised m' fraction and (x, r) samples along one blade section."""
    s = np.linspace(settings.start_hub, settings.end, 513)
    xr = passage.sample(s, span)
    dm = np.linalg.norm(np.diff(xr, axis=0), axis=1)
    mp = np.r_[0., np.cumsum(dm/(0.5*(xr[:-1,1]+xr[1:,1])))]
    return mp, xr


class Impeller:
    """Size a radial impeller; call build_geometry() to generate pyTurbo surfaces.

    Required lengths are in metres. ``blade_thickness`` is total section-normal
    body thickness: a float or function of chord_fraction, span_fraction, x, r.
    Angle and prewhirl functions accept span_fraction, x, r. Span runs hub to
    blade tip, including when tip clearance is present. ``phi1`` refers to the
    unblocked inlet annulus, and d1 overrides its sizing equation.

    ``H`` is the total head at the impeller outlet and ``hydraulic_efficiency``
    the impeller's own (see the module notes); it is a supplied assumption,
    never predicted here. The design is closed on one quantity. By default
    ``beta2_deg`` is solved and the outlet diameter comes from one of ``psi``,
    ``d2``, ``u2`` or a supplied ``meridional``, or failing all of them from
    Eq. (3.26) at ``fT``. Give ``beta2_deg`` instead and ``d2`` is solved. Giving
    both over-specifies the wheel. The solved angle is the one below the peak
    of the swirl it can deliver: at low specific speed that peak is broad and
    flat, and a head above it cannot be had at any angle, which is reported
    with the head coefficient that can.

    ``blade_count`` is always an input. Gülich rejects empirical blade-count
    formulas as likely to violate the criteria that actually set it: pressure
    pulsations against the collector, blade loading, head-curve stability.
    Splitters do not count toward it in the slip prediction.

    ``incidence_deg`` means metal angle minus unblocked relative inlet flow
    angle. No leakage is assumed for volumetric_efficiency=1. ``b2`` defaults
    to Eq. (7.1).

    ``rotation`` is which way the wheel turns about its own +x shaft, and only
    mirrors the blades: every velocity triangle is a magnitude in the rotating
    frame and is the same either way. It has to match the sense of the volute
    collecting the discharge, once both are placed.
    Create a new Impeller to change design inputs; build_geometry can be repeated
    at different numerical resolutions using the same resolved definition.
    """
    def __init__(self, Q, H, n, *, shaft_diameter, inlet_hub_diameter,
                 blade_count, blade_thickness, hydraulic_efficiency,
                 beta2_deg=None, psi=None, d2=None, u2=None, fT=1.0, b2=None,
                 phi1=0.20, d1=None, beta1_deg=None, incidence_deg=0.0, inlet_prewhirl=0.0,
                 wrap_deg=None, volumetric_efficiency=1.0,
                 rho=None, tip_clearance=0.0, meridional=None,
                 axial_length=None, blade_settings=None, splitter=None,
                 rotation=COUNTERCLOCKWISE):
        self.Q = _positive(Q, "Q")
        self.H = _positive(H, "H")
        self.n = _positive(n, "n")
        dh = _positive(inlet_hub_diameter, "inlet_hub_diameter")
        ds = _positive(shaft_diameter, "shaft_diameter")
        _positive(phi1, "phi1")
        _positive(fT, "fT")
        if dh < ds:
            raise ValueError("inlet_hub_diameter must be at least shaft_diameter")
        if isinstance(blade_count, bool) or not isinstance(blade_count, (int, np.integer)) or blade_count < 3:
            raise ValueError("blade_count must be an integer >= 3, the range of Gülich's slip correlation")
        for label, efficiency in (("volumetric_efficiency", volumetric_efficiency),
                                  ("hydraulic_efficiency", hydraulic_efficiency)):
            if not np.isscalar(efficiency) or not 0 < efficiency <= 1:
                raise ValueError(f"{label} must lie in (0,1]")
        self.volumetric_efficiency = float(volumetric_efficiency)
        self.hydraulic_efficiency = float(hydraulic_efficiency)
        self.rho = None if rho is None else _positive(rho, "rho")
        tip_clearance = _positive(tip_clearance, "tip_clearance", zero=True)
        settings = blade_settings or BladeSettings()
        if not isinstance(settings, BladeSettings):
            raise TypeError("blade_settings must be BladeSettings")
        if splitter is not None and not isinstance(splitter, SplitterSettings):
            raise TypeError("splitter must be SplitterSettings")
        if meridional is not None and not isinstance(meridional, MeridionalContours):
            raise TypeError("meridional must be MeridionalContours")
        if meridional is not None and axial_length is not None:
            raise ValueError("Specify meridional or axial_length, not both")
        sizes = [k for k, v in (("psi", psi), ("d2", d2), ("u2", u2), ("meridional", meridional))
                 if v is not None]
        if len(sizes) > 1:
            raise ValueError(f"{' and '.join(sizes)} each set the outlet diameter; give at most one")
        if beta2_deg is not None and sizes:
            raise ValueError(f"beta2_deg and {sizes[0]} together over-specify the wheel: the head "
                             "is closed on one of them, so give the blade angle or the outlet size")
        qi = Q/volumetric_efficiency
        nq = specific_speed(n, Q, H)
        d1 = inlet_diameter_from_phi(Q=qi, n=n, phi1=phi1, inlet_hub_diameter=dh) if d1 is None else _positive(d1, "d1")
        if b2 is not None:
            _positive(b2, "b2")
        # Every blade body is checked before anything reads it.
        thickness = lambda **ctx: design_value(blade_thickness, label="blade_thickness", **ctx)

        def layout(d2):
            """Passage, inlet triangles and slip inputs for one outlet diameter."""
            width = (b2 if b2 is not None else
                     float(meridional.hub[-1,0]-meridional.shroud[-1,0]) if meridional is not None else
                     outlet_width_ratio(nq)*d2)
            _positive(d2, "calculated d2")
            _positive(width, "calculated b2")
            if not dh < d1 < d2:
                raise ValueError("Radial impeller requires inlet_hub_diameter < d1 < d2; change duty, phi1, or diameters")
            mer = meridional if meridional is not None else meridional_contours(d1, dh, d2, width, axial_length=axial_length)
            if not np.allclose([mer.hub[0,1], mer.shroud[0,1], mer.hub[-1,1], mer.shroud[-1,1],
                                mer.hub[-1,0]-mer.shroud[-1,0]], [dh/2, d1/2, d2/2, d2/2, width], rtol=1e-6, atol=1e-10):
                raise ValueError("Meridional inlet radii, outlet radii, and outlet width must match sizing")
            if not np.isclose(mer.hub[0,0], mer.shroud[0,0], atol=1e-10):
                raise ValueError("Meridional inlet must be in one axial plane")
            mer = mer.with_leading_edge(settings.start_hub, settings.start_tip)
            active = mer.with_clearance(tip_clearance)
            # The blade inlet triangle belongs on the leading edge, which a slanted
            # one puts well downstream of the eye. Its meridional velocity follows
            # from unblocked continuity through that line, counting the component of
            # the through-flow that crosses it: a slanted edge stands obliquely to
            # the meridional direction and passes less than the area it sweeps.
            # Eye sizing keeps using the inlet plane and phi1.
            # Take the through-flow direction from the two wall tangents rather than
            # from the span lines: a wall tangent is a property of the passage, so
            # the area a leading edge passes does not depend on how the stations
            # were paired to place it.
            spans = np.linspace(0., 1., 65)
            edge = mer.sample(settings.start_hub, spans)
            walls = mer.wall_tangents(settings.start_hub)
            flow = (1-spans)[:,None]*walls[0]+spans[:,None]*walls[1]
            flow /= np.linalg.norm(flow, axis=1)[:,None]
            span_line = edge[-1]-edge[0]
            normal = np.array([span_line[1], -span_line[0]])/np.linalg.norm(span_line)
            cm1 = qi/(2*pi*np.linalg.norm(span_line)*np.trapezoid(edge[:,1]*abs(flow@normal), spans))
            inlet = []
            for span in (0.0, 0.5, 1.0):
                x, r = active.sample(settings.start_hub, span)
                ctx = dict(span_fraction=span, x=float(x), r=float(r))
                u1 = 2*pi*n/60*r
                cu1 = design_value(inlet_prewhirl, label="inlet_prewhirl", **ctx)
                angle = float(np.degrees(np.arctan2(cm1, u1-cu1)))
                metal = (angle + design_value(incidence_deg, label="incidence_deg", **ctx)
                         if beta1_deg is None else design_value(beta1_deg, label="beta1_deg", **ctx))
                if not 0 < metal < 90:
                    raise ValueError("Supported backward-swept metal angles: 0 < beta1 < 90, 0 < beta2 <= 90 degrees")
                inlet.append(VelocityTriangle(span, r, u1, cm1, cu1, np.hypot(cm1,u1-cu1), angle, metal))
            mp, xr = _section_path(active, settings, .5)
            body = max(_positive(thickness(chord_fraction=float(c), span_fraction=.5,
                                           x=float(np.interp(c, mp/mp[-1], xr[:,0])),
                                           r=float(np.interp(c, mp/mp[-1], xr[:,1]))), "blade_thickness")
                       for c in np.linspace(0, 1, 17))
            return dict(d2=d2, b2=width, mer=mer, active=active, inlet=inlet, thickness=body,
                        # Gülich's mean inlet diameter: the RMS of the leading
                        # edge's hub and tip, so it follows a slanted edge out.
                        d1m=2*np.sqrt((edge[0,1]**2+edge[-1,1]**2)/2))

        def outlet_angle(geometry, span):
            x2, r2 = geometry["active"].sample(1., span)
            return design_value(beta2_deg, label="beta2_deg", span_fraction=span, x=float(x2), r=float(r2))

        def slip(geometry, angle):
            d2 = geometry["d2"]
            if not 0 < angle <= 90:
                raise ValueError("Supported backward-swept metal angles: 0 < beta1 < 90, 0 < beta2 <= 90 degrees")
            return outlet_slip(beta2_deg=angle, blade_count=int(blade_count), d2=d2, b2=geometry["b2"],
                               u2=pi*n*d2/60, flow=qi, blade_thickness=geometry["thickness"],
                               d1m=geometry["d1m"])

        def required_swirl(geometry):
            # Eq. (T3.2.8): the swirl this head takes, on the mean streamline.
            u2 = pi*n*geometry["d2"]/60
            mean = geometry["inlet"][1]
            return G0*H/(hydraulic_efficiency*u2)+mean.u*mean.cu/u2

        eta = hydraulic_efficiency
        if beta2_deg is None:
            solved = "beta2"
            outer = (2*meridional.shroud[-1,1] if meridional is not None else
                     _positive(d2, "d2") if d2 is not None else
                     60*_positive(u2, "u2")/(pi*n) if u2 is not None else
                     outlet_diameter(_positive(psi, "psi") if psi is not None else
                                     pressure_coefficient(nq, fT=fT), n, H))
            geometry = layout(float(outer))
            target = required_swirl(geometry)
            speed = pi*n*geometry["d2"]/60
            lowest = np.degrees(np.arcsin(min(1., geometry["thickness"]*blade_count/(pi*geometry["d2"]))))
            if lowest >= 89.9:
                raise ValueError("The blades are too thick to fit round this outlet at any angle")
            angles = np.linspace(lowest+1e-3*(90-lowest), 90., 1801)
            swirl = np.array([slip(geometry, a).cu for a in angles])
            peak = int(np.argmax(swirl))
            if target > swirl[peak]:
                # psi = 2 eta (cu2 - u1 cu1/u2)/u2 at the most this wheel can do.
                mean = geometry["inlet"][1]
                reach = 2*eta*(swirl[peak]-mean.u*mean.cu/speed)/speed
                raise ValueError(
                    f"No blade outlet angle makes this head: with {blade_count} blades and a hydraulic "
                    f"efficiency of {eta:.3g} the wheel reaches at most psi = {reach:.4g}, near "
                    f"beta2 = {angles[peak]:.3g} degrees, and psi = {2*G0*H/speed**2:.4g} was asked. "
                    "Raise blade_count, lower psi (a larger wheel), or revisit hydraulic_efficiency")
            beta2_deg = float(brentq(lambda a: slip(geometry, a).cu-target, angles[0], angles[peak], xtol=1e-10))
        else:
            solved = "d2"
            def residual(outer):
                geometry = layout(outer)
                return slip(geometry, outlet_angle(geometry, .5)).cu-required_swirl(geometry)
            low = d1*1.05
            if residual(low) > 0:
                raise ValueError("Even an outlet just outside the eye makes more than this head; "
                                 "lower beta2_deg or blade_count")
            high = outlet_diameter(pressure_coefficient(nq, fT=fT), n, H)
            for _ in range(40):
                if high > low and residual(high) > 0:
                    break
                high *= 1.25
            else:
                raise ValueError(f"No outlet diameter makes this head at beta2 = {beta2_deg} degrees "
                                 f"with {blade_count} blades and a hydraulic efficiency of {eta:.3g}")
            geometry = layout(brentq(residual, low, high, xtol=1e-12))
        d2, b2 = geometry["d2"], geometry["b2"]
        mer, active = geometry["mer"], geometry["active"]
        u2 = pi*n*d2/60
        self.dimensions = ImpellerDimensions(nq, 2*G0*H/u2**2, d1, dh, ds, d2, b2,
                                             pi*n*d1/60, u2, qi, int(blade_count))
        inlet, outlet, sections = geometry["inlet"], [], []
        for triangle in inlet:
            span = triangle.span_fraction
            ctx = dict(span_fraction=span, x=float(active.sample(settings.start_hub, span)[0]), r=triangle.radius)
            angle = outlet_angle(geometry, span)
            predicted = slip(geometry, angle)
            x2, r2 = active.sample(1., span)
            flow_angle = float(np.degrees(np.arctan2(predicted.cm, u2-predicted.cu)))
            outlet.append(VelocityTriangle(span, r2, u2, predicted.cm, predicted.cu,
                                           np.hypot(predicted.cm, u2-predicted.cu), flow_angle, angle))
            mp, xr = _section_path(active, settings, span)
            f = mp/mp[-1]
            f = 0.5*(1-np.cos(pi*f)) if settings.angle_law == "cosine" else f
            metal = triangle.beta_metal_deg
            beta = np.radians(metal+(angle-metal)*f)
            wrap = (np.degrees(np.trapezoid(1/np.tan(beta), mp)) if wrap_deg is None
                    else design_value(wrap_deg, label="wrap_deg", **ctx))
            _positive(wrap, "wrap_deg")
            sections.append(BladeSection(span, metal, angle, wrap))
            # Validate the thickness callback over the whole section before
            # loading the optional backend.
            for chord in np.linspace(0, 1, 17):
                _positive(thickness(chord_fraction=float(chord), span_fraction=span,
                                    x=float(np.interp(chord, f, xr[:,0])),
                                    r=float(np.interp(chord, f, xr[:,1]))), "blade_thickness")
        euler = float((outlet[1].u*outlet[1].cu-inlet[1].u*inlet[1].cu)/G0)
        self.velocity_triangles = ImpellerVelocityTriangles(tuple(inlet), tuple(outlet), euler, eta*euler)
        mean = slip(geometry, outlet[1].beta_metal_deg)
        self.diagnostics = ImpellerDiagnostics(
            solved, self.dimensions.psi, float(pressure_coefficient(nq, fT=fT)), mean,
            float(outlet[1].beta_metal_deg-outlet[1].beta_flow_deg), float(outlet[1].w/inlet[-1].w))
        self.rotation = parse_rotation(rotation)
        self.blade_definition = ImpellerBladeDefinition(mer, tuple(sections), int(blade_count),
                                                        blade_thickness, settings, splitter,
                                                        tip_clearance, self.rotation)
        self.geometry = None

    @classmethod
    def from_pressure_rise(cls, *, mdot, rho, dp, n, **kwargs):
        for label, v in (("mdot",mdot), ("rho",rho), ("dp",dp)):
            _positive(v, label)
        return cls(Q=mdot/rho, H=dp/(rho*G0), n=n, rho=rho, **kwargs)

    @classmethod
    def from_duty(cls, duty, **kwargs):
        """Size for a :class:`~pyskyfire.pump.PumpDuty`: its inlet flow, head and speed.

        The inlet density, where the duty knows one, is what the blade power is
        reported at.
        """
        return cls(Q=duty.Q, H=duty.H, n=duty.n, rho=duty.rho_in, **kwargs)

    @property
    def blade_power(self):
        """Power the blades put into the fluid [W], or None without ``rho``.

        Every kilogram through the impeller, leakage included, takes the full
        blade work g*H/eta_h. Disk friction and the bearings are charged by the
        components that own them, not here.
        """
        if self.rho is None:
            return None
        return self.rho*G0*self.dimensions.internal_flow*self.H/self.hydraulic_efficiency

    def build_geometry(self, *, span_points=33, chord_points=129):
        """Generate surfaces lazily. Returns and stores an ImpellerGeometry."""
        from ._pyturbo import build_geometry
        geometry = build_geometry(self.blade_definition, span_points=span_points, chord_points=chord_points)
        self.geometry = geometry
        return geometry

    def __repr__(self):
        d = self.dimensions
        return (f"Impeller(Q={self.Q:g}, H={self.H:g}, n={self.n:g}, "
                f"d1={d.d1:.6g}, d2={d.d2:.6g}, b2={d.b2:.6g}, blades={d.blade_count})")
