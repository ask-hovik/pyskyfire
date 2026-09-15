"""Preliminary radial impeller sizing and physical inputs for pyTurbo Aero.

No surface generator or efficiency predictor is run during construction.
Metal angles are measured from the backward circumferential direction toward
the meridional direction. Positive wrap denotes backward sweep for +x rotation.
"""
from dataclasses import dataclass
from inspect import Parameter, signature
from math import pi

import numpy as np
from scipy.optimize import brentq

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


def pressure_coefficient(nq, *, fT=1.10):
    """Preliminary correlation retained from the former implementation."""
    return 1.21*fT*np.exp(-0.77*nq/100)


def outlet_diameter(psi, n, H):
    return 60/(pi*n)*np.sqrt(2*G0*H/psi)


def outlet_width_ratio(nq):
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
    hydraulic_head: float | None
    # Outlet triangle is blade-congruent (no slip). These are diagnostics,
    # not an independent performance prediction of the requested pump duty.


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


class Impeller:
    """Size a radial impeller; call build_geometry() to generate pyTurbo surfaces.

    Required lengths are in metres. ``blade_thickness`` is total section-normal
    body thickness: a float or function of chord_fraction, span_fraction, x, r.
    Angle and prewhirl functions accept span_fraction, x, r. Span runs hub to
    blade tip, including when tip clearance is present. ``phi1`` refers to the
    unblocked inlet annulus, and d1 overrides its sizing equation.

    ``incidence_deg`` means metal angle minus unblocked relative inlet flow
    angle. No leakage is assumed for volumetric_efficiency=1. Hydraulic and
    mechanical efficiencies are supplied assumptions, never predicted here.
    Create a new Impeller to change design inputs; build_geometry can be repeated
    at different numerical resolutions using the same resolved definition.
    """
    def __init__(self, Q, H, n, *, shaft_diameter, inlet_hub_diameter,
                 blade_count, beta2_deg, blade_thickness, phi1=0.20,
                 d1=None, beta1_deg=None, incidence_deg=0.0, inlet_prewhirl=0.0,
                 wrap_deg=None, volumetric_efficiency=1.0,
                 hydraulic_efficiency=None, mechanical_efficiency=None,
                 rho=None, tip_clearance=0.0, meridional=None,
                 axial_length=None, blade_settings=None, splitter=None):
        self.Q = _positive(Q, "Q")
        self.H = _positive(H, "H")
        self.n = _positive(n, "n")
        dh = _positive(inlet_hub_diameter, "inlet_hub_diameter")
        ds = _positive(shaft_diameter, "shaft_diameter")
        _positive(phi1, "phi1")
        if dh < ds:
            raise ValueError("inlet_hub_diameter must be at least shaft_diameter")
        if isinstance(blade_count, bool) or not isinstance(blade_count, (int, np.integer)) or blade_count < 2:
            raise ValueError("blade_count must be an integer >= 2")
        for label, efficiency in (("volumetric_efficiency", volumetric_efficiency),
                                  ("hydraulic_efficiency", hydraulic_efficiency),
                                  ("mechanical_efficiency", mechanical_efficiency)):
            if efficiency is not None and not 0 < efficiency <= 1:
                raise ValueError(f"{label} must lie in (0,1]")
        if volumetric_efficiency is None:
            raise ValueError("volumetric_efficiency must lie in (0,1]")
        self.volumetric_efficiency = volumetric_efficiency
        self.hydraulic_efficiency = hydraulic_efficiency
        self.mechanical_efficiency = mechanical_efficiency
        self.rho = None if rho is None else _positive(rho, "rho")
        tip_clearance = _positive(tip_clearance, "tip_clearance", zero=True)
        settings = blade_settings or BladeSettings()
        if not isinstance(settings, BladeSettings):
            raise TypeError("blade_settings must be BladeSettings")
        if splitter is not None and not isinstance(splitter, SplitterSettings):
            raise TypeError("splitter must be SplitterSettings")
        qi = Q/volumetric_efficiency
        nq = specific_speed(n, Q, H)
        psi = pressure_coefficient(nq)
        d2 = outlet_diameter(psi, n, H)
        b2 = outlet_width_ratio(nq)*d2
        _positive(d2, "calculated d2")
        _positive(b2, "calculated b2")
        d1 = inlet_diameter_from_phi(Q=qi, n=n, phi1=phi1, inlet_hub_diameter=dh) if d1 is None else _positive(d1, "d1")
        if not dh < d1 < d2:
            raise ValueError("Radial impeller requires inlet_hub_diameter < d1 < d2; change duty, phi1, or diameters")
        self.dimensions = ImpellerDimensions(nq, psi, d1, dh, ds, d2, b2,
                                             pi*n*d1/60, pi*n*d2/60, qi, blade_count)
        if meridional is not None and axial_length is not None:
            raise ValueError("Specify meridional or axial_length, not both")
        mer = meridional if meridional is not None else meridional_contours(d1, dh, d2, b2, axial_length=axial_length)
        if not isinstance(mer, MeridionalContours):
            raise TypeError("meridional must be MeridionalContours")
        if not np.allclose([mer.hub[0,1], mer.shroud[0,1], mer.hub[-1,1], mer.shroud[-1,1],
                            mer.hub[-1,0]-mer.shroud[-1,0]], [dh/2, d1/2, d2/2, d2/2, b2], rtol=1e-6, atol=1e-10):
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
        cm2 = qi/(pi*d2*b2)
        inlet, outlet, sections, section_paths = [], [], [], []
        for span in (0.0, 0.5, 1.0):
            x, r = active.sample(settings.start_hub, span)
            ctx = dict(span_fraction=span, x=float(x), r=float(r))
            u1 = 2*pi*n/60*r
            cu1 = design_value(inlet_prewhirl, label="inlet_prewhirl", **ctx)
            flow = float(np.degrees(np.arctan2(cm1, u1-cu1)))
            metal = (flow + design_value(incidence_deg, label="incidence_deg", **ctx)
                     if beta1_deg is None else design_value(beta1_deg, label="beta1_deg", **ctx))
            x2, r2 = active.sample(1., span)
            beta2 = design_value(beta2_deg, label="beta2_deg", span_fraction=span, x=float(x2), r=float(r2))
            if not 0 < metal < 90 or not 0 < beta2 <= 90:
                raise ValueError("Supported backward-swept metal angles: 0 < beta1 < 90, 0 < beta2 <= 90 degrees")
            cu2 = self.dimensions.u2-cm2/np.tan(np.radians(beta2))
            inlet.append(VelocityTriangle(span, r, u1, cm1, cu1, np.hypot(cm1,u1-cu1), flow, metal))
            outlet.append(VelocityTriangle(span, r2, self.dimensions.u2, cm2, cu2,
                                            np.hypot(cm2,self.dimensions.u2-cu2), beta2, beta2))
            s = np.linspace(settings.start_hub, settings.end, 513)
            xr = active.sample(s, span)
            dm = np.linalg.norm(np.diff(xr, axis=0), axis=1)
            mp = np.r_[0., np.cumsum(dm/(0.5*(xr[:-1,1]+xr[1:,1])))]
            section_paths.append((mp/mp[-1], xr))
            f = mp/mp[-1]
            f = 0.5*(1-np.cos(pi*f)) if settings.angle_law == "cosine" else f
            beta = np.radians(metal+(beta2-metal)*f)
            wrap = (np.degrees(np.trapezoid(1/np.tan(beta), mp)) if wrap_deg is None
                    else design_value(wrap_deg, label="wrap_deg", **ctx))
            _positive(wrap, "wrap_deg")
            sections.append(BladeSection(span, metal, beta2, wrap))
        euler = float((outlet[1].u*outlet[1].cu-inlet[1].u*inlet[1].cu)/G0)
        self.velocity_triangles = ImpellerVelocityTriangles(tuple(inlet), tuple(outlet), euler,
                                                             None if hydraulic_efficiency is None else hydraulic_efficiency*euler)
        # Validate the callback contract before loading the optional backend.
        for section, (fraction, path) in zip(sections, section_paths, strict=True):
            for chord in np.linspace(0, 1, 17):
                xr = [np.interp(chord, fraction, path[:,k]) for k in range(2)]
                t = design_value(blade_thickness, label="blade_thickness", chord_fraction=float(chord),
                                 span_fraction=section.span_fraction, x=float(xr[0]), r=float(xr[1]))
                _positive(t, "blade_thickness")
        self.blade_definition = ImpellerBladeDefinition(mer, tuple(sections), blade_count,
                                                        blade_thickness, settings, splitter, tip_clearance)
        self.geometry = None

    @classmethod
    def from_pressure_rise(cls, *, mdot, rho, dp, n, **kwargs):
        for label, v in (("mdot",mdot), ("rho",rho), ("dp",dp)):
            _positive(v, label)
        return cls(Q=mdot/rho, H=dp/(rho*G0), n=n, rho=rho, **kwargs)

    @property
    def overall_efficiency(self):
        if self.hydraulic_efficiency is None or self.mechanical_efficiency is None:
            return None
        return self.volumetric_efficiency*self.hydraulic_efficiency*self.mechanical_efficiency

    @property
    def shaft_power(self):
        if self.rho is None or self.overall_efficiency is None:
            return None
        return self.rho*G0*self.Q*self.H/self.overall_efficiency

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
