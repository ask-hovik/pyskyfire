"""Translate physical pump inputs to pyTurbo 1.3.8 construction controls.

The conformal (m',theta) plane uses dm'=dm/r. Metal angles are from +m',
and upstream wrap_angle means atan(dtheta/m'), not total blade wrap.
Thickness controls are section-normal offsets in this plane: locally t/(2r).
"""
from importlib.metadata import version, PackageNotFoundError
import warnings

import numpy as np

from pyskyfire.common.rotation import rotation_sign

from .geometry import BladeSurface, ImpellerGeometry, cylindrical_to_cartesian
from .impeller import design_value


def build_geometry(definition, *, span_points, chord_points):
    for label, count, minimum in (("span_points", span_points, 3), ("chord_points", chord_points, 33)):
        if isinstance(count, bool) or not isinstance(count, (int, np.integer)) or count < minimum:
            raise ValueError(f"{label} must be an integer >= {minimum}")
    try:
        installed = version("pyturbo-aero")
    except PackageNotFoundError as error:
        raise ImportError("Blade generation requires pip install 'pyskyfire[impeller]'") from error
    if installed != "1.3.8":
        raise ImportError(f"The impeller adapter supports pyturbo-aero==1.3.8; found {installed}. "
                          "Install 'pyskyfire[impeller]'.")
    from pyturbo.aero import Centrif, CentrifProfile, TrailingEdgeProperties
    from pyturbo.helper import xr_to_mprime
    from shapely.affinity import translate
    from shapely.geometry import Polygon

    class IsolatedCentrif(Centrif):
        def __init__(self, **kwargs):
            # The upstream constructor appends to a class-level mutable list.
            # Shadow it before initialization so independent builds stay local.
            self.patterns = []
            super().__init__(**kwargs)

    settings = definition.settings
    passage = definition.meridional.with_clearance(definition.tip_clearance)
    # Upstream scales span parameters without moving profile xr coordinates.
    # Supply the actual blade-tip boundary to enforce physical tip clearance.
    model = IsolatedCentrif(blade_position=(settings.start_hub, settings.end), use_mid_wrap_angle=False)
    model.add_hub(passage.hub[:,0], passage.hub[:,1])
    model.add_shroud(passage.shroud[:,0], passage.shroud[:,1])
    model.tip_clearance = 0.0
    split_profiles, profile_maps = [], []
    for section in definition.sections:
        span = section.span_fraction
        xr = passage.sample(np.linspace(settings.start_hub, settings.end, chord_points), span)
        mp = xr_to_mprime(xr)[0]
        extent = mp[-1]
        profile_maps.append((mp/extent, xr))

        def thickness(t, scale=1.0):
            x = float(np.interp(t, mp/extent, xr[:,0]))
            r = float(np.interp(t, mp/extent, xr[:,1]))
            physical = scale*design_value(definition.thickness, label="blade_thickness",
                                         chord_fraction=float(t), span_fraction=span, x=x, r=r)
            if physical <= 0:
                raise ValueError("blade_thickness must be positive throughout the blade")
            return physical/(2*r)

        def profile(start=0., scale=1.0):
            ts = np.linspace(start+(1-start)*0.01, start+(1-start)*0.99, 33)
            # Upstream j-3 indexing repeats its last thickness entry at the
            # first body point. Preserve the TE and all intervening positions;
            # check the generated body thickness independently below.
            offsets = [thickness(t, scale) for t in ts[1:]]
            return CentrifProfile(
                percent_span=span, LE_Thickness=thickness(start, scale),
                LE_Metal_Angle=section.beta1_deg-90,
                TE_Metal_Angle=section.beta2_deg-90,
                # These locations make m'(t) linear in the camber parameter.
                LE_Metal_Angle_Loc=1/3, TE_Metal_Angle_Loc=2/3,
                wrap_angle=float(np.degrees(np.arctan(-np.radians(section.wrap_deg)/extent))),
                ss_thickness=offsets, ps_thickness=offsets,
                thickness_start=0.01, thickness_end=0.99, camber_follow_density=40,
                trailing_edge_properties=TrailingEdgeProperties(TE_Cut=settings.trailing_edge == "cut"),
            )
        model.add_profile(profile())
        if definition.splitter is not None:
            split_profiles.append(profile(definition.splitter.start, definition.splitter.thickness_scale))
    if split_profiles:
        model.add_splitter(split_profiles, splitter_start=definition.splitter.start)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Method 'bounded' does not support relative tolerance.*", category=RuntimeWarning)
        model.build(npts_span=span_points, npts_chord=chord_points, nblades=1,
                    nsplitters=int(bool(split_profiles)))
    # Upstream stacks every section on a common leading edge. Each span wraps by
    # its own angle law over its own conformal length, so that leaves the
    # trailing edge trailing across tens of degrees instead of lying in one
    # axial plane. Stack on the trailing edge instead: a constant theta offset
    # per span is a rigid rotation about the shaft, so metal angles, wrap,
    # thickness and the meridional passage are untouched, and the leading edge
    # carries the spanwise wrap difference as sweep.
    stack = model.mainblade.ss_cyl_pts[:,-1,2].copy()
    # pyTurbo builds in the counterclockwise sense. A clockwise wheel is the
    # mirror image, so negating theta after stacking reverses it without
    # touching the meridional passage, the metal angles or the thickness law.
    sign = rotation_sign(definition.rotation)

    def stacked(blade):
        sides = []
        for cyl in (blade.ss_cyl_pts, blade.ps_cyl_pts):
            shifted = np.array(cyl, dtype=float)
            shifted[:,:,2] = sign*(shifted[:,:,2] - stack[:,None])
            sides.append(shifted)
        return BladeSurface(*sides)

    # A splitter follows the main camber to the same trailing edge, so the same
    # offset keeps it in its passage.
    main = stacked(model.mainblade)
    splitter = stacked(model.splitterblade) if split_profiles else None
    errors = []
    cases = [(model.profiles_debug, 0., 1.)]
    if split_profiles:
        cases.append((model.splitter_debug, definition.splitter.start, definition.splitter.thickness_scale))
    for debug_list, start, scale in cases:
        for section, debug, (mfrac, xr) in zip(definition.sections, debug_list, profile_maps, strict=True):
            ss, ps = debug.ss_mp_pts, debug.ps_mp_pts
            frac = (ss[:,0]+ps[:,0])/(2*debug.camber_mp_th[-1,0])
            normal_gap = np.hypot(ss[:,0]-ps[:,0], ss[:,1]-ps[:,1])
            # The closure regions taper to zero by construction.
            for j in np.flatnonzero((frac > start+0.12*(1-start)) & (frac < start+0.88*(1-start))):
                x = float(np.interp(frac[j], mfrac, xr[:,0]))
                r = float(np.interp(frac[j], mfrac, xr[:,1]))
                target = scale*design_value(definition.thickness, label="blade_thickness",
                                            chord_fraction=float(frac[j]), span_fraction=section.span_fraction, x=x, r=r)
                if target <= 0:
                    raise ValueError("blade_thickness must be positive throughout the blade")
                errors.append(abs(normal_gap[j]*r/target-1))
    max_error = max(errors, default=0.)
    if max_error > 0.20:
        raise ValueError(f"pyTurbo cannot reproduce this thickness law within 20% in the blade body "
                         f"(maximum error {max_error:.1%}); use a smoother thickness law")
    for blade in (main, splitter):
        if blade is None:
            continue
        for ss, ps in zip(blade.suction, blade.pressure, strict=True):
            for side in (ss, ps):
                xyz = cylindrical_to_cartesian(side)
                if np.any(np.linalg.norm(np.diff(xyz, axis=0), axis=1) <= 1e-14):
                    raise ValueError("pyTurbo generated duplicate adjacent surface points")
    # Compare closed sections in the conformal plane. Rotation is translation
    # in theta, allowing a pitch check at corresponding physical span sections
    # even when the two sides have different meridional coordinates.
    pitch = 2*np.pi/definition.blade_count
    for i in range(span_points):
        def polygon(blade):
            ss, ps = blade.ss_mp_pts[i], blade.ps_mp_pts[i]
            shape = Polygon(np.concatenate((ss, ps[-2:0:-1])))
            if not shape.is_valid or shape.area <= 0:
                raise ValueError("pyTurbo generated a self-intersecting blade section")
            if shape.intersects(translate(shape, yoff=pitch)):
                raise ValueError("Blade sections overlap adjacent blades")
            return shape

        main_section = polygon(model.mainblade)
        if splitter is not None:
            split_section = translate(polygon(model.splitterblade),
                                      yoff=definition.splitter.pitch_fraction*pitch)
            if (main_section.intersects(split_section) or
                    translate(main_section, yoff=pitch).intersects(split_section)):
                raise ValueError("Splitter overlaps a main blade; change pitch_fraction or thickness")
    return ImpellerGeometry(definition.meridional.hub, definition.meridional.shroud,
                            main, definition.blade_count, splitter,
                            0.5 if definition.splitter is None else definition.splitter.pitch_fraction,
                            installed, max_error, definition.rotation)
