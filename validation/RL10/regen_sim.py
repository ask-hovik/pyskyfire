import json
import os

import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator
import pyskyfire as psf


script_dir = os.path.dirname(os.path.abspath(__file__))
reference_dir = os.path.join(script_dir, "reference_data")


# ====================
# Silver throat insert
# ====================

def silver_insert_thickness_from_json(path):
    """Load the RL10 silver insert and return its radial thickness function.

    The hot-gas-facing inner profile has a sharp minimum at ``x=0``, so it is
    represented by independent chamber- and nozzle-side cubic splines. The
    tube-facing outer profile is represented by one cubic spline. The returned
    callable evaluates ``r_outer - r_inner`` over the insert span and zero
    everywhere else.
    """
    with open(path, encoding="utf-8") as insert_file:
        data = json.load(insert_file)

    if data.get("units") != {"x": "m", "r": "m"}:
        raise ValueError("silver-insert coordinates must be stored in metres")

    def curve_coordinates(name):
        try:
            curve = data["curves"][name]
            x_curve = np.asarray(curve["x"], dtype=float)
            r_curve = np.asarray(curve["r"], dtype=float)
        except (KeyError, TypeError) as exc:
            raise ValueError(f"silver insert is missing the {name} curve") from exc
        if x_curve.ndim != 1 or r_curve.ndim != 1:
            raise ValueError(f"silver-insert {name} coordinates must be 1-D")
        if x_curve.size < 3 or x_curve.size != r_curve.size:
            raise ValueError(
                f"silver-insert {name} coordinates must have equal length >= 3"
            )
        if np.any(~np.isfinite(x_curve)) or np.any(~np.isfinite(r_curve)):
            raise ValueError(f"silver-insert {name} curve must be finite")
        if np.any(np.diff(x_curve) <= 0.0):
            raise ValueError(
                f"silver-insert {name} axial coordinates must increase strictly"
            )
        if np.any(r_curve <= 0.0):
            raise ValueError(f"silver-insert {name} radii must be positive")
        return x_curve, r_curve

    inner_x, inner_r = curve_coordinates("inner")
    outer_x, outer_r = curve_coordinates("outer")
    throat_indices = np.flatnonzero(inner_x == 0.0)
    if throat_indices.size != 1:
        raise ValueError("silver-insert inner curve must contain one x=0 throat")
    throat_index = int(throat_indices[0])
    if throat_index == 0 or throat_index == inner_x.size - 1:
        raise ValueError("silver-insert throat must have points on both sides")
    if throat_index != int(np.argmin(inner_r)):
        raise ValueError("silver-insert x=0 point must be its minimum radius")
    if not np.allclose(
        [inner_x[0], inner_r[0], inner_x[-1], inner_r[-1]],
        [outer_x[0], outer_r[0], outer_x[-1], outer_r[-1]],
    ):
        raise ValueError("silver-insert inner and outer endpoints must match")

    inner_chamber_spline = CubicSpline(
        inner_x[: throat_index + 1],
        inner_r[: throat_index + 1],
        extrapolate=False,
    )
    inner_nozzle_spline = CubicSpline(
        inner_x[throat_index:],
        inner_r[throat_index:],
        extrapolate=False,
    )
    outer_spline = CubicSpline(outer_x, outer_r, extrapolate=False)
    x_start = float(inner_x[0])
    x_end = float(inner_x[-1])

    def thickness(x):
        values = np.asarray(x, dtype=float)
        result = np.zeros_like(values, dtype=float)
        inside = (values >= x_start) & (values <= x_end)
        chamber = inside & (values <= 0.0)
        nozzle = inside & (values > 0.0)
        inner_values = np.empty_like(values, dtype=float)
        inner_values[chamber] = inner_chamber_spline(values[chamber])
        inner_values[nozzle] = inner_nozzle_spline(values[nozzle])
        result[inside] = np.maximum(
            outer_spline(values[inside]) - inner_values[inside],
            0.0,
        )
        return float(result) if result.ndim == 0 else result

    sample_x = np.linspace(x_start, x_end, 1001)
    sample_inner_r = np.empty_like(sample_x)
    chamber_sample = sample_x <= 0.0
    sample_inner_r[chamber_sample] = inner_chamber_spline(
        sample_x[chamber_sample]
    )
    sample_inner_r[~chamber_sample] = inner_nozzle_spline(
        sample_x[~chamber_sample]
    )
    if np.any(outer_spline(sample_x) - sample_inner_r < -1.0e-12):
        raise ValueError("silver-insert outer curve crosses its inner curve")

    thickness.throat_radius = float(inner_r[throat_index])
    thickness.inner_x = inner_x
    thickness.inner_r = inner_r
    thickness.outer_x = outer_x
    thickness.outer_r = outer_r
    thickness.x_span = (x_start, x_end)
    thickness.source = data.get("source")
    return thickness


# ==================
# Contour from Table
# ==================

def contour_from_area_ratios(
    path,
    r_t,
    exit_x,
    exit_area_ratio,
):
    """Build a contour directly from the paper's tabulated area ratios.

    Binder et al. Table E1 lists the hot-gas property node locations as area
    ratios ``A/A_t`` rather than radii, so the throat radius has to be supplied.
    The physical throat radius is supplied from the detailed silver-insert
    profile. Do not use ``aerothermodynamics.r_t``; that is derived from the
    ideal ``c_star`` and lands near P&W's effective flow area instead.

    Parameters
    ----------
    path : str
        Path to ``reference_area_ratios.json``.
    r_t : float
        Throat radius [m]. Every station is placed at ``r_t * sqrt(A/A_t)``.
    exit_x : float
        Nozzle exit-plane location relative to the throat [m].
    exit_area_ratio : float
        Nozzle exit area ratio ``A_e/A_t``.
    Returns
    -------
    psf.regen.Contour
        Contour whose input points are the published stations plus the physical
        insert throat and nozzle exit.

    Notes
    -----
    The table has no rows at the throat or exit, so ``(0, 1)`` and
    ``(exit_x, exit_area_ratio)`` are inserted before converting each ratio to
    a radius. The detailed insert curve supplies the throat radius only; its
    other points are deliberately not added to the contour. No intermediate
    grid is created.
    """
    with open(path, encoding="utf-8") as f:
        data = json.load(f)

    stations = sorted(data["stations"], key=lambda s: s["x"])
    x_station = np.array([s["x"] for s in stations], dtype=float)
    area_ratio = np.array(
        [s["area_ratio"] for s in stations],
        dtype=float,
    )

    if r_t <= 0.0 or not np.isfinite(r_t):
        raise ValueError("r_t must be finite and positive")
    if exit_x <= x_station[-1] or not np.isfinite(exit_x):
        raise ValueError("exit_x must be finite and beyond the last table station")
    if exit_area_ratio <= 0.0 or not np.isfinite(exit_area_ratio):
        raise ValueError("exit_area_ratio must be finite and positive")
    if np.any(~np.isfinite(x_station)) or np.any(~np.isfinite(area_ratio)):
        raise ValueError("area-ratio table contains non-finite values")
    if np.any(np.diff(x_station) <= 0.0):
        raise ValueError("area-ratio x stations must increase strictly")
    if np.any(area_ratio <= 0.0):
        raise ValueError("area ratios must be positive")
    if np.any(x_station == 0.0):
        raise ValueError("area-ratio table already contains a throat station")

    # Binder et al. Table E1 omits both of these physical boundary stations.
    x_station = np.concatenate([x_station, [0.0, exit_x]])
    area_ratio = np.concatenate([area_ratio, [1.0, exit_area_ratio]])
    order = np.argsort(x_station)
    x_station = x_station[order]
    area_ratio = area_ratio[order]
    r_station = r_t * np.sqrt(area_ratio)

    return psf.regen.Contour(
        x_station,
        r_station,
        name="Area-ratio stations",
    )


# =============
# Engine inputs
# =============

LBM = 0.45359237  # kg
PSI = 6894.757    # Pa

silver_insert_path = os.path.join(reference_dir, "reference_silver_insert.json")
silver_insert_thickness = silver_insert_thickness_from_json(silver_insert_path)

params = dict(
    # All values from Binder et. al. unless stated otherwise.
    #
    # Operating point: Tables 2.4.1 and 2.5.1 share the footnote "values taken at
    # typical engine operating point as predicted by the model", and 2.5.1's
    # implied fuel flow (37.36/(1+5.26) = 5.968 lbm/s) matches 2.4.1's stated
    # cooling-jacket flow (5.973 lbm/s) to 0.08%. They are one operating point,
    # and it is the one the jacket reference data in 2.4.1 is quoted at.
    p_c = 482.0 * PSI,              # Table 2.5.1, injector-face static
    mdot_fu = 5.973 * LBM,          # Table 2.4.1, cooling-jacket mass flow
    mdot_ox = (37.36 - 5.973) * LBM,  # Table 2.5.1 total (37.36) minus the fuel
    # -> MR = 5.255, against the 5.26 stated in Table 2.5.1.
    #
    # The alternative is the engine station data, which is a different balance:
    # fuel 2.7587, ox 13.9480 kg/s at the injector (MR 5.056). Note its MR at the
    # *pump inlets* is 4.993, the 1.28% difference being fuel lost overboard to
    # the shaft-seal / gear-box leakage that cools the bearings (sec. 2.2).
    eps = 61.0,   # Table 2.5.1, with the silver throat insert (57.0 without)
    # Figure E1 places the nozzle exit at approximately 46 in downstream of
    # the throat in the same coordinate system used by Table E1.
    x_e = 46.0 * 0.0254,

    # Fuel/oxidizer parameters
    cea_fu = psf.common.Fluid(type="fuel", propellants=["H2"], fractions=[1.0]),
    cea_ox = psf.common.Fluid(type="oxidizer", propellants=["O2(L)"], fractions=[1.0]),
    coolant_fu =  psf.common.Fluid(type="fuel", propellants=["Hydrogen"], fractions=[1.0]),
    T_gas_fu_in = 200,  # Adjusted by a few kelvin to fall within NASA CEA tables
    T_gas_ox_in = 100,  # Also adjusted
    minimum_cea_temperature = 200.0,
    T_coolant_in = 33,
    p_coolant_in = 69e5, 

    # Chamber-nozzle parameters, contour from Binder et. al.
    #theta_conv = 25,
    #R_1f = 1.5,
    #R_2f = 3,
    #R_3f = 0.5,
    #length_fraction = 0.713,
    L_star = 0.95,  # I can't remember where I read this (not Binder)

    # Cooling channel
    roughness_height = 1.1684e-6,  # found in Binder et al 4.6e-5 inches. Conventional SS tubing is often cited as 0.0015e-3
    wall_thickness = 0.3302e-3,
)


# =====
# Setup
# =====

aerothermodynamics = psf.skycea.Aerothermodynamics.from_mass_flows_eps_Lstar(fu=params["cea_fu"],
                                                                    ox=params["cea_ox"],
                                                                    T_fu_in=params["T_gas_fu_in"],
                                                                    T_ox_in=params["T_gas_ox_in"],
                                                                    mdot_fu=params["mdot_fu"],
                                                                    mdot_ox=params["mdot_ox"],
                                                                    p_c=params["p_c"],
                                                                    eps=params["eps"],
                                                                    L_star=params["L_star"],
                                                                    p_amb=1e5, # ca atmospheric pressure, optional
                                                                    transport="frozen",
                                                                    minimum_cea_temperature=params["minimum_cea_temperature"],
                                                                    )

def _reference_model_x(filename, model="new_system_model"):
    """Load the digitised axial stations for one RL10 reference curve."""
    path = os.path.join(reference_dir, filename)
    with open(path, encoding="utf-8") as reference_file:
        return np.asarray(json.load(reference_file)[model]["x"], dtype=float)


# Use the RTE wall-temperature stations for both circuits. RTE's outermost
# nozzle station lies beyond the modeled cooling-jacket span and is removed
# below once the physical pass boundaries have been established.
wall_reference_x = _reference_model_x(
    "reference_wall_temperature.json",
    model="RTE_model",
)
heat_flux_reference_x = _reference_model_x("reference_heat_flux.json")
system_coolant_reference_x = _reference_model_x(
    "reference_coolant_static_temperature.json"
)
pressure_reference_x = _reference_model_x(
    "reference_coolant_static_pressure.json"
)
rte_pressure_reference_x = _reference_model_x(
    "reference_coolant_static_pressure.json",
    model="RTE_model",
)

# The first New System Model coolant-temperature point is the 360-to-180 tube
# interleaving transition. Retain its full digitized precision. Keep this as an
# explicit geometry setting because the low-resolution RECOP width profile may
# justify moving the transition after the two datasets have been compared.
reference_interlacing_x = float(system_coolant_reference_x[0])
interlacing_x = reference_interlacing_x
half_coolant_x = np.array(
    [interlacing_x, system_coolant_reference_x[1]],
    dtype=float,
)

# RTE treated the short-tube half pass differently from the physical RL10
# circuit, so use the reference interlacing location and the New System Model's
# physical turnaround endpoint for that pass.
# On the long-tube full pass, use the distributed RTE pressure stations while
# retaining the New System Model's turnaround and injector-end boundaries.
full_pass_start = system_coolant_reference_x[1]
full_pass_end = system_coolant_reference_x[-1]
wall_reference_x = wall_reference_x[
    (wall_reference_x <= full_pass_start)
    & (wall_reference_x >= full_pass_end)
]
full_coolant_x = np.unique(
    np.concatenate(
        (
            [full_pass_start],
            rte_pressure_reference_x[
                (rte_pressure_reference_x < full_pass_start)
                & (rte_pressure_reference_x > full_pass_end)
            ],
            [full_pass_end],
        )
    )
)[::-1]


def _wall_nodes_for_coolant_segments(reference_x, coolant_x):
    """Add a wall solve only where an RTE coolant interval has none."""
    wall_x = list(np.asarray(reference_x, dtype=float))
    direction = -1.0 if coolant_x[0] > coolant_x[-1] else 1.0
    s_coolant = direction * np.asarray(coolant_x, dtype=float)
    s_wall = direction * np.asarray(wall_x, dtype=float)
    owners = np.clip(
        np.searchsorted(s_coolant, s_wall, side="right") - 1,
        0,
        len(coolant_x) - 2,
    )
    for segment in range(len(coolant_x) - 1):
        if not np.any(owners == segment):
            # Bias the support point toward the segment outlet. The broad
            # nozzle intervals otherwise place a midpoint near x=0.48 m,
            # where the RL10 CEA wall balance is poorly conditioned.
            wall_x.append(
                0.10 * coolant_x[segment]
                + 0.90 * coolant_x[segment + 1]
            )
    return np.unique(wall_x)


full_wall_x = _wall_nodes_for_coolant_segments(
    wall_reference_x,
    full_coolant_x,
)


def _thermal_nodes_in_span(reference_x, coolant_x):
    lower, upper = np.min(coolant_x), np.max(coolant_x)
    tolerance = 1.0e-3 * (upper - lower)
    return reference_x[
        (reference_x >= lower - tolerance)
        & (reference_x <= upper + tolerance)
    ]


half_wall_x = _thermal_nodes_in_span(wall_reference_x, half_coolant_x)
# RTE did not represent the physical short-tube half pass. Its first station
# above the interlacing manifold (x=0.3029 m) is therefore retained on the
# long-tube full pass but not projected onto the short tubes; doing so places a
# wall balance directly in the geometric transition and is poorly conditioned.
half_wall_x = half_wall_x[half_wall_x > np.min(half_wall_x)]

half_pass_nodes = {
    "wall": half_wall_x,
    "heat_flux": _thermal_nodes_in_span(
        heat_flux_reference_x, half_coolant_x
    ),
    "coolant": half_coolant_x,
}
full_pass_nodes = {
    "wall": full_wall_x,
    "heat_flux": heat_flux_reference_x,
    "coolant": full_coolant_x,
}

area_ratio_path = os.path.join(reference_dir, "reference_area_ratios.json")
contour = contour_from_area_ratios(
    area_ratio_path,
    r_t=silver_insert_thickness.throat_radius,
    exit_x=params["x_e"],
    exit_area_ratio=params["eps"],
)

wall = psf.regen.Wall(material = psf.common.solids.StainlessSteel347, thickness = params["wall_thickness"]) 
silver_wall = psf.regen.Wall(
    material=psf.common.solids.PureSilver,
    thickness=silver_insert_thickness,
    name="Silver throat insert",
)
#wall = psf.regen.Wall(material = psf.common.solids.copper, thickness = params["wall_thickness"] )

analytic_channel_height_fn = psf.regen.make_channel_height_fn(
    contour=contour, 
    region_fractions=[-1.0, 0.25, 1.0], 
    flat_heights= [0.0032, 0.00134], 
    pinch_factors= [0.8, -5.0], #baseline had [0.6, -5.0]
    transition_widths=[0.01]
) # total coolant volume should be ca 0.015831543m3

channel_height_path = os.path.join(
    reference_dir,
    "reference_channel_height.json",
)

# Approximate transition in the RECOP Figure 10 coordinate system. The profile
# loader maps this point to ``interlacing_x`` while holding the throat and exit
# fixed. Both values remain explicit so the transition can be revisited without
# editing the source data.
recop_profile_interlacing_x = 0.3336


def align_positive_channel_height_x(
    x,
    source_interleaving_x,
    target_interleaving_x,
    *,
    source_end_x=None,
    target_end_x=None,
):
    """Align a positive-side channel profile to a new interleaving point.

    Coordinates from the throat through the source interleaving point are
    scaled with the throat fixed at zero. Downstream coordinates retain their
    original nozzle spacing unless source and target endpoints are supplied.
    With endpoints, the downstream interval is scaled so the transition and
    nozzle exit are both fixed. Negative chamber coordinates are unchanged.
    """
    values = np.asarray(x, dtype=float)
    source_interleaving_x = float(source_interleaving_x)
    target_interleaving_x = float(target_interleaving_x)
    if not np.all(np.isfinite(values)):
        raise ValueError("channel-height coordinates must be finite")
    if not np.isfinite(source_interleaving_x) or source_interleaving_x <= 0.0:
        raise ValueError("source interleaving position must be finite and positive")
    if not np.isfinite(target_interleaving_x) or target_interleaving_x <= 0.0:
        raise ValueError("target interleaving position must be finite and positive")
    if (source_end_x is None) != (target_end_x is None):
        raise ValueError("source_end_x and target_end_x must be supplied together")
    if source_end_x is not None:
        source_end_x = float(source_end_x)
        target_end_x = float(target_end_x)
        if not np.isfinite(source_end_x) or source_end_x <= source_interleaving_x:
            raise ValueError("source_end_x must be beyond the source interleaving point")
        if not np.isfinite(target_end_x) or target_end_x <= target_interleaving_x:
            raise ValueError("target_end_x must be beyond the target interleaving point")

    aligned = values.copy()
    before_interleaving = (values > 0.0) & (
        values <= source_interleaving_x
    )
    aligned[before_interleaving] = (
        values[before_interleaving]
        * target_interleaving_x
        / source_interleaving_x
    )
    downstream = values > source_interleaving_x
    downstream_scale = 1.0
    if source_end_x is not None:
        downstream_scale = (
            (target_end_x - target_interleaving_x)
            / (source_end_x - source_interleaving_x)
        )
    aligned[downstream] = target_interleaving_x + downstream_scale * (
        values[downstream] - source_interleaving_x
    )
    return float(aligned) if aligned.ndim == 0 else aligned


def channel_height_from_json(
    path,
    target_interleaving_x,
    *,
    source_interleaving_x,
    wall_thickness,
):
    """Load RECOP tube OD height and return internal coolant height.

    RECOP Figure 10 reports ``TUBE HEIGHT (OD)``. The rounded cross-section
    input is the full internal coolant height, so twice the tube wall thickness
    is removed before constructing the PCHIP profile.
    """
    with open(path, encoding="utf-8") as profile_file:
        profile = json.load(profile_file)

    try:
        source_series = profile["series"]["design_points"]
    except (KeyError, TypeError) as exc:
        raise ValueError("channel-height reference is missing design_points") from exc

    original_x_profile = np.asarray(source_series["x"], dtype=float)
    x_profile = align_positive_channel_height_x(
        original_x_profile,
        source_interleaving_x,
        target_interleaving_x,
        source_end_x=original_x_profile[-1],
        target_end_x=contour.xs[-1],
    )
    tube_height_od = np.asarray(source_series["tube_height_od"], dtype=float)
    if callable(wall_thickness):
        wall_profile = np.asarray(wall_thickness(x_profile), dtype=float)
    else:
        wall_profile = np.full_like(x_profile, float(wall_thickness))
    height_profile = tube_height_od - 2.0 * wall_profile

    if x_profile.ndim != 1 or height_profile.ndim != 1:
        raise ValueError("channel-height arrays must be one-dimensional")
    if x_profile.size < 2 or x_profile.size != height_profile.size:
        raise ValueError("channel-height arrays must have equal length >= 2")
    if wall_profile.shape != x_profile.shape:
        raise ValueError("wall-thickness profile must match channel-height x")
    if not np.all(np.isfinite(x_profile)) or not np.all(np.isfinite(height_profile)):
        raise ValueError("channel-height profile contains non-finite values")
    if not np.all(np.diff(x_profile) > 0.0):
        raise ValueError("channel-height x_m values must increase strictly")
    if np.any(height_profile <= 0.0):
        raise ValueError("channel heights must be positive")
    if x_profile[0] > contour.xs[0] or x_profile[-1] < contour.xs[-1]:
        raise ValueError("channel-height profile must cover the contour")

    interpolator = PchipInterpolator(
        x_profile,
        height_profile,
        extrapolate=False,
    )

    def channel_height(x):
        values = interpolator(x)
        if np.any(~np.isfinite(values)):
            raise ValueError(f"channel height requested outside profile at x={x}")
        return values

    aligned_profile = dict(profile)
    aligned_profile["original_x_m"] = original_x_profile.tolist()
    aligned_profile["x_m"] = x_profile.tolist()
    aligned_profile["tube_height_od_m"] = tube_height_od.tolist()
    aligned_profile["wall_thickness_m"] = wall_profile.tolist()
    aligned_profile["channel_height_internal_m"] = height_profile.tolist()
    aligned_profile["source_interleaving_x_m"] = float(source_interleaving_x)
    aligned_profile["target_interleaving_x_m"] = float(
        target_interleaving_x
    )
    channel_height.profile = aligned_profile
    channel_height.x_profile = x_profile
    channel_height.source_path = path
    return channel_height


channel_height_fn = channel_height_from_json(
    channel_height_path,
    target_interleaving_x=interlacing_x,
    source_interleaving_x=recop_profile_interlacing_x,
    wall_thickness=params["wall_thickness"],
)

def helix_fn(x):
    return 0#45*3.14*180

cross_section = psf.regen.CrossSectionRounded()
#cross_section_squared = psf.regen.CrossSectionSquared(blockage_ratio=0.1)
LH2_transport = psf.skycea.CoolantTransport(params["coolant_fu"])

def build_thrust_chamber(channel_height=channel_height_fn, *, compute_gas=False):
    """Construct the RL10 chamber for a supplied channel-height function."""
    half_pass = psf.regen.CoolingCircuit(
        name="Half Pass",
        contour=contour,
        coolant_transport=LH2_transport,
        cross_section=cross_section,
        span=[
            half_coolant_x[0] / contour.xs[-1],
            1.0, #half_coolant_x[-1] / contour.xs[-1],
        ],
        placement=psf.regen.SurfacePlacement(
            n_channel_positions=180,
            helix_angle=helix_fn,
        ),
        walls=[wall],
        roughness=params["roughness_height"],
        channel_height=channel_height,
    )
    
    full_pass = psf.regen.CoolingCircuit(
        name="Full Pass",
        contour=contour,
        coolant_transport=LH2_transport,
        cross_section=cross_section,
        span=[
            1.0, #full_coolant_x[0] / contour.xs[-1],
            -1.0,
        ],
        placement=psf.regen.SurfacePlacement(
            n_channel_positions=180,
            helix_angle=helix_fn,
        ),
        walls=[silver_wall, wall],
        roughness=params["roughness_height"],
        channel_height=channel_height,
    )
    
    return psf.regen.ThrustChamber(
        contour=contour,
        combustion_transport=aerothermodynamics,
        cooling_circuits=[half_pass, full_pass],
        h_gas_corr=1.08,
        h_cold_corr=1.0,
        n_nodes=150,
        enable_fin=True,
        compute_gas=compute_gas,
    )


# Populate the gas-property tables once. Optimization candidates reuse the
# same contour and combustion state and therefore set ``compute_gas=False``.
thrust_chamber = build_thrust_chamber(compute_gas=True)


def run_simulation(
    chamber=None,
    *,
    half_nodes=30,
    full_nodes=100,
    use_reference_nodes=False,
    output=True,
):
    """Run both RL10 cooling passes and return ``(chamber, [half, full])``.

    Set ``use_reference_nodes=True`` to reproduce the distinct wall,
    heat-flux, and coolant station locations digitised from Binder et al. In
    that mode, ``half_nodes`` and ``full_nodes`` are ignored. This is the
    default; set ``use_reference_nodes=False`` to use the uniform grids given
    by ``half_nodes`` and ``full_nodes`` instead.
    """
    chamber = thrust_chamber if chamber is None else chamber
    if use_reference_nodes:
        half_nodes = half_pass_nodes
        full_nodes = full_pass_nodes

    mdot_fu = aerothermodynamics.mdot_fu
    boundary_conditions_a = psf.regen.BoundaryConditions(
        T_coolant_in=params["T_coolant_in"],
        p_coolant_in=params["p_coolant_in"],
        mdot_coolant=mdot_fu,
    )
    cooling_data_a = psf.regen.coupled_steady_heating_analysis(
        chamber,
        nodes=50,#half_nodes,
        circuit_index=0,
        boundary_conditions=boundary_conditions_a,
        solver="newton",
        output=output,
        heat_curvature_correction=True,
    )

    boundary_conditions_b = psf.regen.BoundaryConditions(
        T_coolant_in=cooling_data_a.T_stagnation[-1],
        p_coolant_in=cooling_data_a.p_stagnation[-1],
        mdot_coolant=mdot_fu,
    )
    cooling_data_b = psf.regen.coupled_steady_heating_analysis(
        chamber,
        nodes=150,#full_nodes,
        circuit_index=1,
        boundary_conditions=boundary_conditions_b,
        solver="newton",
        output=output,
        heat_curvature_correction=True,
    )
    return chamber, [cooling_data_a, cooling_data_b]


def save_simulation(chamber, cooling_data, path=None):
    """Save a completed RL10 cooling simulation."""
    path = os.path.join(script_dir, "regen_results.pkl") if path is None else path
    results = psf.common.Results()
    results.add(name="params", obj=params)
    results.add(name="thrust_chamber", obj=chamber)
    results.add(name="cooling_data", obj=cooling_data)
    results.save(path)
    return path


def main(*, use_reference_nodes=True):
    chamber, cooling_data = run_simulation(
        use_reference_nodes=use_reference_nodes,
    )
    save_simulation(chamber, cooling_data)
    return chamber, cooling_data


if __name__ == "__main__":
    # Use ``main(use_reference_nodes=False)`` for the uniform 30/100-node run.
    main()
