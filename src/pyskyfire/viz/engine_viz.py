import time
from math import gcd
from functools import reduce
from pathlib import Path
import base64
import tempfile
import os

import numpy as np
import pyvista as pv

from pyskyfire.regen.cross_section import (
    CrossSectionSquared,
    CrossSectionRounded,
    CrossSectionRoundedInternal,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_EPS = 1e-12
TWO_PI = 2.0 * np.pi


# ---------------------------------------------------------------------------
# Small utilities
# ---------------------------------------------------------------------------

def gcd_list(ints):
    """Greatest common divisor of a sequence of integers."""
    return reduce(gcd, ints) if len(ints) else 1


def wrap_pi(delta: float) -> float:
    """Wrap an angle to (-pi, pi]."""
    return (delta + np.pi) % (2.0 * np.pi) - np.pi


def wedge_info(thrust_chamber):
    """Return (g, wedge_angle, n_list) for the chamber's cooling circuits."""
    n_list = [
        int(c.placement.n_channel_positions)
        for c in thrust_chamber.cooling_circuits
    ]
    g = gcd_list(n_list)
    wedge_angle = 2.0 * np.pi / g
    return g, wedge_angle, n_list


def circular_pattern(centerline: np.ndarray, n: int):
    """Replicate a single centerline n times around the circumference."""
    thetas = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    out = []
    for th in thetas:
        cl = centerline.copy()
        cl[:, 2] = cl[:, 2] + th
        out.append(cl)
    return out


# ---------------------------------------------------------------------------
# Channel collection / tiling
# ---------------------------------------------------------------------------

def collect_channels(thrust_chamber, *, jitter_divisor: float = 100.0):
    """
    Collect centerlines and widths for a minimal tiling wedge (if gcd>1)
    or the full circumference, plus circuit ids and meta info.

    Parameters
    ----------
    thrust_chamber
        Thrust chamber object with `cooling_circuits`.
    jitter_divisor : float
        Controls tiny angular jitter used to break symmetry.

    Returns
    -------
    centerlines : list[np.ndarray]
        List of (Ni,3) [x, r, theta] centerlines within one wedge.
    widths : list[np.ndarray]
        List of (Ni,) arrays for channel angular width vs x.
    circuit_ids : list[int]
        Circuit index for each centerline.
    meta : dict
        Contains "g", "wedge_angle", "n_list" (channel counts per circuit).
    """
    g, wedge_angle, n_list = wedge_info(thrust_chamber)

    def _jit_step(w_base: np.ndarray) -> float:
        w_base = np.asarray(w_base, float)
        w_mean = np.nanmean(w_base[w_base > 0]) if np.any(w_base > 0) else 0.0
        return (w_mean / jitter_divisor) if w_mean > 0.0 else 1e-6

    cl_out, w_out, cid_out = [], [], []
    k_global = 0

    for icirc, circuit in enumerate(thrust_chamber.cooling_circuits):
        base = np.asarray(circuit.centerlines[0], float)  # (Ni,3)
        n = int(circuit.placement.n_channel_positions)
        w_base = np.asarray(circuit.channel_width, float)
        step = _jit_step(w_base)

        # number from this circuit to include in one wedge
        m = n // g if g > 1 else n

        pats = circular_pattern(base, n)
        for k in range(m):
            cl = pats[k].copy()
            cl[:, 2] += k_global * step  # tiny asymmetry to avoid random “pop”
            cl_out.append(cl)
            w_out.append(w_base.copy())
            cid_out.append(icirc)
            k_global += 1

    meta = {"g": g, "wedge_angle": wedge_angle, "n_list": n_list}
    return cl_out, w_out, cid_out, meta


def tile_centerlines(
    centerlines_subset, widths_subset, circuit_ids_subset, g: int, wedge_angle: float
):
    """Replicate a wedge subset g times around the circumference."""
    full_cls, full_ws, full_cids = [], [], []
    for rot in range(g):
        dth = rot * wedge_angle
        for cl, w, cid in zip(centerlines_subset, widths_subset, circuit_ids_subset):
            cl2 = cl.copy()
            cl2[:, 2] = cl2[:, 2] + dth
            full_cls.append(cl2)
            full_ws.append(w.copy())
            full_cids.append(cid)
    return full_cls, full_ws, full_cids


# ---------------------------------------------------------------------------
# Redistribution (packing solver)
# ---------------------------------------------------------------------------


def redistribute_channels(
    centerlines,
    widths,
    *,
    k_pack: float = 1.0,
    k_bend: float = 1.0,
    k_slope: float = 3.0,          # NEW: soft anchoring strength
    desired_slope=None,            # NEW: list of arrays (one per centerline), same length as theta
    dt: float = 0.25,
    tol: float = 1e-3,
    max_iters: int = 800,
    gamma_damp: float = 0.4,
    wedge_angle: float | None = None,
    recorder: list | None = None,
    # NEW: per-centerline endpoint slopes dθ/dx [rad/m]
    endpoint_slope0=None,   # list[float] | None
    endpoint_slopeN=None,   # list[float] | None
):
    """
    Relax channel angles θ(x) so channels pack nicely without overlap.

    If endpoint_slope0/endpoint_slopeN are provided, they are enforced as:
      (θ[1] - θ[0]) / (x[1] - x[0]) = endpoint_slope0[i]
      (θ[-1] - θ[-2]) / (x[-1] - x[-2]) = endpoint_slopeN[i]
    """
    M = len(centerlines)
    xs_list = [cl[:, 0] for cl in centerlines]
    rs_list = [cl[:, 1] for cl in centerlines]
    th_list = [cl[:, 2].copy() for cl in centerlines]
    W_list = [w.copy() for w in widths]
    v_list = [np.zeros_like(th) for th in th_list]

    # Normalize slope inputs
    if endpoint_slope0 is None:
        endpoint_slope0 = [None] * M
    if endpoint_slopeN is None:
        endpoint_slopeN = [None] * M
    if len(endpoint_slope0) != M or len(endpoint_slopeN) != M:
        raise ValueError("endpoint_slope0/endpoint_slopeN must have length == number of centerlines")

    # Normalize distributed desired slope input
    if desired_slope is None:
        desired_slope = [None] * M
    if len(desired_slope) != M:
        raise ValueError("desired_slope must have length == number of centerlines")

    def _record():
        if recorder is not None:
            snapshot = [
                np.column_stack([xs_list[i], rs_list[i], th_list[i]]).astype(float, copy=True)
                for i in range(M)
            ]
            recorder.append(snapshot)

    def _apply_hard_endpoint_slopes():
        """Hard-set theta endpoints to match desired endpoint slopes."""
        for i in range(M):
            s0 = endpoint_slope0[i]
            sN = endpoint_slopeN[i]
            if s0 is None and sN is None:
                continue

            th = th_list[i]
            x = xs_list[i]
            if th.size < 2:
                continue

            dx0 = float(x[1] - x[0])
            dxN = float(x[-1] - x[-2])
            if abs(dx0) < 1e-14 or abs(dxN) < 1e-14:
                continue

            if s0 is not None:
                # (th[1]-th[0])/dx0 = s0  -> th[0] = th[1] - s0*dx0
                th[0] = th[1] - float(s0) * dx0

            if sN is not None:
                # (th[-1]-th[-2])/dxN = sN -> th[-1] = th[-2] + sN*dxN
                th[-1] = th[-2] + float(sN) * dxN

            th_list[i] = th

    x_master = np.unique(np.concatenate(xs_list))
    K = x_master.size

    active_masks = []
    for i in range(M):
        xmin, xmax = xs_list[i][0], xs_list[i][-1]
        if xmin > xmax:
            xmin, xmax = xmax, xmin
        active_masks.append((x_master >= xmin) & (x_master <= xmax))

    _record()

    # Ensure initial condition respects your endpoint intent
    _apply_hard_endpoint_slopes()

    for _ in range(int(max_iters)):
        prev_concat = np.concatenate(th_list)

        # θ, W on master grid (inactive -> 0)
        theta_m = np.zeros((M, K), float)
        width_m = np.zeros((M, K), float)
        for i in range(M):
            xmin, xmax = xs_list[i][0], xs_list[i][-1]
            if xmin > xmax:
                xmin, xmax = xmax, xmin

            th_i = np.interp(x_master, xs_list[i], th_list[i],
                             left=th_list[i][0], right=th_list[i][-1])
            w_i = np.interp(x_master, xs_list[i], W_list[i],
                            left=W_list[i][0], right=W_list[i][-1])

            mask = (x_master >= xmin) & (x_master <= xmax)
            theta_m[i, mask] = th_i[mask]
            width_m[i, mask] = w_i[mask]

        # Packing forces slice-by-slice (unchanged)
        Fm = np.zeros_like(theta_m)
        for j in range(K):
            active = [i for i in range(M) if active_masks[i][j]]
            if len(active) < 2:
                continue

            if wedge_angle is None:
                thj = np.array([theta_m[i, j] for i in active])
                order = np.argsort(thj)
                idx = [active[k] for k in order]
                for k in range(len(idx)):
                    a = idx[k]
                    b = idx[(k + 1) % len(idx)]
                    err = ((theta_m[a, j] + 0.5 * width_m[a, j])
                         - (theta_m[b, j] - 0.5 * width_m[b, j]))
                    err = (err + np.pi) % (2.0 * np.pi) - np.pi
                    Fm[a, j] += -k_pack * err
                    Fm[b, j] += +k_pack * err
            else:
                delta = float(wedge_angle)
                #print(delta)
                th_real = np.array([theta_m[i, j] for i in active])
                #print(th_real)
                w_real  = np.array([width_m[i, j] for i in active])
                ph = np.mod(th_real, delta)
                #input()

                ph_all = np.concatenate([ph - delta, ph, ph + delta])
                w_all  = np.concatenate([w_real, w_real, w_real])
                id_all = np.concatenate([active, active, active])

                order = np.argsort(ph_all)
                mid_mask = (ph_all >= 0.0) & (ph_all < delta)
                mid_positions = [p for p in order if mid_mask[p]]

                def wrap_delta(x):
                    return (x + 0.5 * delta) % delta - 0.5 * delta

                for pos in mid_positions:
                    a = int(id_all[pos])
                    pos_idx = np.where(order == pos)[0][0]
                    pos_next = order[(pos_idx + 1) % len(order)]
                    b = int(id_all[pos_next])

                    pha, phb = ph_all[pos], ph_all[pos_next]
                    wa, wb = w_all[pos], w_all[pos_next]
                    err_phase = wrap_delta((pha + 0.5 * wa) - (phb - 0.5 * wb))
                    Fm[a, j] += -k_pack * err_phase
                    Fm[b, j] += +k_pack * err_phase

        # Back to native + bending + damping
        for i in range(M):
            Fi = np.interp(xs_list[i], x_master, Fm[i], left=0.0, right=0.0)

            th = th_list[i]
            lap = np.zeros_like(th)
            if th.size > 2:
                lap[1:-1] = th[:-2] - 2.0 * th[1:-1] + th[2:]

            # IMPORTANT: do not bias endpoints toward zero slope
            if th.size > 1:
                lap[0] = 0.0
                lap[-1] = 0.0

            Fi += k_bend * lap

            # ------------------------------------------------------------
            # NEW: Soft helix-angle anchoring everywhere (slope bias)
            # Bias local dθ/dx toward the intended placement profile without
            # pinning absolute θ(x), so packing can still shift channels.
            # ------------------------------------------------------------
            if k_slope > 0.0 and desired_slope[i] is not None and th.size > 2:
                s_des_nodes = np.asarray(desired_slope[i], dtype=float)
                if s_des_nodes.shape[0] != th.shape[0]:
                    raise ValueError(
                        f"desired_slope[{i}] must have same length as theta array"
                    )

                x = xs_list[i]
                dx = np.diff(x)

                # segment slopes and desired slopes at segment midpoints
                s = np.diff(th) / (dx + 1e-16)
                s_des = 0.5 * (s_des_nodes[:-1] + s_des_nodes[1:])
                e = s - s_des

                # Discrete gradient of E = 0.5*Σ (dx_k * e_k^2) with w_k=dx_k
                g = np.zeros_like(th)
                dx = np.diff(x)
                w = dx / (np.mean(dx) + 1e-16)  # normalized segment weights ~ O(1)
                # e is length N-1, w is length N-1
                g[1:-1] = w[:-1] * e[:-1] - w[1:] * e[1:]


                # Do not fight hard endpoint slope constraints
                g[0] = 0.0
                g[-1] = 0.0

                # Force is negative gradient
                Fi += -k_slope * g

            v = v_list[i]
            v = (1.0 - gamma_damp) * v + dt * Fi
            th_new = th + dt * v
            
            v_list[i] = v
            th_list[i] = th_new

        # Enforce endpoint slopes after each update
        _apply_hard_endpoint_slopes()

        _record()

        if np.max(np.abs(np.concatenate(th_list) - prev_concat)) < tol:
            break

    out = []
    for i in range(M):
        cl = centerlines[i].copy()
        cl[:, 2] = th_list[i]
        out.append(cl)
    return out


# ---------------------------------------------------------------------------
# Geometry helpers for cross-sections
# ---------------------------------------------------------------------------

def cyl_centerline_to_cartesian(cl: np.ndarray) -> np.ndarray:
    """Convert a (N,3) [x, r, theta] centerline to Cartesian [X, Y, Z]."""
    x = cl[:, 0]
    r = cl[:, 1]
    th = cl[:, 2]
    y = r * np.cos(th)
    z = r * np.sin(th)
    return np.column_stack([x, y, z])


def estimate_tangents_cartesian(cl_cart: np.ndarray) -> np.ndarray:
    """Estimate unit tangents along a polyline in Cartesian coordinates."""
    N = cl_cart.shape[0]
    t = np.zeros_like(cl_cart)
    if N == 1:
        t[0] = np.array([1.0, 0.0, 0.0])
        return t

    for i in range(N):
        if i < N - 1:
            v = cl_cart[i + 1] - cl_cart[i]
        else:
            v = cl_cart[i] - cl_cart[i - 1]
        nrm = np.linalg.norm(v)
        if nrm < 1e-12:
            p = cl_cart[i]
            vr = np.array([0.0, p[1], p[2]])
            nrm = np.linalg.norm(vr)
            v = vr if nrm > 1e-12 else np.array([1.0, 0.0, 0.0])
        t[i] = v / (np.linalg.norm(v) + 1e-16)
    return t


def _blockage_at_x(circuit, x: float) -> float:
    """Return blockage_ratio at axial coordinate x for a circuit."""
    br = getattr(circuit, "blockage_ratio", None)
    if br is None:
        return 0.0

    br_arr = np.asarray(br, dtype=float)
    if br_arr.ndim == 0:
        return float(br_arr)

    x_dom = np.asarray(circuit.x_domain, dtype=float)
    return float(np.interp(x, x_dom, br_arr))


def get_section_scalars(circuit, x: float, r_center: float):
    """
    For a given circuit and axial location x, return scalar geometry fields.

    Returns a dict with:
    - r_center  : centerline radius (already stored in your centerline)
    - r_inner   : inner coolant radius (centerline + wall thickness)
    - r_outer   : outer coolant radius (r_inner + h(x))
    - theta_eff : effective wedge angle (blockage applied)
    plus the raw h, t_wall, and theta_raw.
    """
    h = float(circuit.channel_height(x))
    t_wall = float(circuit.total_thickness(x))
    #theta_raw = float(circuit.wedge_angle(x))
    theta_raw = float(circuit.local_sector_angle(x))
    br_here = _blockage_at_x(circuit, x)

    theta_eff = theta_raw * (1.0 - br_here)
    r_inner = r_center + t_wall
    r_outer = r_inner + h

    return {
        "h": h,
        "t_wall": t_wall,
        "theta_raw": theta_raw,
        "theta_eff": theta_eff,
        "r_center": r_center,
        "r_inner": r_inner,
        "r_outer": r_outer,
    }


def local_section_basis(cl_cyl_point: np.ndarray, t_hat: np.ndarray):
    """Build an orthonormal basis (e_r, e_phi) in the plane ⟂ tangent."""
    x_c, r_c, th_c = cl_cyl_point
    r_hat_world = np.array([0.0, np.cos(th_c), np.sin(th_c)], dtype=float)

    r_perp = r_hat_world - np.dot(r_hat_world, t_hat) * t_hat
    nr = np.linalg.norm(r_perp)
    if nr < 1e-12:
        r_perp = np.array([0.0, 1.0, 0.0]) - np.dot([0.0, 1.0, 0.0], t_hat) * t_hat
        nr = np.linalg.norm(r_perp)
    e_r = r_perp / nr

    e_phi = np.cross(t_hat, e_r)
    e_phi /= np.linalg.norm(e_phi) + 1e-16

    return e_r, e_phi

import numpy as np
import math
from typing import List, Tuple

Point3 = Tuple[float, float, float]


def _wrap_to_pi(a: float) -> float:
    """Wrap angle to (-pi, pi]."""
    a = (a + math.pi) % (2.0 * math.pi) - math.pi
    return a


def _arc_points_tangent_to_rays_yz(
    r: float,
    theta_center: float,
    theta_eff: float,
    n: int,
    *,
    use_major_arc: bool,
) -> List[Tuple[float, float]]:
    """
    Same construction as _arc_points_tangent_to_rays, but returns (y,z) points
    in the global y–z plane (no x).
    """
    if n < 2:
        raise ValueError("n must be >= 2")
    if not (0.0 < theta_eff < math.pi):
        raise ValueError("theta_eff must satisfy 0 < theta_eff < pi for a well-behaved convex construction.")

    a = 0.5 * theta_eff
    ca = math.cos(a)
    sa = math.sin(a)
    if abs(ca) < 1e-12:
        raise ValueError("theta_eff too close to pi; cos(theta_eff/2) ~ 0.")

    # Rotated 2D frame: bisector is +X axis. Endpoints at angles ±a, radius r.
    P_right = (r * ca, +r * sa)
    P_left  = (r * ca, -r * sa)

    d = r / ca
    C = (d, 0.0)

    R = r * math.tan(a)

    ang_left  = math.atan2(P_left[1]  - C[1], P_left[0]  - C[0])
    ang_right = math.atan2(P_right[1] - C[1], P_right[0] - C[0])

    delta = _wrap_to_pi(ang_right - ang_left)
    if use_major_arc:
        delta = (delta - 2.0 * math.pi) if delta >= 0 else (delta + 2.0 * math.pi)

    pts_rot: List[Tuple[float, float]] = []
    for i in range(n):
        t = i / (n - 1)
        ang = ang_left + t * delta
        xr = C[0] + R * math.cos(ang)
        yr = C[1] + R * math.sin(ang)
        pts_rot.append((xr, yr))

    # Rotate back by theta_center into global (y,z)
    cc = math.cos(theta_center)
    sc = math.sin(theta_center)

    yz: List[Tuple[float, float]] = []
    for (yr, zr) in pts_rot:
        y = yr * cc - zr * sc
        z = yr * sc + zr * cc
        yz.append((y, z))

    return yz

def rounded_slice_arcs_yz(
    r_inner: float,
    r_outer: float,
    theta_center: float,
    theta_eff: float,
    n: int = 8,
) -> Tuple[List[Tuple[float, float]], List[Tuple[float, float]]]:
    if r_inner <= 0 or r_outer <= 0:
        raise ValueError("r_inner and r_outer must be > 0.")
    if r_outer <= r_inner:
        raise ValueError("r_outer must be > r_inner.")

    arc_inner_yz = _arc_points_tangent_to_rays_yz(
        r=r_inner,
        theta_center=theta_center,
        theta_eff=theta_eff,
        n=n,
        use_major_arc=False,
    )
    arc_outer_yz = _arc_points_tangent_to_rays_yz(
        r=r_outer,
        theta_center=theta_center,
        theta_eff=theta_eff,
        n=n,
        use_major_arc=True,
    )
    return arc_inner_yz, arc_outer_yz

def build_rounded_section_points(
    P0: np.ndarray,
    t_hat: np.ndarray,
    cl_cyl_point: np.ndarray,
    geom: dict,
    n_arc_pts: int,
):
    """
    Build rounded wedge-like ring segment by rigidly embedding the original
    cross-section geometry into the plane normal to t_hat.

    No projection to cylinder-plane intersection.
    """
    x_c, r_center, theta_center = cl_cyl_point
    r_inner = float(geom["r_inner"])
    r_outer = float(geom["r_outer"])
    theta_eff = float(geom["theta_eff"])

    # plane basis (perpendicular to tangent)
    e_r, e_phi = local_section_basis(cl_cyl_point, t_hat)

    inner_uv, outer_uv = rounded_slice_arcs_uv(
        r_inner=r_inner,
        r_outer=r_outer,
        r_center=float(r_center),
        theta_center=float(theta_center),
        theta_eff=float(theta_eff),
        n=n_arc_pts,
    )

    # Winding: inner forward, outer reversed
    outer_uv = outer_uv[::-1]

    inner_pts = plane_uv_to_xyz(P0, e_r, e_phi, inner_uv[:, 0], inner_uv[:, 1])
    outer_pts = plane_uv_to_xyz(P0, e_r, e_phi, outer_uv[:, 0], outer_uv[:, 1])

    return np.vstack([inner_pts, outer_pts])


def build_sections_for_centerlines(
    centerlines,
    cid_full,
    thrust_chamber,
    *,
    stride: int = 1,
    n_arc_pts: int = 8,
):
    """
    Build 3D cross-section polygons for each centerline, using the actual
    CoolingCircuit + cross_section definitions.

    Parameters
    ----------
    centerlines :
        List of (Ni,3) centerlines in cylindrical [x, r, theta].
    cid_full :
        List/array of circuit indices for each centerline.
    thrust_chamber :
        Thrust chamber providing cooling circuit geometry.
    stride : int
        Subsampling stride along centerlines (increase to reduce mesh size).
    n_arc_pts : int
        Number of angular samples on inner and outer arcs.

    Returns
    -------
    sections_per_channel : list[list[np.ndarray]]
        For each channel, a list of (Nk,3) section polygons in XYZ.
    """
    circuits = thrust_chamber.cooling_circuits
    cid_full = np.asarray(cid_full, dtype=int)

    sections_per_channel: list[list[np.ndarray]] = []

    for i_cl, cl in enumerate(centerlines):
        cl_cyl = np.asarray(cl, float)
        cl_cart = cyl_centerline_to_cartesian(cl_cyl)
        T = estimate_tangents_cartesian(cl_cart)

        circuit_idx = cid_full[i_cl]
        circuit = circuits[circuit_idx]

        Ni = cl_cyl.shape[0]
        sections: list[np.ndarray] = []

        for i in range(0, Ni, stride):
            P0 = cl_cart[i]
            t_vec = T[i]
            t_norm = np.linalg.norm(t_vec)
            if t_norm < 1e-12:
                continue
            t_hat = t_vec / t_norm

            x_c = float(cl_cyl[i, 0])
            r_center = float(cl_cyl[i, 1])
            theta_center = float(cl_cyl[i, 2])

            geom = get_section_scalars(circuit, x_c, r_center)
            cl_cyl_point = np.array([x_c, r_center, theta_center], dtype=float)
            cs = circuit.cross_section

            if isinstance(cs, CrossSectionSquared):
                section_xyz = build_squared_section_points(
                    P0, t_hat, cl_cyl_point, geom, n_arc_pts=n_arc_pts
                )
            elif isinstance(cs, (CrossSectionRounded, CrossSectionRoundedInternal)):
                section_xyz = build_rounded_section_points(
                    P0, t_hat, cl_cyl_point, geom, n_arc_pts=n_arc_pts
                )
            else:
                section_xyz = build_squared_section_points(
                    P0, t_hat, cl_cyl_point, geom, n_arc_pts=n_arc_pts
                )

            sections.append(section_xyz)

        sections_per_channel.append(sections)

    return sections_per_channel


# ---------------------------------------------------------------------------
# PyVista lofting
# ---------------------------------------------------------------------------

def loft_sections_to_polydata(sections, *, cap_ends: bool = False) -> pv.PolyData:
    """
    Loft a list of cross-section polygons into a surface.

    Parameters
    ----------
    sections : list of (Ni, 3) arrays
        All sections must have the same number of points N,
        ordered consistently around the perimeter.
    cap_ends : bool
        If True, add a polygon at each end to close the channel.

    Returns
    -------
    pv.PolyData
        Lofted surface mesh.
    """
    sections = [np.asarray(sec, float) for sec in sections]
    if len(sections) < 2:
        raise ValueError("Need at least 2 sections to loft")

    N = sections[0].shape[0]
    for sec in sections:
        if sec.shape != (N, 3):
            raise ValueError("All sections must have shape (N,3) with the same N")

    M = len(sections)
    points = np.vstack(sections)

    def vid(i_sec, j_pt):
        return i_sec * N + j_pt

    faces = []
    for i_sec in range(M - 1):
        for j in range(N):
            j_next = (j + 1) % N
            v0 = vid(i_sec, j)
            v1 = vid(i_sec, j_next)
            v2 = vid(i_sec + 1, j_next)
            v3 = vid(i_sec + 1, j)
            faces.extend([4, v0, v1, v2, v3])

    if cap_ends:
        start_face = [N]
        start_face.extend(vid(0, j) for j in range(N))
        faces.extend(start_face)

        end_face = [N]
        end_face.extend(vid(M - 1, j) for j in reversed(range(N)))
        faces.extend(end_face)

    faces = np.asarray(faces, dtype=np.int64)
    return pv.PolyData(points, faces)


# ---------------------------------------------------------------------------
# PyVista entry point
# ---------------------------------------------------------------------------

def make_engine_3d(
    thrust_chamber,
    *,
    stride: int = 2,
    n_arc_pts: int = 8,
    show: bool = True,
):
    """
    Build a PyVista visualization of cooling channels using real cross-sections.

    Pipeline:
    1. Collect a minimal wedge of channels from the thrust chamber.
    2. Redistribute θ(x) to pack channels without overlap.
    3. Tile the wedge around the full circumference.
    4. For each centerline, build 3D cross-section polygons.
    5. Loft each channel to a PolyData and add to a PyVista plotter.

    Parameters
    ----------
    thrust_chamber
        Thrust chamber object with cooling circuit definitions.
    stride : int
        Subsampling stride along each centerline (increase to reduce mesh size).
    n_arc_pts : int
        Number of angular samples used for inner/outer arcs.
    show : bool
        If True, open the interactive PyVista window.

    Returns
    -------
    cl_full : list[np.ndarray]
        Full set of cylindrical centerlines [x, r, theta].
    w_full : list[np.ndarray]
        Channel width arrays for each centerline.
    cid_full : list[int]
        Circuit index per centerline.
    meshes : list[pv.PolyData]
        Lofted meshes for each channel.
    """

    print(f"\nStarting cooling channel visualisation")
    cl_sub, w_sub, cid_sub, meta = collect_channels(thrust_chamber, jitter_divisor=100.0)    

    g = meta["g"]
    wedge_angle = meta["wedge_angle"]

    start = time.time()
    recorder: list[list[np.ndarray]] = []

    # --- compute endpoint slopes from the *placement* for each sub-centerline ---
    slope0 = []
    slopeN = []
    slope_profile = []   # NEW
    contour = thrust_chamber.contour
    circuits = thrust_chamber.cooling_circuits

    for cl, cid in zip(cl_sub, cid_sub):
        circuit = circuits[int(cid)]
        placement = circuit.placement

        x0 = float(cl[0, 0])
        xN = float(cl[-1, 0])

        slope0.append(float(placement.dtheta_dx(x0, contour)))
        slopeN.append(float(placement.dtheta_dx(xN, contour)))

        xs = np.asarray(cl[:, 0], float)
        slope_profile.append(
            np.asarray([placement.dtheta_dx(float(x), contour) for x in xs], dtype=float)
        )

    cl_sub_solved = redistribute_channels(
        cl_sub,
        w_sub,
        k_pack=2.0,
        k_bend=0.1,
        k_slope=0.05,              
        desired_slope=slope_profile,  
        dt=0.15,
        tol=1e-5,
        max_iters=800,
        gamma_damp=0.4,
        wedge_angle=(wedge_angle if g > 1 else None),
        recorder=recorder,
        endpoint_slope0=slope0,
        endpoint_slopeN=slopeN,
    )


    cl_full, w_full, cid_full = tile_centerlines(
        cl_sub_solved, w_sub, cid_sub, max(g, 1), wedge_angle
    )


    sections_per_channel = build_sections_for_centerlines(
        cl_full,
        cid_full,
        thrust_chamber,
        stride=stride,
        n_arc_pts=n_arc_pts,
    )

    end = time.time()
    print(f"Channel redistribution and tiling finished in {end - start:.3f} s")
    print(f"Starting plotter")
    start = time.time()
    
    plotter = pv.Plotter()
    meshes = []
    cid_full = np.asarray(cid_full, int)

    palette = [
        "#1f77b4",
        "#ff7f0e",
        "#2ca02c",
        "#d62728",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#7f7f7f",
        "#bcbd22",
        "#17becf",
    ]

    for i_cl, sections in enumerate(sections_per_channel):
        if len(sections) < 2:
            continue

        mesh = loft_sections_to_polydata(sections, cap_ends=False)
        meshes.append(mesh)

        ci = int(cid_full[i_cl])
        color = palette[ci % len(palette)]

        plotter.add_mesh(
            mesh,
            color=color,
            smooth_shading=True,
            show_edges=False,
            name=f"channel_{i_cl}_circ_{ci}",
        )
    plotter.hide_axes()
    #plotter.show_bounds(xtitle="x", ytitle="y", ztitle="z")
    plotter.set_background("white")
    # Export scene to a glTF file with inline data (buffers embedded) :contentReference[oaicite:1]{index=1}
    with tempfile.TemporaryDirectory() as td:
        gltf_path = Path(td) / "engine.gltf"
        plotter.export_gltf(str(gltf_path), inline_data=True)
        gltf_text = gltf_path.read_text(encoding="utf-8")

    script_dir = Path(__file__).resolve().parent
    assets_dir = script_dir / "assets"
    viewer_src = _threejs_gltf_viewer_data_url(gltf_text, assets_dir=assets_dir, title="Engine 3D (glTF/three.js)")
    end = time.time()

    print(f"Plotter finished in {end - start:.3f} s")

    return plotter, viewer_src #cl_full, w_full, cid_full, meshes

def plane_uv_to_xyz(
    P0: np.ndarray,
    e_r: np.ndarray,
    e_phi: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
) -> np.ndarray:
    """Map section-plane offsets (u,v) into XYZ: X = P0 + u*e_r + v*e_phi."""
    u = np.asarray(u, float)
    v = np.asarray(v, float)
    return P0[None, :] + u[:, None] * e_r[None, :] + v[:, None] * e_phi[None, :]


def yz_to_uv_offsets(
    y: np.ndarray,
    z: np.ndarray,
    y0: float,
    z0: float,
    theta_center: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Convert world (y,z) points into local (u,v) offsets relative to the
    centerline point (y0,z0), expressed in the (r_hat, phi_hat) basis at theta_center.

    u: radial offset (along r_hat)
    v: circumferential offset (along phi_hat)
    """
    y = np.asarray(y, float)
    z = np.asarray(z, float)

    dy = y - float(y0)
    dz = z - float(z0)

    c = np.cos(theta_center)
    s = np.sin(theta_center)

    # r_hat at theta_center is [0, c, s]
    # phi_hat at theta_center is [0, -s, c]
    u = dy * c + dz * s
    v = -dy * s + dz * c
    return u, v

def rounded_slice_arcs_uv(
    r_inner: float,
    r_outer: float,
    r_center: float,
    theta_center: float,
    theta_eff: float,
    n: int = 8,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Build inner/outer rounded arcs as (u,v) offsets in a local 2D cross-section frame,
    relative to the centerline point at radius r_center and angle theta_center.

    Returns
    -------
    inner_uv : (n,2) array
    outer_uv : (n,2) array
    """
    inner_yz, outer_yz = rounded_slice_arcs_yz(
        r_inner=r_inner,
        r_outer=r_outer,
        theta_center=theta_center,
        theta_eff=theta_eff,
        n=n,
    )

    # Centerline point in world yz
    y0 = float(r_center * np.cos(theta_center))
    z0 = float(r_center * np.sin(theta_center))

    inner_yz = np.asarray(inner_yz, float)
    outer_yz = np.asarray(outer_yz, float)

    u_in, v_in = yz_to_uv_offsets(inner_yz[:, 0], inner_yz[:, 1], y0, z0, theta_center)
    u_out, v_out = yz_to_uv_offsets(outer_yz[:, 0], outer_yz[:, 1], y0, z0, theta_center)

    inner_uv = np.column_stack([u_in, v_in])
    outer_uv = np.column_stack([u_out, v_out])
    return inner_uv, outer_uv

def build_squared_section_points(
    P0: np.ndarray,
    t_hat: np.ndarray,
    cl_cyl_point: np.ndarray,
    geom: dict,
    n_arc_pts: int,
):
    """
    Build squared wedge-like ring segment, but by rigidly embedding the original
    x-normal cross-section shape into the plane normal to t_hat.

    No cylinder-plane intersection, no "solve x".
    """
    x_c, r_center, theta_center = cl_cyl_point
    r_inner = float(geom["r_inner"])
    r_outer = float(geom["r_outer"])
    theta_eff = float(geom["theta_eff"])

    # plane basis (perpendicular to tangent)
    e_r, e_phi = local_section_basis(cl_cyl_point, t_hat)

    # Centerline point in world yz (same x as P0)
    y0 = float(r_center * np.cos(theta_center))
    z0 = float(r_center * np.sin(theta_center))

    # Angular samples about axis (original definition)
    dth_in = np.linspace(-theta_eff / 2.0, +theta_eff / 2.0, n_arc_pts)
    dth_out = dth_in[::-1]

    th_in = theta_center + dth_in
    th_out = theta_center + dth_out

    # Absolute yz points on circles (original x-normal section geometry)
    y_in = r_inner * np.cos(th_in)
    z_in = r_inner * np.sin(th_in)

    y_out = r_outer * np.cos(th_out)
    z_out = r_outer * np.sin(th_out)

    # Convert to (u,v) offsets relative to the centerline point
    u_in, v_in = yz_to_uv_offsets(y_in, z_in, y0, z0, theta_center)
    u_out, v_out = yz_to_uv_offsets(y_out, z_out, y0, z0, theta_center)

    inner_pts = plane_uv_to_xyz(P0, e_r, e_phi, u_in, v_in)
    outer_pts = plane_uv_to_xyz(P0, e_r, e_phi, u_out, v_out)

    return np.vstack([inner_pts, outer_pts])


def _read_text(path: Path) -> str:
    txt = path.read_text(encoding="utf-8")
    # Prevent accidental </script> termination inside inline scripts
    return txt.replace("</script>", "</scr" + "ipt>")


def _threejs_gltf_viewer_data_url(
    gltf_text: str,
    *,
    assets_dir: Path,
    title: str = "3D Viewer",
    background: str = "#ffffff",
) -> str:
    """
    Offline, single-file: embeds three.js + OrbitControls + GLTFLoader inline,
    and embeds glTF JSON inline (base64) inside the iframe HTML.
    """
    import base64

    three_js = _read_text(assets_dir / "three.min.js")
    orbit_js = _read_text(assets_dir / "OrbitControls.js")
    gltf_loader_js = _read_text(assets_dir / "GLTFLoader.js")

    gltf_b64 = base64.b64encode(gltf_text.encode("utf-8")).decode("ascii")

    viewer_html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width,initial-scale=1" />
  <title>{title}</title>
  <style>
    html, body {{ height:100%; margin:0; background:{background}; overflow:hidden; }}
    #c {{ width:100%; height:100%; display:block; }}
    #hud {{
      position:absolute; top:10px; left:10px;
      padding:6px 10px; border-radius:8px;
      font:12px/1.2 system-ui, -apple-system, Segoe UI, Roboto, sans-serif;
      color:rgba(0,0,0,0.75);
      background:rgba(255,255,255,0.75);
      user-select:none;
      z-index: 10;
      white-space: pre;
    }}
  </style>
</head>
<body>
  <div id="hud">Drag: rotate • Right-drag: pan • Wheel: zoom</div>
  <canvas id="c"></canvas>

  <script id="gltf-b64" type="text/plain">{gltf_b64}</script>

  <!-- Offline-inlined libraries -->
  <script>{three_js}</script>
  <script>{orbit_js}</script>
  <script>{gltf_loader_js}</script>

  <script>
    (function () {{
      const hud = document.getElementById("hud");
      function fail(msg, err) {{
        hud.textContent = msg + (err ? ("\\n" + (err.message || err)) : "");
        console.error(err || msg);
      }}

      if (!window.THREE) return fail("three.js not available (inline load failed).");
      if (!THREE.OrbitControls) return fail("OrbitControls not available (inline load failed).");
      if (!THREE.GLTFLoader) return fail("GLTFLoader not available (inline load failed).");

      const canvas = document.getElementById("c");
      const renderer = new THREE.WebGLRenderer({{ canvas: canvas, antialias: true }});
      renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
      renderer.setSize(window.innerWidth, window.innerHeight, false);

      const scene = new THREE.Scene();
      scene.background = new THREE.Color("{background}");

      const camera = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 0.01, 1e6);
      camera.position.set(2, 1.2, 2);

      const controls = new THREE.OrbitControls(camera, renderer.domElement);
      controls.enableDamping = true;
      controls.dampingFactor = 0.08;

      scene.add(new THREE.AmbientLight(0xffffff, 0.65));
      const dir = new THREE.DirectionalLight(0xffffff, 0.85);
      dir.position.set(5, 8, 6);
      scene.add(dir);

      window.addEventListener("resize", () => {{
        renderer.setSize(window.innerWidth, window.innerHeight, false);
        camera.aspect = window.innerWidth / window.innerHeight;
        camera.updateProjectionMatrix();
      }});

        const b64 = document.getElementById("gltf-b64").textContent.trim();

        async function decodeGltfJsonFromB64(b64) {{
        // Offline: this is not a network fetch; it's decoding a data: URL in-memory.
        const dataUrl = "data:application/json;base64," + b64;
        const resp = await fetch(dataUrl);
        if (!resp.ok) throw new Error("Failed to decode glTF data URL: " + resp.status);
        return await resp.text(); // GLTFLoader.parse expects JSON text for .gltf
        }}

        const loader = new THREE.GLTFLoader();

        (async () => {{
        try {{
            hud.textContent = "Decoding glTF…";
            // Give the HUD a chance to paint
            await new Promise(r => setTimeout(r, 0));

            const gltfJsonText = await decodeGltfJsonFromB64(b64);

            hud.textContent = "Parsing model…";
            loader.parse(
            gltfJsonText,
            "",
            (gltf) => {{
                const root = gltf.scene || (gltf.scenes && gltf.scenes[0]);
                if (!root) return fail("Parsed glTF, but no scene found.");

                scene.add(root);

                root.traverse((obj) => {{
                if (!obj.isMesh) return;
                const mats = Array.isArray(obj.material) ? obj.material : [obj.material];
                for (const m of mats) {{
                    if (!m) continue;
                    m.side = THREE.DoubleSide;
                    m.needsUpdate = true;
                }}
                }});

                const box = new THREE.Box3().setFromObject(root);
                const size = box.getSize(new THREE.Vector3());
                const center = box.getCenter(new THREE.Vector3());
                const maxDim = Math.max(size.x, size.y, size.z);
                if (!isFinite(maxDim) || maxDim <= 0) return fail("Model bounds invalid (empty geometry?).");

                const fov = camera.fov * Math.PI / 180;
                let dist = (maxDim / 2) / Math.tan(fov / 2);
                dist *= 1.35;

                camera.position.set(center.x + dist, center.y + dist * 0.5, center.z + dist);
                camera.near = Math.max(dist / 1000, 0.001);
                camera.far  = dist * 1000;
                camera.updateProjectionMatrix();

                controls.target.copy(center);
                controls.update();

                hud.textContent = "Drag: rotate • Right-drag: pan • Wheel: zoom";
            }},
            (err) => fail("Failed to parse glTF.", err)
            );
        }} catch (e) {{
            return fail("Failed to decode glTF.", e);
        }}
        }})();


      function animate() {{
        controls.update();
        renderer.render(scene, camera);
        requestAnimationFrame(animate);
      }}
      animate();
    }})();
  </script>
</body>
</html>
"""
    viewer_b64 = base64.b64encode(viewer_html.encode("utf-8")).decode("ascii")
    return f"data:text/html;base64,{viewer_b64}"
