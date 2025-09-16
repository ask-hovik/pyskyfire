from __future__ import annotations
from dataclasses import dataclass
from abc import ABC, abstractmethod
import numpy as np

@dataclass(frozen=True)
class SectionProfiles:
    """All geometry profiles along the centerline (length N)."""
    h: np.ndarray            # channel height [N]
    theta: np.ndarray        # included angle [rad, N]
    t_wall: np.ndarray       # wall thickness [N]
    centerline: np.ndarray   # (N,3) [x, r, z] or whatever convention you use
    local_coords: np.ndarray # (N,3,3) orthonormal frames [t, n, b]
    blockage_ratio: np.ndarray

    def __post_init__(self):
        N = self.centerline.shape[0]
        for name in ("h", "theta", "t_wall", "blockage_ratio"):
            v = getattr(self, name)
            if v.shape[0] != N:
                raise ValueError(f"{name} length {v.shape[0]} != N {N}")
        if self.centerline.shape != (N, 3):
            raise ValueError("centerline must be (N,3)")
        if self.local_coords.shape != (N, 3, 3):
            raise ValueError("local_coords must be (N,3,3)")

class ChannelSection(ABC):
    """Common interface all cross-section shapes must implement."""
    def __init__(self, n_points: int = 16):
        self.n_points = int(n_points)

    # ---- Required API (vectorized over stations) -----------------------------
    @abstractmethod
    def A_coolant(self, prof: SectionProfiles) -> np.ndarray: ...
    @abstractmethod
    def Dh_coolant(self, prof: SectionProfiles) -> np.ndarray: ...
    @abstractmethod
    def P_thermal(self, prof: SectionProfiles) -> np.ndarray: ...
    @abstractmethod
    def P_coolant(self, prof: SectionProfiles) -> np.ndarray: ...
    @abstractmethod
    def compute_cross_section(self, prof: SectionProfiles, i: int) -> int: ...
    # -------------------------------------------------------------------------

# ===================== Squared implementation ===============================
class CrossSectionSquared(ChannelSection):
    def __init__(self, n_points: int = 8):
        super().__init__(n_points=n_points)

    def _theta_real(self, prof: SectionProfiles) -> np.ndarray:
        return prof.theta * (1 - prof.blockage_ratio)

    def A_coolant(self, prof: SectionProfiles) -> np.ndarray:
        r = prof.centerline[:, 1]
        r_inner = r + prof.t_wall
        r_outer = r + prof.t_wall + prof.h
        th = self._theta_real(prof)
        A_sector = (th/2) * (r_outer**2 - r_inner**2)
        return A_sector

    def Dh_coolant(self, prof: SectionProfiles) -> np.ndarray:
        A = self.A_coolant(prof)
        r = prof.centerline[:, 1]
        r_inner = r + prof.t_wall
        r_outer = r + prof.t_wall + prof.h
        th = self._theta_real(prof)
        P = r_inner*th + r_outer*th + 2*prof.h
        return 4.0 * A / P

    def P_thermal(self, prof: SectionProfiles) -> np.ndarray:
        r = prof.centerline[:, 1]
        th = self._theta_real(prof)
        return r * th

    def P_coolant(self, prof: SectionProfiles) -> np.ndarray:
        r = prof.centerline[:, 1]
        th = self._theta_real(prof)
        r_inner = r + prof.t_wall
        return r_inner * th

    def compute_cross_section(self, prof: SectionProfiles, i: int):
        """
        Build a closed OCC wire (int tag) for station i by:
        1) creating a circle EDGE in XY@origin,
        2) applying an affine transform to place it at (P_i, t_i, n_i, b_i),
        3) wrapping the transformed edge into a wire.
        NOTE: gmsh must be initialized and a model added by the caller.
        """
        import gmsh # lazy import gmsh
        # -------- station data --------
        x_i, r_i, th_i = map(float, prof.centerline[i])   # cylindrical (x, r, theta)
        t_i = np.asarray(prof.local_coords[i, 0], float)  # tangent
        n_i = np.asarray(prof.local_coords[i, 1], float)  # normal (will be re-orthonormalized)
        b_i = np.asarray(prof.local_coords[i, 2], float)  # binormal

        # world-space center point for the section
        P_i = np.array([x_i, r_i * np.sin(th_i), r_i * np.cos(th_i)], float)
        

        # -------- affine matrix, transformation from world to local coordinate system --------
        bx, by, bz = map(float, b_i)
        nx, ny, nz = map(float, n_i)
        tx, ty, tz = map(float, t_i)
        Px, Py, Pz = map(float, P_i)
        A = [
            bx, nx, tx, Pz,
            by, ny, ty, Py,
            bz, nz, tz, Px,
        ]

        A = [
            1.0, 0.0, 0.0, Pz, # TODO: this rotation stuff is just too hard to figure out in an evening man
            0.0, 1.0, 0.0, Py, # I had to switch out x and z to get them oriented right for some reason??
            0.0, 0.0, 1.0, Px, ## Revisit this and get the rotation matrix right!
        ]

        th_val  = prof.theta[i]
        br      = prof.blockage_ratio[i]
        phi     = th_val * (1 - br)
        h       = prof.h[i]
        
        cntr = gmsh.model.occ.addPoint(-r_i, 0, 0)
        Ai   = gmsh.model.occ.addPoint(-1*(r_i - r_i*np.cos(phi/2)), -r_i*np.sin(phi/2), 0)
        Bi   = gmsh.model.occ.addPoint(-1*(r_i - r_i*np.cos(phi/2)), r_i*np.sin(phi/2), 0)
        Ao   = gmsh.model.occ.addPoint((r_i + h)*np.cos(phi/2) - r_i, -(r_i + h)*np.sin(phi/2), 0)
        Bo   = gmsh.model.occ.addPoint((r_i + h)*np.cos(phi/2) - r_i, (r_i + h)*np.sin(phi/2), 0)
        
        arc_in  = gmsh.model.occ.addCircleArc(Ai, cntr, Bi)
        arc_out = gmsh.model.occ.addCircleArc(Bo, cntr, Ao)  # reversed to close
        wall_p  = gmsh.model.occ.addLine(Bi, Bo)
        wall_m  = gmsh.model.occ.addLine(Ao, Ai)
        edges = [arc_in, wall_p, arc_out, wall_m]

        # Transform ALL edges (not the wire)
        gmsh.model.occ.affineTransform([(1, e) for e in edges], A)
        gmsh.model.occ.synchronize()

        # -------- wrap transformed edge into a wire and return --------
        w = gmsh.model.occ.addWire([*edges])
        return w


# ===================== Rounded implementation ===============================
class CrossSectionRounded(ChannelSection):
    def __init__(self, n_points: int = 16): # TODO: Does n-points have any funcion left? 
        super().__init__(n_points=n_points)

    def _beta_alpha(self, theta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        beta = (np.pi - theta) * 0.5
        alpha = np.pi - beta
        return beta, alpha

    def A_coolant(self, prof: SectionProfiles) -> np.ndarray:
        r = prof.centerline[:, 1]
        r1 = r * np.sin(prof.theta/2)
        r2 = (r + prof.h) * np.sin(prof.theta/2)
        beta, alpha = self._beta_alpha(prof.theta)
        L_side = prof.h * np.cos(prof.theta/2)

        A1 = 0.5 * (r1 - prof.t_wall)**2 * 2 * beta
        A2 = 0.5 * (r2 - prof.t_wall)**2 * 2 * alpha
        S1 = (r1 - prof.t_wall) * L_side
        S2 = (r2 - prof.t_wall) * L_side
        return A1 + A2 + S1 + S2

    def Dh_coolant(self, prof: SectionProfiles) -> np.ndarray:
        A = self.A_coolant(prof)
        r = prof.centerline[:, 1]
        r1 = r * np.sin(prof.theta/2)
        r2 = (r + prof.h) * np.sin(prof.theta/2)
        beta, alpha = self._beta_alpha(prof.theta)
        L_side = prof.h * np.cos(prof.theta/2)
        arc1 = (r1 - prof.t_wall) * 2 * beta
        arc2 = (r2 - prof.t_wall) * 2 * alpha
        P = arc1 + arc2 + 2 * L_side
        return 4.0 * A / P

    def P_thermal(self, prof: SectionProfiles) -> np.ndarray:
        r = prof.centerline[:, 1]
        r1 = r * np.sin(prof.theta/2)
        angle = np.radians(112.0)     # your chosen effective contact angle
        return r1 * angle

    def P_coolant(self, prof: SectionProfiles) -> np.ndarray:
        r = prof.centerline[:, 1]
        r1 = r * np.sin(prof.theta/2)
        angle = np.radians(112.0)
        return (r1 - prof.t_wall) * angle
    
    def compute_cross_section(self, prof: SectionProfiles, i: int) -> int:
        """
        Build a closed OCC wire (int tag) for station i by:
        1) creating a circle EDGE in XY@origin,
        2) applying an affine transform to place it at (P_i, t_i, n_i, b_i),
        3) wrapping the transformed edge into a wire.
        NOTE: gmsh must be initialized and a model added by the caller.
        """
        import gmsh # lazy import
        # -------- station data --------
        x_i, r_i, th_i = map(float, prof.centerline[i])   # cylindrical (x, r, theta)
        t_i = np.asarray(prof.local_coords[i, 0], float)  # tangent
        n_i = np.asarray(prof.local_coords[i, 1], float)  # normal (will be re-orthonormalized)
        b_i = np.asarray(prof.local_coords[i, 2], float)  # binormal

        # world-space center point for the section
        P_i = np.array([x_i, r_i * np.sin(th_i), r_i * np.cos(th_i)], float)
        

        # -------- affine matrix, transformation from world to local coordinate system --------
        bx, by, bz = map(float, b_i)
        nx, ny, nz = map(float, n_i)
        tx, ty, tz = map(float, t_i)
        Px, Py, Pz = map(float, P_i)
        A = [
            bx, nx, tx, Pz, # this is the real rotation matrix. But I need to redo the 
            by, ny, ty, Py,
            bz, nz, tz, Px,
        ]

        A = [
            1.0, 0.0, 0.0, Pz, # TODO: this rotation stuff is just too hard to figure out in an evening man
            0.0, 1.0, 0.0, Py, # I had to switch out x and z to get them oriented right for some reason??
            0.0, 0.0, 1.0, Px, ## Revisit this and get the rotation matrix right!
        ]

        th_val  = prof.theta[i]
        t       = prof.t_wall[i]
        phi     = th_val - t/r_i
        h       = prof.h[i]

        # calculate some distances
        v = r_i + h
        q = v*np.cos(phi/2)
        xo = q-r_i
        yo = (r_i + h)*np.sin(phi/2)
        rotr = np.sqrt((v - q)**2 + yo**2)
        xotr = h + rotr

        # set points
        p_cntr_i = [0, 0, 0]
        p_cntr_o = [h, 0, 0]
        p_Ai   = [-1*(r_i - r_i*np.cos(phi/2)), -r_i*np.sin(phi/2), 0]
        p_Bi   = [-1*(r_i - r_i*np.cos(phi/2)), r_i*np.sin(phi/2), 0]
        p_Ao   = [xo, -yo, 0]
        p_Bo   = [xo, yo, 0]
        p_outr = [xotr, 0, 0]

        # Make gmsh points
        cntr_i = gmsh.model.occ.addPoint(*p_cntr_i)
        cntr_o = gmsh.model.occ.addPoint(*p_cntr_o)
        Ai   = gmsh.model.occ.addPoint(*p_Ai)
        Bi   = gmsh.model.occ.addPoint(*p_Bi)
        Ao   = gmsh.model.occ.addPoint(*p_Ao)
        Bo   = gmsh.model.occ.addPoint(*p_Bo)
        outr = gmsh.model.occ.addPoint(*p_outr)        

        # make gmsh elements
        arc_in  = gmsh.model.occ.addCircleArc(Ai, cntr_i, Bi)
        arc_out1 = gmsh.model.occ.addCircleArc(Bo, cntr_o, outr)  # reversed to close
        arc_out2 = gmsh.model.occ.addCircleArc(outr, cntr_o, Ao)  # reversed to close
        wall_p  = gmsh.model.occ.addLine(Bi, Bo)
        wall_m  = gmsh.model.occ.addLine(Ao, Ai)
        edges = [arc_in, wall_p, arc_out1, arc_out2, wall_m]

        # Transform to local coordinate system
        gmsh.model.occ.affineTransform([(1, e) for e in edges], A)
        gmsh.model.occ.synchronize()

        # -------- wrap transformed edge into a wire and return --------
        w = gmsh.model.occ.addWire([*edges])
        return w
