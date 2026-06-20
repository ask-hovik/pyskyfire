from __future__ import annotations

import os
import math
import warnings
from typing import Optional

# Set proper trans.lib for CEA_Wrap
script_dir = os.path.dirname(os.path.abspath(__file__))
trans_path = os.path.join(script_dir, "data", r"trans.lib")
os.environ["CEA_TRANS_LIB"] = trans_path
import CEA_Wrap as cea

from dataclasses import dataclass, asdict
import numpy as np
from bisect import bisect_left


class Aerothermodynamics:  
    """Aerothermodynamic property precomputation and lookup along an engine contour.

    This class computes 2-D thermo/fluid property maps along a prescribed
    nozzle/combustor contour using **CEA_Wrap** and exposes fast interpolating
    getters. At each axial position ``x`` the equilibrium state (column 0) is
    evaluated at the local area ratio; additional columns tabulate thermodynamic
    properties over a temperature grid at **fixed local static pressure**. Lookups
    support:
    - equilibrium at a given ``x`` (no ``T``/``h`` provided),
    - TP evaluation at specified ``T`` and local pressure ``p(x)``,
    - HP evaluation at specified mixture enthalpy ``h`` and local pressure ``p(x)``
    (using an adjustable-enthalpy reactant representation).

    Parameters
    ----------
    optimum : dict[str, float | str]
        Calibrated/derived inputs (e.g., ``fu``, ``ox``, ``MR``, ``p_c``, ``F``,
        ``eps``, ``L_star``, ``c_star``, inlet temperatures, etc.). Typically
        produced by :meth:`from_F_eps_Lstar`.
    chemrep_map : dict[str, str], optional
        Mapping from reactant *names* to ``chemical_representation`` strings used by
        CEA for HP evaluations with enthalpy offsets (exploded formula like
        ``"C 2 H 6 O 1"``). Required only when calling HP-based getters (``h=...``).

    Attributes
    ----------
    # Inputs / design-point scalars
    fu : Any
        User’s fuel specification (mixture container understood by CEA_Wrap).
    ox : Any
        User’s oxidizer specification (mixture container understood by CEA_Wrap).
    MR : float
        Mixture ratio ``m_ox / m_fu`` [-].
    p_c : float
        Chamber pressure [Pa].
    F : float
        Target thrust [N].
    eps : float
        Exit-to-throat area ratio ``A_e/A_t`` [-].
    L_star : float
        Characteristic chamber length [m].
    T_fu_in, T_ox_in : float
        Inlet temperatures for fuel/oxidizer [K].
    p_amb : float
        Ambient/static back pressure used for CF/Isp_amb calculations [Pa].
    npts : int
        Number of axial nodes used for precomputation along the contour [-].

    # Derived performance at design point
    c_star : float
        Characteristic velocity [m/s].
    Isp_vac, Isp_amb, Isp_SL, Isp_ideal_amb : float
        Vacuum / ambient / sea-level / perfectly expanded specific impulses [s].
    CF_vac, CF_amb, CF_SL : float
        Thrust coefficients [-].
    mdot, mdot_fu, mdot_ox : float
        Total/fuel/oxidizer mass flows [kg/s].
    A_t, A_e : float
        Throat and exit areas [m²].
    r_t, r_e : float
        Throat and exit radii [m].
    V_c : float
        Chamber volume estimated from ``L_star`` [m³].
    t_stay : float
        Mean residence time in chamber [s].

    # Precomputed grids (after :meth:`compute_aerothermodynamics`)
    x_nodes : ndarray
        Monotone axial grid spanning ``contour.xs[0]`` → ``contour.xs[-1]`` [m].
    Nt : int
        Number of temperature samples per axial station [-].
    T_grid : ndarray, shape (Nx, Nt)
        Temperature grid per axial row [K] (col 0 equals local equilibrium T).
    M_map, T_map, p_map, rho_map, cp_map, gamma_map, h_map, a_map, mu_map, k_map, Pr_map : ndarray
        Property maps; column 0 is the local equilibrium; remaining columns are TP
        samples at fixed local pressure. Units: M [-], T [K], p [bar], rho [kg/m³],
        cp [kJ/kg-K], gamma [-], h [kJ/kg], a [m/s], mu [Pa·s], k [W/m-K], Pr [-].
    X_map : list[dict[str, float]]
        Equilibrium product mole fractions per axial node (unnormalized dicts).

    See Also
    --------
    CEA_Wrap.RocketProblem
        Equilibrium nozzle/chamber solves used for the equilibrium column.
    CEA_Wrap.TPProblem
        Temperature-pressure state solves used for the TP columns and fallbacks.
    CEA_Wrap.HPProblem
        Enthalpy-pressure equilibrium solves used when ``h`` is specified.

    Notes
    -----
    - **Units:** Inputs use SI (Pa, K, N). CEA_Wrap expects/returns some fields in
    *psi* / *bar*; all conversions are handled internally. Stored ``p_map`` is in
    **bar** (matching CEA); :meth:`get_p` returns **Pa**.
    - **Equilibrium column:** column 0 in each map corresponds to the equilibrium
    state at the local area ratio (subsonic for ``x < 0``; supersonic for
    ``x >= 0`` with your sign convention).
    - **Interpolation:** Lookups first try a bilinear map interpolation on
    ``(x, T)``. If the requested ``T`` falls outside the tabulated range, a live
    TP solve at ``p(x)`` is performed. For ``h`` queries an HP solve at ``p(x)``
    is performed and requires valid ``chemrep_map`` entries for all reactants.
    """

    def __init__(self, optimum: dict[str, float | str], chemrep_map: Optional[dict[str, str]] = None):
        """Add all values in the optimum dict to self"""
        for key, value in optimum.items():
            setattr(self, key, value)
        self.optimum = optimum
        self.chemrep_map = chemrep_map or {}

    @classmethod
    def from_F_eps_Lstar(cls, fu, ox, MR, p_c, F, eps, L_star, T_fu_in=298.15, T_ox_in=298.15, p_amb=1.013e5, npts=15):
        """Construct from thrust, area ratio, and L* at a given chamber pressure.

        This helper solves a CEA rocket problem at the specified design point and
        assembles the ``optimum`` dict used to initialize :class:`Aerothermodynamics`.

        Parameters
        ----------
        fu : Any
            Fuel mixture descriptor (must expose ``propellants`` and ``fractions`` and
            be consumable by ``CEA_Wrap.Fuel``).
        ox : Any
            Oxidizer mixture descriptor (same structural expectations as ``fu``).
        MR : float
            Mixture ratio ``m_ox / m_fu`` [-].
        p_c : float
            Chamber pressure [Pa].
        F : float
            Target thrust [N].
        eps : float
            Nozzle area ratio ``A_e / A_t`` [-].
        L_star : float
            Characteristic chamber length [m].
        T_fu_in : float, default 298.15
            Fuel inlet temperature [K].
        T_ox_in : float, default 298.15
            Oxidizer inlet temperature [K].
        p_amb : float, default 1.013e5
            Ambient/static pressure for CF/Isp_amb calculations [Pa].
        npts : int, default 15
            Number of axial stations to precompute along the contour.

        Returns
        -------
        Aerothermodynamics
            Initialized instance with populated design-point performance and all fields
            necessary for :meth:`compute_aerothermodynamics`.

        Notes
        -----
        - The CEA call assumes perfect expansion for ``Isp_vac`` and uses standard
        thrust-coefficient relations to compute ``Isp_amb`` and ``Isp_SL``.
        - Areas and chamber volume are derived from ``c_star``, mass flow, ``L_star``,
        and mixture density at chamber conditions.
        """

        fus = []
        oxs = []
        for prop, frac in zip(fu.propellants, fu.fractions):
            fus.append(cea.Fuel(prop, wt=frac, temp=T_fu_in))

        for prop, frac in zip(ox.propellants, ox.fractions):
            oxs.append(cea.Oxidizer(prop, wt=frac, temp=T_ox_in))
        
        rp = cea.RocketProblem(o_f = MR, 
                               pressure=p_c*0.000145038, # Convert Pa to psi
                               materials=[*oxs, *fus], 
                               sup=eps) 
        R = rp.run()
        
        #names = cea.ThermoInterface.keys()
        #print(names)

        c_star = float(getattr(R, "cstar"))           
        Isp_ideal_amb = float(getattr(R, "isp"))      # perfectly expanded Isp
        Isp_vac = float(getattr(R, "ivac"))           
        rho_c = float(getattr(R, "c_rho"))
        T_c = float(getattr(R, "c_t"))
        T_t = float(getattr(R, "t_t"))
        p_t = float(getattr(R, "t_p"))

        g = 9.81
        p_SL = 1.01325e5 # sea level pressure
        mdot = F/(Isp_vac*g)
        mdot_fu = mdot/(1+MR)
        mdot_ox = mdot - mdot_fu

        A_t = c_star*mdot/p_c
        r_t = np.sqrt(A_t/np.pi)
        t_stay = L_star*A_t*rho_c/mdot
        V_c = mdot*t_stay/rho_c

        A_e = A_t*eps
        r_e = np.sqrt(A_e/np.pi)

        # calculate ambient Isp
        CF_vac = Isp_vac * g / c_star
        CF_amb = CF_vac - (p_amb / p_c) * (A_e/A_t)
        Isp_amb = CF_amb*c_star/g

        CF_SL = CF_vac - (p_SL / p_c) * (A_e/A_t)
        Isp_SL = CF_SL*c_star/g

        optimum = dict(
            fu=fu, ox=ox, MR=MR, p_c=p_c, p_t=p_t, T_c=T_c, T_t=T_t, F=F, eps=eps, L_star=L_star, c_star=c_star,
            p_amb=p_amb, Isp_ideal_amb=Isp_ideal_amb, Isp_vac=Isp_vac, Isp_amb=Isp_amb,
            Isp_SL=Isp_SL, CF_vac=CF_vac, CF_amb=CF_amb, CF_SL=CF_SL, mdot=mdot,
            mdot_fu=mdot_fu, mdot_ox=mdot_ox, t_stay=t_stay, A_t=A_t, A_e=A_e,
            r_t=r_t, r_e=r_e, V_c=V_c, T_fu_in=T_fu_in, T_ox_in=T_ox_in, npts=npts,
        )

        return cls(optimum)
    
    @classmethod
    def from_F_pe_Lstar(cls, fu, ox, MR, p_c, F, p_e, L_star, T_fu_in=298.15, T_ox_in=298.15, p_amb=1.013e5, npts=15):
        """Construct from thrust, area ratio, and L* at a given chamber pressure.

        This helper solves a CEA rocket problem at the specified design point and
        assembles the ``optimum`` dict used to initialize :class:`Aerothermodynamics`.

        Parameters
        ----------
        fu : Any
            Fuel mixture descriptor (must expose ``propellants`` and ``fractions`` and
            be consumable by ``CEA_Wrap.Fuel``).
        ox : Any
            Oxidizer mixture descriptor (same structural expectations as ``fu``).
        MR : float
            Mixture ratio ``m_ox / m_fu`` [-].
        p_c : float
            Chamber pressure [Pa].
        F : float
            Target thrust [N].
        eps : float
            Nozzle area ratio ``A_e / A_t`` [-].
        L_star : float
            Characteristic chamber length [m].
        T_fu_in : float, default 298.15
            Fuel inlet temperature [K].
        T_ox_in : float, default 298.15
            Oxidizer inlet temperature [K].
        p_amb : float, default 1.013e5
            Ambient/static pressure for CF/Isp_amb calculations [Pa].
        npts : int, default 15
            Number of axial stations to precompute along the contour.

        Returns
        -------
        Aerothermodynamics
            Initialized instance with populated design-point performance and all fields
            necessary for :meth:`compute_aerothermodynamics`.

        Notes
        -----
        - The CEA call assumes perfect expansion for ``Isp_vac`` and uses standard
        thrust-coefficient relations to compute ``Isp_amb`` and ``Isp_SL``.
        - Areas and chamber volume are derived from ``c_star``, mass flow, ``L_star``,
        and mixture density at chamber conditions.
        """

        fus = []
        oxs = []
        for prop, frac in zip(fu.propellants, fu.fractions):
            fus.append(cea.Fuel(prop, wt=frac, temp=T_fu_in))

        for prop, frac in zip(ox.propellants, ox.fractions):
            oxs.append(cea.Oxidizer(prop, wt=frac, temp=T_ox_in))
        
        rp = cea.RocketProblem(o_f = MR, 
                               pressure=p_c*0.000145038, # Convert Pa to psi
                               materials=[*oxs, *fus], 
                               pip=p_c/p_e) 
        R = rp.run()
        
        #names = cea.ThermoInterface.keys()
        #print(names)
        
        eps = float(getattr(R, "ae"))
        c_star = float(getattr(R, "cstar"))           
        Isp_ideal_amb = float(getattr(R, "isp"))      # perfectly expanded Isp
        Isp_vac = float(getattr(R, "ivac"))           
        rho_c = float(getattr(R, "c_rho"))
        T_c = float(getattr(R, "c_t"))
        T_t = float(getattr(R, "t_t"))
        p_t = float(getattr(R, "t_p"))

        g = 9.81
        p_SL = 1.01325e5 # sea level pressure
        mdot = F/(Isp_vac*g)
        mdot_fu = mdot/(1+MR)
        mdot_ox = mdot - mdot_fu

        A_t = c_star*mdot/p_c
        r_t = np.sqrt(A_t/np.pi)
        t_stay = L_star*A_t*rho_c/mdot
        V_c = mdot*t_stay/rho_c

        A_e = A_t*eps
        r_e = np.sqrt(A_e/np.pi)

        # calculate ambient Isp
        CF_vac = Isp_vac * g / c_star
        CF_amb = CF_vac - (p_amb / p_c) * (A_e/A_t)
        Isp_amb = CF_amb*c_star/g

        CF_SL = CF_vac - (p_SL / p_c) * (A_e/A_t)
        Isp_SL = CF_SL*c_star/g

        optimum = dict(
            fu=fu, ox=ox, MR=MR, p_c=p_c, p_t=p_t, T_c=T_c, T_t=T_t, F=F, eps=eps, L_star=L_star, c_star=c_star,
            p_amb=p_amb, Isp_ideal_amb=Isp_ideal_amb, Isp_vac=Isp_vac, Isp_amb=Isp_amb,
            Isp_SL=Isp_SL, CF_vac=CF_vac, CF_amb=CF_amb, CF_SL=CF_SL, mdot=mdot,
            mdot_fu=mdot_fu, mdot_ox=mdot_ox, t_stay=t_stay, A_t=A_t, A_e=A_e,
            r_t=r_t, r_e=r_e, V_c=V_c, T_fu_in=T_fu_in, T_ox_in=T_ox_in, npts=npts,
        )

        return cls(optimum)

    def compute_aerothermodynamics(
        self,
        contour,
        Nt: int = 64,
        *,
        T_low: float = 500.0,
        T_hi: float | None = None,
    ):
        """Precompute equilibrium states and common-grid TP property fields.

        The RocketProblem solution defines the equilibrium nozzle trajectory.
        The TP fields are evaluated at the local equilibrium pressure using one
        shared temperature axis from T_low to T_hi.

        Notes
        -----
        - ``*_eq`` arrays hold RocketProblem equilibrium values along x.
        - ``*_map`` arrays hold TPProblem values on a common x-T grid.
        - Temperatures below T_low are clamped during property lookup.
        """

        Nt = int(Nt)
        if Nt < 2:
            raise ValueError("Nt must be at least 2.")

        T_low = float(T_low)
        if not np.isfinite(T_low):
            raise ValueError("T_low must be finite.")

        # ------------------------------------------------------------------
        # Axial grid
        # ------------------------------------------------------------------
        self.x_nodes = np.linspace(contour.xs[0], contour.xs[-1], self.npts)

        A_t = float(contour.A_t)
        eps_nodes = np.asarray(
            [contour.A(x) / A_t for x in self.x_nodes],
            dtype=float,
        )

        Nx = len(self.x_nodes)
        self.Nt = Nt

        # ------------------------------------------------------------------
        # Equilibrium path arrays: RocketProblem output along x.
        # ------------------------------------------------------------------
        self.M_eq = np.full(Nx, np.nan)
        self.T_eq = np.full(Nx, np.nan)
        self.p_eq = np.full(Nx, np.nan)
        self.rho_eq = np.full(Nx, np.nan)
        self.cp_eq = np.full(Nx, np.nan)
        self.gamma_eq = np.full(Nx, np.nan)
        self.h_eq = np.full(Nx, np.nan)
        self.a_eq = np.full(Nx, np.nan)
        self.mu_eq = np.full(Nx, np.nan)
        self.k_eq = np.full(Nx, np.nan)
        self.Pr_eq = np.full(Nx, np.nan)
        self.mw_eq = np.full(Nx, np.nan)

        self.X_map = []

        fus, oxs = self._make_cea_reactants()

        # ------------------------------------------------------------------
        # Pass 1: equilibrium RocketProblem path.
        # ------------------------------------------------------------------
        for i, (x, eps) in enumerate(zip(self.x_nodes, eps_nodes)):
            progress = math.ceil((i + 1) / Nx * 100)
            print(
                f"\rPrecomputing equilibrium path: {progress}%",
                end="",
                flush=True,
            )

            if x < 0.0:
                rp = cea.RocketProblem(
                    o_f=self.MR,
                    pressure=self.p_c * 0.000145038,  # Pa -> psi
                    pressure_units="psi",
                    materials=[*fus, *oxs],
                    sub=eps,
                )
                R = rp.run()
                self.X_map.append(getattr(R, "prod_c", {}))
            else:
                rp = cea.RocketProblem(
                    o_f=self.MR,
                    pressure=self.p_c * 0.000145038,  # Pa -> psi
                    pressure_units="psi",
                    materials=[*fus, *oxs],
                    sup=eps,
                )
                R = rp.run()
                self.X_map.append(getattr(R, "prod_e", {}))

            self.M_eq[i] = float(getattr(R, "mach", np.nan))
            self.T_eq[i] = float(getattr(R, "t", np.nan))
            self.p_eq[i] = float(getattr(R, "p", np.nan))
            self.rho_eq[i] = float(getattr(R, "rho", np.nan))
            self.cp_eq[i] = float(getattr(R, "cp", np.nan))
            self.gamma_eq[i] = float(getattr(R, "gamma", np.nan))
            self.h_eq[i] = float(getattr(R, "h", np.nan))
            self.a_eq[i] = float(getattr(R, "son", np.nan))
            self.mu_eq[i] = float(getattr(R, "visc", np.nan))
            self.k_eq[i] = float(getattr(R, "cond", np.nan))
            self.Pr_eq[i] = float(getattr(R, "pran", np.nan))
            self.mw_eq[i] = float(getattr(R, "mw", np.nan))

        print("\rPrecomputing equilibrium path: 100%")

        # ------------------------------------------------------------------
        # Shared temperature grid for every axial station.
        # ------------------------------------------------------------------
        max_T_eq = float(np.nanmax(self.T_eq))

        if T_hi is None:
            T_hi = max_T_eq
        else:
            T_hi = float(T_hi)

            if T_hi < max_T_eq:
                warnings.warn(
                    f"Requested T_hi={T_hi:.1f} K is below the maximum "
                    f"equilibrium temperature {max_T_eq:.1f} K. "
                    "Using the maximum equilibrium temperature instead."
                )
                T_hi = max_T_eq

        if T_hi < T_low:
            raise ValueError(
                f"T_hi={T_hi:.1f} K must be greater than T_low={T_low:.1f} K."
            )

        self.T_low = T_low
        self.T_hi = T_hi

        # Ascending grid is convenient for np.interp.
        self.T_values = np.linspace(T_low, T_hi, Nt)

        shape = (Nx, Nt)

        # Retained as a 2-D array for plotting and compatibility.
        self.T_grid = np.tile(self.T_values, (Nx, 1))

        # ------------------------------------------------------------------
        # TP fields. Every value in these maps comes from TPProblem.
        # ------------------------------------------------------------------
        self.M_map = np.full(shape, np.nan)
        self.T_map = np.full(shape, np.nan)
        self.p_map = np.full(shape, np.nan)
        self.rho_map = np.full(shape, np.nan)
        self.cp_map = np.full(shape, np.nan)
        self.gamma_map = np.full(shape, np.nan)
        self.h_map = np.full(shape, np.nan)
        self.a_map = np.full(shape, np.nan)
        self.mu_map = np.full(shape, np.nan)
        self.k_map = np.full(shape, np.nan)
        self.Pr_map = np.full(shape, np.nan)
        self.mw_map = np.full(shape, np.nan)

        # ------------------------------------------------------------------
        # Pass 2: TP field at p_eq(x).
        # ------------------------------------------------------------------
        for i, p_bar in enumerate(self.p_eq):
            progress = math.ceil((i + 1) / Nx * 100)
            print(
                f"\rPrecomputing TP property field: {progress}%",
                end="",
                flush=True,
            )

            p_pa = p_bar * 1e5

            for j, Tj in enumerate(self.T_values):
                tp = cea.TPProblem(
                    pressure=p_pa * 0.000145038,  # Pa -> psi
                    pressure_units="psi",
                    temperature=float(Tj),
                    temperature_units="k",
                    materials=[*fus, *oxs],
                    o_f=self.MR,
                )
                Rt = tp.run()

                # TPProblem has no meaningful nozzle-flow Mach number.
                self.T_map[i, j] = float(getattr(Rt, "t", np.nan))
                self.p_map[i, j] = float(getattr(Rt, "p", np.nan))
                self.rho_map[i, j] = float(getattr(Rt, "rho", np.nan))
                self.cp_map[i, j] = float(getattr(Rt, "cp", np.nan))
                self.gamma_map[i, j] = float(getattr(Rt, "gamma", np.nan))
                self.h_map[i, j] = float(getattr(Rt, "h", np.nan))
                self.a_map[i, j] = float(getattr(Rt, "son", np.nan))
                self.mu_map[i, j] = float(getattr(Rt, "visc", np.nan))

                # CEA_Wrap TPProblem conductivity is ten times the RocketProblem
                # conductivity convention. Normalize it here to W/(m K).
                self.k_map[i, j] = float(getattr(Rt, "cond", np.nan)) / 10.0

                self.Pr_map[i, j] = float(getattr(Rt, "pran", np.nan))
                self.mw_map[i, j] = float(getattr(Rt, "mw", np.nan))

        print("\rPrecomputing TP property field: 100%")

    def _bilinear_map(
        self,
        prop_map: np.ndarray,
        x: float,
        T: float,
    ) -> float | None:
        """Interpolate a TP property field in x and T.

        Temperatures below T_low clamp to the value at T_low.
        Temperatures above T_hi return None so a live TP calculation can be used.
        """

        if not (self.x_nodes[0] <= x <= self.x_nodes[-1]):
            return None

        T = float(T)

        if T < self.T_low:
            T = self.T_low

        if T > self.T_hi:
            return None

        i1 = bisect_left(self.x_nodes, x)

        if i1 == 0:
            i0 = 0
            wx = 0.0
        elif i1 >= len(self.x_nodes):
            i0 = len(self.x_nodes) - 1
            wx = 0.0
        else:
            i0 = i1 - 1
            x0 = self.x_nodes[i0]
            x1 = self.x_nodes[i1]
            wx = (x - x0) / (x1 - x0)

        value_0 = float(np.interp(T, self.T_values, prop_map[i0, :]))
        value_1 = float(np.interp(T, self.T_values, prop_map[i1, :]))

        return (1.0 - wx) * value_0 + wx * value_1

    # ---------- Getter functionality ----------
    def _interp_scalar(self, x: float, xs: np.ndarray, ys: np.ndarray) -> float:
        """Linear interpolation with endpoint clamping."""
        if x <= xs[0]:  return float(ys[0])
        if x >= xs[-1]: return float(ys[-1])
        i = bisect_left(xs, x)
        x0, x1 = xs[i-1], xs[i]
        y0, y1 = ys[i-1], ys[i]
        w = (x - x0) / (x1 - x0)
        return float(y0 * (1.0 - w) + y1 * w)

    def _interp_X_dict(self, x: float) -> dict[str, float]:
        """
        Interpolate mole-fraction dict at position x from self.X_map (list[dict]).
        Uses the union of species at the two bracketing nodes, fills missing with 0,
        clips negatives, and renormalizes to sum = 1.
        """
        xs = self.x_nodes
        Xlist = self.X_map
        if not len(xs) or not len(Xlist):
            return {}

        # Clamp OOB to nearest endpoint
        if x <= xs[0]:
            d = Xlist[0] or {}
            S = sum(d.values()) or 1.0
            return {k: max(0.0, v)/S for k, v in d.items()}
        if x >= xs[-1]:
            d = Xlist[-1] or {}
            S = sum(d.values()) or 1.0
            return {k: max(0.0, v)/S for k, v in d.items()}

        # Bracket + blend
        i = bisect_left(xs, x)
        x0, x1 = xs[i-1], xs[i]
        w = (x - x0) / (x1 - x0)
        d0 = Xlist[i-1] or {}
        d1 = Xlist[i]   or {}
        species = set(d0) | set(d1)

        out = {}
        for s in species:
            v = (1.0 - w) * d0.get(s, 0.0) + w * d1.get(s, 0.0)
            out[s] = 0.0 if v < 0.0 else v

        # Renormalize
        S = sum(out.values())
        if S > 0.0:
            for s in out:
                out[s] /= S
        return out

    def _evaluate_tp(self, *, T: float, p: float):
        """Run a live TP calculation at temperature T [K] and pressure p [Pa]."""

        fus, oxs = self._make_cea_reactants()

        tp = cea.TPProblem(
            pressure=p * 0.000145038,  # Pa -> psi
            pressure_units="psi",
            temperature=float(T),
            temperature_units="k",
            materials=[*fus, *oxs],
            o_f=self.MR,
        )

        warnings.warn("Computation outside precomputed range")
        return tp.run()

    def _interp_eq_column(self, map2d: np.ndarray, x: float) -> float:
        """Interpolate along x using the equilibrium column (col 0)."""
        return self._interp_scalar(x, self.x_nodes, map2d[:, 0])

    def _evaluate_hp(self, *, h_target_Jkg: float, p: float):
        # STILL UNDER CONSTRUCTION
        """
        Evaluate an equilibrium HP state at local pressure p [Pa] and mixture specific enthalpy h [J/kg].
        Uses the 'chemical_representation + hf' override path.
        Strategy: assign zero enthalpy to all reactants EXCEPT one 'adjuster' fuel species,
        whose molar enthalpy is chosen so that the mixture enthalpy equals h_target_Jkg.

        Requirements:
          - self.chemrep_map must provide exploded formulas for any species we override.
          - We'll pick the FIRST fuel in self.fu as the 'adjuster'.
        """
        if not self.fu:
            raise RuntimeError("No fuel materials found.")
        adjuster = self.fu[0]
        adj_name = getattr(adjuster, "name", None) or getattr(adjuster, "NAME", None)
        if adj_name is None:
            # CEA_Wrap objects typically have .name
            raise RuntimeError("Cannot determine name of first fuel for enthalpy adjustment.")

        # ---- build overall reactant mass fractions (normalized) ----
        # Group-internal weights:
        fu_w = np.array([getattr(m, "wt", 0.0) for m in self.fu], dtype=float)
        ox_w = np.array([getattr(m, "wt", 0.0) for m in self.ox], dtype=float)
        if fu_w.sum() <= 0 or ox_w.sum() <= 0:
            raise RuntimeError("Fuel/oxidizer component weights must be positive.")

        fu_w /= fu_w.sum()
        ox_w /= ox_w.sum()

        # Overall stream mass fractions given MR = m_ox / m_fu:
        w_fu_stream = 1.0 / (1.0 + self.MR)
        w_ox_stream = self.MR / (1.0 + self.MR)
        w_overall = []
        names = []
        is_fuel_flags = []
        for (m, wf) in zip(self.fu, fu_w):
            nm = getattr(m, "name", None) or getattr(m, "NAME", None)
            names.append(nm); is_fuel_flags.append(True)
            w_overall.append(w_fu_stream * wf)
        for (m, wo) in zip(self.ox, ox_w):
            nm = getattr(m, "name", None) or getattr(m, "NAME", None)
            names.append(nm); is_fuel_flags.append(False)
            w_overall.append(w_ox_stream * wo)
        w_overall = np.asarray(w_overall, dtype=float)

        # ---- compute adjuster molar mass from exploded formula ----
        if adj_name not in self.chemrep_map:
            raise KeyError(f"Missing chemical_representation for adjuster '{adj_name}'. "
                           f"Add to self.chemrep_map.")
        mw_adjuster = _mw_from_exploded(self.chemrep_map[adj_name])  # kg/mol

        # Mass fraction of the adjuster in the overall mixture:
        try:
            idx_adj = names.index(adj_name)  # first occurrence among fuels
        except ValueError:
            raise RuntimeError("Adjuster fuel name not found in assembled reactant list.")
        w_adj = w_overall[idx_adj]
        if w_adj <= 0.0:
            raise RuntimeError("Adjuster species has zero mass fraction; cannot carry enthalpy.")

        # ---- choose assigned molar enthalpies H_i (J/mol) ----
        # All others = 0; adjuster gets H_adj so that h_mix = w_adj * H_adj / mw_adj:
        H_adj = h_target_Jkg * mw_adjuster / w_adj  # J/mol
        assigned_H = np.zeros_like(w_overall)
        assigned_H[idx_adj] = H_adj

        # ---- construct CEA_Wrap materials with chemical_representation + hf ----
        # We preserve inlet temperatures from self.fu/self.ox objects.
        materials = []
        # fuels first, keeping original order
        for (m, wf) in zip(self.fu, fu_w):
            nm = getattr(m, "name", None) or getattr(m, "NAME", None)
            Tm = getattr(m, "temp", 298.15)
            if nm not in self.chemrep_map:
                raise KeyError(f"Missing chemical_representation for fuel '{nm}'.")
            hf = H_adj if (nm == adj_name) else 0.0
            materials.append(cea.Fuel(
                name=nm, wt=float(wf*100.0), temp=float(Tm),
                chemical_representation=self.chemrep_map[nm], hf=float(hf)
            ))
        # oxidizers next
        for (m, wo) in zip(self.ox, ox_w):
            nm = getattr(m, "name", None) or getattr(m, "NAME", None)
            Tm = getattr(m, "temp", 298.15)
            if nm not in self.chemrep_map:
                raise KeyError(f"Missing chemical_representation for oxidizer '{nm}'.")
            materials.append(cea.Oxidizer(
                name=nm, wt=float(wo*100.0), temp=float(Tm),
                chemical_representation=self.chemrep_map[nm], hf=0.0
            ))

        # ---- run HP at local pressure ----
        hp = cea.HPProblem(
            pressure=p * 0.000145038,   # Pa -> psi
            materials=materials,
            o_f=self.MR
        )
        return hp.run()

    def _make_cea_reactants(self):
        """Create fresh CEA_Wrap reactant objects from the stored mixtures."""
        fus = [
            cea.Fuel(prop, wt=frac, temp=self.T_fu_in)
            for prop, frac in zip(self.fu.propellants, self.fu.fractions)
        ]
        oxs = [
            cea.Oxidizer(prop, wt=frac, temp=self.T_ox_in)
            for prop, frac in zip(self.ox.propellants, self.ox.fractions)
        ]
        return fus, oxs

    def _get_prop(
        self,
        *,
        x: float,
        T: float | None,
        h: float | None,
        eq_values: np.ndarray,
        prop_map: np.ndarray,
        res_attr: str,
        live_scale: float = 1.0,
    ) -> float:
        if T is not None and h is not None:
            raise ValueError("Provide only one of T or h.")

        # Equilibrium RocketProblem trajectory.
        if T is None and h is None:
            return self._interp_scalar(x, self.x_nodes, eq_values)

        p_pa = self.get_p(x)

        if T is not None:
            T = max(float(T), self.T_low)

            value = self._bilinear_map(prop_map, x, T)
            if value is not None:
                return value

            R = self._evaluate_tp(T=T, p=p_pa)
            return live_scale * float(getattr(R, res_attr, np.nan))

        R = self._evaluate_hp(h_target_Jkg=h, p=p_pa)
        return live_scale * float(getattr(R, res_attr, np.nan))


    # ----------- public getters (maps + composition) -----------
    def get_T(self, x: float, T: float | None = None, h: float | None = None) -> float:
        return self._get_prop(
            x=x,
            T=T,
            h=h,
            eq_values=self.T_eq,
            prop_map=self.T_map,
            res_attr="t",
        )


    def get_p(self, x: float, T: float | None = None, h: float | None = None) -> float:
        # Local static pressure is always the RocketProblem equilibrium pressure.
        return self._interp_scalar(x, self.x_nodes, self.p_eq) * 1e5


    def get_M(self, x: float, T: float | None = None, h: float | None = None) -> float:
        if T is not None or h is not None:
            raise ValueError(
                "Mach number is only defined for the equilibrium nozzle-flow path."
            )

        return self._interp_scalar(x, self.x_nodes, self.M_eq)


    def get_rho(self, x: float, T: float | None = None, h: float | None = None) -> float:
        return self._get_prop(
            x=x,
            T=T,
            h=h,
            eq_values=self.rho_eq,
            prop_map=self.rho_map,
            res_attr="rho",
        )


    def get_cp(self, x: float, T: float | None = None, h: float | None = None) -> float:
        return self._get_prop(
            x=x,
            T=T,
            h=h,
            eq_values=self.cp_eq,
            prop_map=self.cp_map,
            res_attr="cp",
        ) * 1e3


    def get_gamma(self, x: float, T: float | None = None, h: float | None = None) -> float:
        return self._get_prop(
            x=x,
            T=T,
            h=h,
            eq_values=self.gamma_eq,
            prop_map=self.gamma_map,
            res_attr="gamma",
        )


    def get_h(self, x: float, T: float | None = None, h: float | None = None) -> float:
        return self._get_prop(
            x=x,
            T=T,
            h=h,
            eq_values=self.h_eq,
            prop_map=self.h_map,
            res_attr="h",
        ) * 1e3


    def get_H(self, x: float, T: float | None = None, h: float | None = None) -> float:
        return self.get_h(x=x, T=T, h=h)


    def get_a(self, x: float, T: float | None = None, h: float | None = None) -> float:
        return self._get_prop(
            x=x,
            T=T,
            h=h,
            eq_values=self.a_eq,
            prop_map=self.a_map,
            res_attr="son",
        )


    def get_mu(self, x: float, T: float | None = None, h: float | None = None) -> float:
        return self._get_prop(
            x=x,
            T=T,
            h=h,
            eq_values=self.mu_eq,
            prop_map=self.mu_map,
            res_attr="visc",
        )


    def get_k(self, x: float, T: float | None = None, h: float | None = None) -> float:
        return self._get_prop(
            x=x,
            T=T,
            h=h,
            eq_values=self.k_eq,
            prop_map=self.k_map,
            res_attr="cond",
            live_scale=0.1,
        )


    def get_Pr(self, x: float, T: float | None = None, h: float | None = None) -> float:
        return self._get_prop(
            x=x,
            T=T,
            h=h,
            eq_values=self.Pr_eq,
            prop_map=self.Pr_map,
            res_attr="pran",
        )


    def get_molecular_weight(
        self,
        x: float,
        T: float | None = None,
        h: float | None = None,
    ) -> float:
        return self._get_prop(
            x=x,
            T=T,
            h=h,
            eq_values=self.mw_eq,
            prop_map=self.mw_map,
            res_attr="mw",
        )

    def get_X(self, x: float) -> dict[str, float]:
        """Interpolated mole-fraction dict at position x."""
        return self._interp_X_dict(x)




# Map "public" property keys -> (CEA result attribute, unit_scale)
# unit_scale applies only when we are extracting from a TP run (not from maps).
_PROP_SPEC = {
    "T":     ("t",    1.0),      # K
    "rho":   ("rho",  1.0),      # kg/m^3
    "cp":    ("cp",   1e3),      # kJ/kg-K -> J/kg-K
    "gamma": ("gamma",1.0),
    "h":     ("h",    1e3),      # kJ/kg -> J/kg
    "H":     ("h",    1e3),      # alias of h
    "a":     ("son",  1.0),      # m/s
    "mu":    ("visc", 1.0),      # Pa·s
    "k":     ("cond", 1.0),      # W/m-K
    "Pr":    ("pran", 1.0),
    # "p" is handled specially (we take it from the precomputed map, in Pa)
}

# --- minimal periodic table for MW (kg/mol). Extend as needed. ---
_PERIODIC = {
    "H": 1.00794e-3, "C": 12.0107e-3, "N": 14.0067e-3, "O": 15.9994e-3,
    "F": 18.9984e-3, "Cl": 35.453e-3, "Ar": 39.948e-3, "He": 4.0026e-3,
    "Ne": 20.1797e-3, "S": 32.065e-3, "Si": 28.0855e-3, "B": 10.811e-3,
}
def _mw_from_exploded(chemrep: str) -> float:
    """
    Compute molar mass [kg/mol] from an exploded formula like 'C 2 H 6 O 1'.
    Supports one- or two-letter symbols; coefficients must be integers.
    """
    toks = chemrep.split()
    if len(toks) % 2 != 0:
        raise ValueError(f"Bad chemical_representation: {chemrep!r}")
    i = 0; mw = 0.0
    while i < len(toks):
        sym = toks[i]; n = int(float(toks[i+1])); i += 2
        if sym not in _PERIODIC:
            raise KeyError(f"Element {sym!r} not in periodic table; extend _PERIODIC.")
        mw += _PERIODIC[sym] * n
    return mw


    from __future__ import annotations
