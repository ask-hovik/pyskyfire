"""Friction on the outside of a spinning impeller: its disk friction loss.

A shrouded impeller's two outer faces turn in the fluid filling the gaps
between them and the casing. The shear on them costs shaft power and heats the
fluid. The flow through the blades never sees any of it, so it is charged here,
not in the hydraulic efficiency.

The correlation is Gülich's Eqs. (T3.6.2) and (T3.6.3), with his roughness
factor (T3.6.6). It holds with no leakage through the side gaps: the leakage
factor f_L of (T3.6.7) is taken as one, which matches the no-leakage state of
the pump model until a seal model supplies that flow.

References
----------
Gülich, J. F., *Centrifugal Pumps*, 4th ed., Springer, 2020. Table 3.6.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import cos, log10, pi, radians

__all__ = ["DiskFriction", "disk_friction_power"]


def disk_friction_power(*, rho, mu, n, radius, axial_gap, inner_radius=0.0,
                        roughness=0.0, cone_angle_deg=0.0):
    """Friction power on one face of a rotating disk [W], Gülich (T3.6.2).

    Parameters
    ----------
    rho, mu : float
        Density [kg/m³] and dynamic viscosity [Pa s] of the fluid in the gap.
    n : float
        Shaft speed [rpm].
    radius : float
        Outer radius of the face, normally r2 [m].
    axial_gap : float
        Axial distance between the face and the casing wall, s_ax [m].
    inner_radius : float, optional
        Radius the wetted face starts at, R_n [m]: the hub on the back, the
        seal on the front.
    roughness : float, optional
        Equivalent sand roughness of the face [m].
    cone_angle_deg : float, optional
        Inclination of the face from the radial plane, delta. Below 45 degrees.
    """
    for label, value in (("rho", rho), ("mu", mu), ("n", n), ("radius", radius),
                         ("axial_gap", axial_gap)):
        if not value > 0:
            raise ValueError(f"{label} must be positive")
    if not 0 <= inner_radius < radius:
        raise ValueError("inner_radius must lie in [0, radius)")
    if not roughness >= 0:
        raise ValueError("roughness must be nonnegative")
    if not 0 <= cone_angle_deg < 45:
        raise ValueError("The disk correlation holds for cone angles below 45 degrees")
    omega = 2*pi*n/60
    reynolds = omega*radius**2*rho/mu
    if reynolds <= 10:
        raise ValueError("The disk correlation holds for Re > 10")
    # (T3.6.6): the rotating face's own roughness; the casing's is not counted.
    smooth = log10(12.5/reynolds)
    rough = 1.0 if roughness == 0 else (smooth/log10(0.2*roughness/radius+12.5/reynolds))**2.15
    gap = axial_gap/radius
    k = pi*radius/(2*reynolds*axial_gap) + 0.02/reynolds**0.2*(1+gap)/(1+gap/2)*rough
    return k/cos(radians(cone_angle_deg))*rho*omega**3*radius**5*(1-(inner_radius/radius)**5)


@dataclass(frozen=True)
class DiskFriction:
    """The disk friction of one impeller, all its wetted faces together.

    ``inner_radii`` holds one entry per face: a shrouded impeller has two, the
    front shroud starting at its eye seal and the back at the hub. An open
    impeller has only the back one.
    """
    radius: float
    axial_gap: float
    inner_radii: tuple[float, ...] = (0.0, 0.0)
    roughness: float = 0.0

    def __post_init__(self):
        object.__setattr__(self, "inner_radii", tuple(float(r) for r in self.inner_radii))
        if not self.inner_radii:
            raise ValueError("A disk has at least one wetted face")

    def power(self, *, rho, mu, n):
        """Friction power of every face [W] at one fluid state and speed."""
        return sum(disk_friction_power(rho=rho, mu=mu, n=n, radius=self.radius,
                                       axial_gap=self.axial_gap, inner_radius=inner,
                                       roughness=self.roughness)
                   for inner in self.inner_radii)
