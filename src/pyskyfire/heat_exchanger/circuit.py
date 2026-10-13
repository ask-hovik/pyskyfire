"""Reusable channel banks, finite-volume thermal links, and prescribed gas."""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from pyskyfire.ducts.path import Path, PathPlacement, Straight, Continuous
from pyskyfire.regen.cross_section import SectionProfiles
from pyskyfire.common.physics import T_aw, h_gas_bartz_enthalpy_driven

__all__ = [
    "HeatExchangerCircuit",
    "ThermalLink",
    "FlowConnection",
    "PrescribedGas",
    "ChamberGas",
    "Inlet",
]


def profile(value, s, name):
    values = [value(float(x)) for x in s] if callable(value) else value
    result = np.broadcast_to(np.asarray(values, float), np.shape(s)).copy()
    invalid = (
        np.isnan(result).any() or np.isneginf(result).any()
        if name == "radius"
        else not np.isfinite(result).all()
    )
    if invalid:
        raise ValueError(f"{name} must be finite at every station")
    return result


@dataclass(frozen=True)
class Inlet:
    """Bank inlet: temperature [K], static pressure [Pa], TOTAL mass flow [kg/s]."""

    T: float
    p: float
    mass_flow: float

    def __post_init__(self):
        if (
            not np.isfinite([self.T, self.p, self.mass_flow]).all()
            or min(self.T, self.p, self.mass_flow) <= 0
        ):
            raise ValueError(
                "inlet temperature, pressure and mass_flow must be positive and finite"
            )


class HeatExchangerCircuit:
    """Identical parallel copies of one path, each receiving equal mass flow.

    ``dimensions`` is a mapping accepted by SectionProfiles.from_dimensions;
    each value may be a scalar, station array, or function of path distance s.
    Wall.thickness(s) supplies the floor thickness (and rounded tube thickness).
    Cross-section roof and sidewall settings retain their own definitions.
    Alternatively supply a complete ``section_profiles`` object.

    Path stations anchor the passage area centroid. The local section +u axis
    points floor-to-roof. One thermal solution represents all identical copies;
    prescribed boundaries are checked for equality across placements.

    The default internal coefficient is Gnielinski with bulk properties,
    corrected for a heated fluid by ``(T_b/T_w)^0.45``. ``property_correction=False``
    leaves the bulk-property value uncorrected. Neither applies to a supplied
    ``heat_transfer`` callback.
    """

    def __init__(
        self,
        *,
        name,
        path,
        cross_section,
        wall,
        transport,
        dimensions=None,
        section_profiles=None,
        placements=None,
        roughness=0.0,
        heat_transfer=None,
        laminar_nusselt=3.66,
        property_correction=True,
    ):
        if not isinstance(path, Path):
            raise TypeError("path must be a ducts.Path")
        if not all(isinstance(c, (Straight, Continuous)) for c in path.connections):
            raise ValueError(
                "heat exchangers require Straight or Continuous connections; RoundedBend is unsupported"
            )
        if not name:
            raise ValueError("a circuit requires a nonempty name")
        self.name, self.path, self.cross_section = name, path, cross_section
        self.wall, self.transport = wall, transport
        self.placements = (
            tuple(placements) if placements is not None else (PathPlacement(),)
        )
        if not self.placements or not all(
            isinstance(p, PathPlacement) for p in self.placements
        ):
            raise ValueError("placements must contain at least one PathPlacement")
        self.placed_paths = tuple(p.place(path) for p in self.placements)
        self.frames = tuple(p.frames(path) for p in self.placements)
        self.s = path.s
        if (dimensions is None) == (section_profiles is None):
            raise ValueError("supply dimensions or section_profiles, exclusively")
        if dimensions is not None:
            dims = dict(dimensions)
            if "stations" in dims or "wall_thickness" in dims:
                raise ValueError(
                    "stations come from Path; wall_thickness comes from Wall"
                )
            resolved = {k: profile(v, self.s, k) for k, v in dims.items()}
            # Infinite curvature radius is supported in direct-width mode.
            self.profiles = SectionProfiles.from_dimensions(
                **resolved, stations=self.s, wall_thickness=wall.thickness(self.s)
            )
        else:
            self.profiles = section_profiles
        if len(self.profiles.h) != len(path.stations):
            raise ValueError("section profiles require one entry per path station")
        if not np.allclose(
            self.profiles.t_wall, wall.thickness(self.s), rtol=1e-8, atol=1e-12
        ):
            raise ValueError(
                "section floor thickness must agree with Wall.thickness(s)"
            )
        self.area = cross_section.A_coolant(self.profiles)
        self.diameter = cross_section.Dh_coolant(self.profiles)
        self.outer_area = cross_section.A_exterior(self.profiles)
        if (
            np.any(self.area <= 0)
            or not np.isfinite(self.area + self.diameter + self.outer_area).all()
        ):
            raise ValueError(
                "channel areas and hydraulic diameters must be finite and positive"
            )
        self.inner_perimeters = cross_section.thermal_perimeters(self.profiles)
        self.outer_perimeters = cross_section.thermal_perimeters(
            self.profiles, exterior=True
        )
        self.thicknesses = cross_section.surface_thicknesses(self.profiles)
        if cross_section.enable_fin:
            cross_section.validate_fin(self.profiles)
        self.roughness = profile(roughness, self.s, "roughness")
        if np.any(self.roughness < 0):
            raise ValueError("roughness must be nonnegative")
        self.heat_transfer, self.laminar_nusselt = heat_transfer, laminar_nusselt
        self.property_correction = bool(property_correction)

    @property
    def n_channels(self):
        return len(self.placements)

    @property
    def displaced_volume(self):
        """Sum of swept outer-envelope volumes [m³], assuming no overlap."""
        return float(
            self.n_channels
            * np.sum(
                0.5 * (self.outer_area[:-1] + self.outer_area[1:]) * self.path.lengths
            )
        )

    def section_points(self, placement_index=0, *, exterior=True, samples=12):
        """Closed swept-section vertices [station, perimeter point, xyz]."""
        surfaces = self.cross_section.surfaces(self.profiles, exterior=exterior)
        curves = [
            surfaces[n].evaluate(np.linspace(0, 1, samples, endpoint=False))
            for n in self.cross_section.surface_names
        ]
        # Surfaces are already directed in floor/right/roof/left boundary order.
        local = np.concatenate(curves, axis=1)
        local = local - self.cross_section.centroid(self.profiles)[:, None, :]
        frames = self.frames[placement_index][:, :, :2]
        return self.placed_paths[placement_index].centerline[:, None, :] + np.einsum(
            "nij,nkj->nki", frames, local
        )


@dataclass(frozen=True)
class FlowConnection:
    """An ideal series header: all outlet flow enters the next bank.

    The banks must have the same number of copies and use the same fluid.
    The header has no heat or fitting loss; connecting geometry is not swept.
    """

    upstream: str
    downstream: str


@dataclass(frozen=True)
class ThermalLink:
    """Pair finite volumes on two paths through ONE shared wall.

    Each (i,j) pair couples the segment beginning at station i on a to the
    segment beginning at station j on b. Reverse j for counterflow. All copied
    channels are paired one-to-one. Contact areas [m² per pair per copy] may be
    explicit, otherwise use the smaller exposed surface area. Wall thickness
    is evaluated at the midpoint s on circuit a; it is never counted twice.
    """

    circuit_a: str
    circuit_b: str
    pairs: tuple[tuple[int, int], ...]
    wall: object
    surface_a: str = "floor"
    surface_b: str = "floor"
    areas: object = None


@dataclass(frozen=True)
class PrescribedGas:
    """Fixed surrounding flow: callbacks receive assembly position and wall T.

    ``temperature`` may be a constant or ``temperature(position)``; coefficient
    may be a constant or ``coefficient(position, wall_temperature)`` [W/m²/K].
    """

    temperature: object
    coefficient: object

    def sample(self, position, wall_temperature):
        T = (
            self.temperature(position)
            if callable(self.temperature)
            else self.temperature
        )
        h = (
            self.coefficient(position, wall_temperature)
            if callable(self.coefficient)
            else self.coefficient
        )
        if not np.isfinite([T, h]).all() or min(T, h) <= 0:
            raise ValueError(
                "prescribed gas temperature and coefficient must be finite and positive"
            )
        return float(T), float(h)


class ChamberGas:
    """Adapter to an unchanged thrust-chamber aerothermodynamics solution.

    Assembly coordinates map into chamber coordinates by rotation and position.
    A supplied coefficient [W/m²/K] or callback(position, wall_T) is preferred
    for a specific baffle geometry. The default uses the existing Bartz-style
    reference-temperature correlation as an engineering approximation; it is
    not a resolved flow-around-baffles model.
    """

    def __init__(
        self, chamber, *, position=(0.0, 0.0, 0.0), rotation=None, coefficient=None
    ):
        transform = PathPlacement(position, np.eye(3) if rotation is None else rotation)
        self.position, self.rotation = (
            np.asarray(transform.start_position),
            transform.rotation,
        )
        self.chamber, self.coefficient = chamber, coefficient

    def world_position(self, point):
        return self.rotation @ point + self.position

    def sample(self, point, wall_temperature):
        position = self.world_position(point)
        x = float(position[0])
        contour = self.chamber.contour
        if x < np.min(contour.xs) or x > np.max(contour.xs):
            raise ValueError(
                "heat-exchanger station lies outside the chamber axial domain"
            )
        if np.linalg.norm(position[1:]) > float(contour.r(x)):
            raise ValueError("heat-exchanger station lies outside the chamber radius")
        gas = self.chamber.combustion_transport
        state = gas.get_state(x)
        temperature = T_aw(state.gamma, state.M, state.T, state.Pr)
        if self.coefficient is not None:
            h = (
                self.coefficient(position, wall_temperature)
                if callable(self.coefficient)
                else self.coefficient
            )
        else:
            ref = gas.get_state(x, T=0.5 * (state.T + wall_temperature))
            radius = float(contour.r(x))
            h = h_gas_bartz_enthalpy_driven(
                ref.k,
                2 * radius,
                ref.cp,
                ref.mu,
                gas.mdot,
                np.pi * radius**2,
                state.T,
                ref.T,
            )
        if not np.isfinite([temperature, h]).all() or min(temperature, h) <= 0:
            raise ValueError(
                "chamber boundary returned nonpositive temperature or coefficient"
            )
        return float(temperature), float(h)
