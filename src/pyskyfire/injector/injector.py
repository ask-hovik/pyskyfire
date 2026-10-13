"""An injector: sets of elements placed on the face, and the flow through them.

An :class:`Injector` is a list of :class:`ElementSet` objects. Each set is one
element type repeated at positions the user supplies as ``(r, θ)`` arrays,
however they were found. The injector holds no face outline or keep-out
regions: the positions are the only record of where elements are, and
:meth:`Injector.checks` reports what the injector itself can check.

Hydraulics are per feed. Every hole on one feed shares that feed's plenum
state and discharges to the same face pressure, and the ideal mass flux of an
:class:`~pyskyfire.ducts.Orifice` depends only on those two, so all the holes
on a feed run at one flux. A feed is therefore exactly one orifice with the
summed area and effective area (:meth:`Injector.equivalent_orifice`), and each
hole takes the share of the feed's flow its ``Cd·A`` sets. Holes on different
feeds never combine: a film ring fed at a higher pressure than the main
plenum is solved on its own.

:meth:`Injector.solve` answers either question per feed: the flow a plenum
pressure drives (``pressure`` given), or the plenum pressure a flow needs
(``mass_flow`` given). :meth:`Injector.solve_feed` is the engine network's
direction, plenum state and flow in, face pressure out.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import math

import numpy as np
from scipy.spatial import cKDTree

from pyskyfire.ducts.orifice import Orifice, OrificeResult
from pyskyfire.thermo.eos import as_equation_of_state, fluid_name
from pyskyfire.thermo.state import FluidState

from .element import Element, JetFlow

__all__ = [
    "ElementSet",
    "Injector",
    "InjectorHoles",
    "FeedCondition",
    "FeedResult",
    "SetResult",
    "InjectorResult",
    "InjectorChecks",
]


def _rotation_x(angle):
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


def _diameter(orifice: Orifice) -> float:
    if orifice.hydraulic_diameter is not None:
        return orifice.hydraulic_diameter
    return math.sqrt(4.0 * orifice.area / math.pi)


class ElementSet:
    """One element type repeated at supplied face positions.

    Parameters
    ----------
    element : Element
        The element every position carries.
    r, theta : array_like
        Radius [m] and angle [rad] of each element centre on the face. The
        face is the plane x = 0 with the chamber towards +x, and θ is measured
        from +y towards +z.
    orientation : float or array_like, optional
        Angle of each element's own +y axis from the local radial direction,
        towards the local tangential direction [rad]. One value for all, or
        one per element.
    name : str
        Label for results and scenes. Unique within an injector.
    """

    def __init__(self, element, r, theta, *, orientation=0.0, name):
        if not isinstance(element, Element):
            raise TypeError("element must be an injector Element")
        r = np.array(r, dtype=float, ndmin=1)
        theta = np.array(theta, dtype=float, ndmin=1)
        if r.ndim != 1 or r.shape != theta.shape or r.size == 0:
            raise ValueError("r and theta must be 1D arrays of equal, nonzero length")
        if not (np.isfinite(r).all() and np.isfinite(theta).all()) or np.any(r < 0.0):
            raise ValueError("r must be finite and non-negative, theta finite")
        orientation = np.broadcast_to(np.asarray(orientation, dtype=float), r.shape).copy()
        if not np.isfinite(orientation).all():
            raise ValueError("orientation must be finite")
        if not isinstance(name, str) or not name:
            raise ValueError("name must be a nonempty string")
        for array in (r, theta, orientation):
            array.setflags(write=False)
        self.element = element
        self.r = r
        self.theta = theta
        self.orientation = orientation
        self.name = name

    def __repr__(self):
        return f"ElementSet({self.name!r}, {self.count} × {self.element!r})"

    @property
    def count(self) -> int:
        """Number of elements in the set."""
        return int(self.r.size)

    @property
    def centers(self) -> np.ndarray:
        """Element centres on the face, shape ``(count, 3)`` [m]."""
        return np.column_stack((np.zeros(self.count), self.r * np.cos(self.theta),
                                self.r * np.sin(self.theta)))

    def rotations(self) -> np.ndarray:
        """Rotation of each element's frame into the injector's, ``(count, 3, 3)``."""
        return np.array([_rotation_x(t + o) for t, o in zip(self.theta, self.orientation)])


@dataclass(frozen=True)
class InjectorHoles:
    """Every hole of an injector, placed in the injector frame.

    Arrays have one row per hole. ``set_name`` and ``element_index`` say
    which element a hole belongs to; ``feeds`` and ``orifices`` are tuples
    aligned with them. ``footprint`` is the radius of the hole's outline on
    the face [m]: half its diameter, stretched by the inclination.
    """

    set_name: tuple[str, ...]
    element_index: np.ndarray
    feeds: tuple
    orifices: tuple[Orifice, ...]
    exit: np.ndarray
    direction: np.ndarray
    target: np.ndarray
    footprint: np.ndarray

    def __len__(self):
        return len(self.orifices)


@dataclass(frozen=True)
class FeedCondition:
    """The operating point of one feed.

    Give the plenum ``temperature`` [K], the ``fluid`` (anything
    :func:`~pyskyfire.thermo.eos.as_equation_of_state` reads), and either the
    plenum ``pressure`` [Pa], to find the flow it drives, or the
    ``mass_flow`` [kg/s], to find the pressure it needs.
    """

    fluid: object
    temperature: float
    pressure: float | None = None
    mass_flow: float | None = None

    def __post_init__(self):
        if (self.pressure is None) == (self.mass_flow is None):
            raise ValueError("give exactly one of pressure and mass_flow")


@dataclass(frozen=True)
class FeedResult:
    """The flow through every hole on one feed.

    ``orifice`` is the result for the feed's equivalent orifice: its inlet
    is the plenum and its outlet the face. ``set_mass_flows`` gives each
    element set's share [kg/s]. ``fluid`` is the CoolProp name the feed
    ran on, or None when its fluid named none.
    """

    feed: object
    orifice: OrificeResult
    effective_area: float
    hole_count: int
    set_mass_flows: dict
    fluid: str | None = None

    @property
    def mass_flow(self) -> float:
        """Total mass flow of the feed [kg/s]."""
        return self.orifice.mass_flow

    @property
    def plenum_pressure(self) -> float:
        """Plenum (inlet) pressure [Pa]."""
        return self.orifice.inlet_state.p

    @property
    def pressure_drop(self) -> float:
        """Drop from the plenum to the face [Pa]."""
        return self.orifice.pressure_loss

    @property
    def velocity(self) -> float:
        """Ideal jet velocity, the same from every hole on the feed [m/s]."""
        return self.orifice.velocity

    @property
    def mass_flux(self) -> float:
        """Mass flow per unit effective area [kg/(m² s)]."""
        return self.mass_flow / self.effective_area

    def hole_mass_flow(self, orifice: Orifice) -> float:
        """Mass flow through one hole on this feed [kg/s]."""
        return orifice.effective_area * self.mass_flux


@dataclass(frozen=True)
class SetResult:
    """The flow through one element set.

    ``mass_flows`` maps each of the set's feeds to the set's total from it
    [kg/s]. ``diagnostics`` is the element type's own mixing report for one
    element; every element in a set runs identically.
    """

    name: str
    count: int
    mass_flows: dict
    diagnostics: object

    @property
    def mass_flow(self) -> float:
        """Total mass flow through the set [kg/s]."""
        return sum(self.mass_flows.values())

    @property
    def element_mass_flow(self) -> float:
        """Mass flow through one element [kg/s]."""
        return self.mass_flow / self.count


@dataclass(frozen=True)
class InjectorResult:
    """Every feed and set at one operating point, into a face at ``p_face``."""

    p_face: float
    feeds: dict
    sets: dict
    oxidizer: frozenset

    @property
    def mass_flow(self) -> float:
        """Total mass flow [kg/s]."""
        return sum(feed.mass_flow for feed in self.feeds.values())

    @property
    def oxidizer_mass_flow(self) -> float:
        """Mass flow of the oxidizer feeds [kg/s]."""
        return sum(r.mass_flow for f, r in self.feeds.items() if f in self.oxidizer)

    @property
    def fuel_mass_flow(self) -> float:
        """Mass flow of every other feed [kg/s]."""
        return self.mass_flow - self.oxidizer_mass_flow

    @property
    def mixture_ratio(self) -> float:
        """Overall oxidizer over fuel mass flow, NaN without a named oxidizer."""
        if not self.oxidizer or self.fuel_mass_flow <= 0.0:
            return math.nan
        return self.oxidizer_mass_flow / self.fuel_mass_flow


@dataclass(frozen=True)
class InjectorChecks:
    """What the injector can check about its own layout, reported not enforced.

    ``min_web`` is the narrowest ligament between two hole outlines on the
    face [m], between holes in sets ``min_web_sets``; ``overlaps`` counts
    pairs whose outlines intersect. ``outer_radius`` is how far the
    outermost hole outline reaches from the axis, and ``target_radius`` the
    furthest any jet's target lies from it [m]. ``wall_clearance`` is the
    chamber radius less ``outer_radius``: room between the holes and the
    wall. ``target_clearance`` is the chamber radius less ``target_radius``:
    zero where a jet is aimed at the wall on purpose, as a film-cooling jet
    is, and negative where a target lies beyond it. Both are NaN without a
    chamber radius. ``min_length_to_diameter`` is the shortest hole relative
    to its diameter.
    """

    min_web: float
    min_web_sets: tuple[str, str]
    overlaps: int
    outer_radius: float
    target_radius: float
    wall_clearance: float
    target_clearance: float
    min_length_to_diameter: float


class Injector:
    """Element sets on a face, and the hydraulics of the feeds that supply them.

    Parameters
    ----------
    element_sets : sequence of ElementSet
        The sets, with unique names.
    """

    def __init__(self, element_sets):
        sets = tuple(element_sets)
        if not sets:
            raise ValueError("an injector needs at least one element set")
        if not all(isinstance(s, ElementSet) for s in sets):
            raise TypeError("element_sets must contain only ElementSet instances")
        names = [s.name for s in sets]
        if len(set(names)) != len(names):
            raise ValueError("element set names must be unique")
        self.element_sets = sets

    def __repr__(self):
        return f"Injector({', '.join(repr(s) for s in self.element_sets)})"

    @property
    def sets(self) -> dict:
        """Element sets by name."""
        return {s.name: s for s in self.element_sets}

    @property
    def feeds(self) -> tuple:
        """Every feed, in the order the sets first name them."""
        return tuple(dict.fromkeys(f for s in self.element_sets for f in s.element.feeds))

    def _check_feed(self, feed):
        if feed not in self.feeds:
            raise KeyError(f"{feed!r} feeds no hole of this injector")

    def equivalent_orifice(self, feed) -> Orifice:
        """One orifice with the summed area and ``Cd·A`` of every hole on ``feed``."""
        self._check_feed(feed)
        holes = [(s.count, o) for s in self.element_sets for o in s.element.orifices(feed)]
        area = sum(n * o.area for n, o in holes)
        return Orifice(area, cd=sum(n * o.effective_area for n, o in holes) / area)

    def hole_count(self, feed) -> int:
        """Number of holes on ``feed``."""
        self._check_feed(feed)
        return sum(s.count * len(s.element.orifices(feed)) for s in self.element_sets)

    def _feed_result(self, feed, orifice_result, fluid) -> FeedResult:
        equivalent = self.equivalent_orifice(feed)
        flux = orifice_result.mass_flow / equivalent.effective_area
        shares = {s.name: s.count * flux * sum(o.effective_area for o in s.element.orifices(feed))
                  for s in self.element_sets if feed in s.element.feeds}
        return FeedResult(feed, orifice_result, equivalent.effective_area,
                          self.hole_count(feed), shares, fluid_name(fluid))

    def solve_feed(self, feed, *, inlet_state: FluidState, mass_flow: float, fluid) -> FeedResult:
        """The face pressure one feed reaches from its plenum at ``mass_flow``.

        The engine network's direction: see :meth:`Orifice.solve`.
        """
        result = self.equivalent_orifice(feed).solve(
            inlet_state=inlet_state, mass_flow=mass_flow, fluid=fluid)
        return self._feed_result(feed, result, fluid)

    def solve(self, conditions: Mapping, *, p_face: float, oxidizer=()) -> InjectorResult:
        """Every feed at its operating point, discharging into a face at ``p_face``.

        Parameters
        ----------
        conditions : mapping of feed to FeedCondition
            One condition for every feed of the injector. A condition with a
            pressure gives that feed's flow; one with a mass flow gives the
            plenum pressure it needs.
        p_face : float
            Pressure at the injector face [Pa].
        oxidizer : iterable of feeds, optional
            The feeds carrying oxidizer. Mixture ratios need them.
        """
        feeds = self.feeds
        missing = [f for f in feeds if f not in conditions]
        extra = [f for f in conditions if f not in feeds]
        if missing or extra:
            raise ValueError(f"conditions must cover exactly the injector's feeds; "
                             f"missing {missing}, unknown {extra}")
        oxidizer = frozenset(oxidizer)
        if not oxidizer <= set(feeds):
            raise ValueError("oxidizer names a feed the injector does not have")

        results = {}
        for feed in feeds:
            condition = conditions[feed]
            orifice = self.equivalent_orifice(feed)
            if condition.mass_flow is not None:
                result = orifice.solve_inlet_pressure(
                    mass_flow=condition.mass_flow, p_out=p_face,
                    T_in=condition.temperature, fluid=condition.fluid)
            else:
                state = as_equation_of_state(condition.fluid).state_pT(
                    condition.pressure, condition.temperature)
                result = orifice.mass_flow(inlet_state=state, p_out=p_face, fluid=condition.fluid)
            results[feed] = self._feed_result(feed, result, condition.fluid)

        sets = {}
        for s in self.element_sets:
            flows = tuple(JetFlow(results[h.feed].hole_mass_flow(h.orifice),
                                  results[h.feed].velocity,
                                  results[h.feed].orifice.throat_state.rho)
                          for h in s.element.holes())
            sets[s.name] = SetResult(
                s.name, s.count,
                {f: results[f].set_mass_flows[s.name] for f in s.element.feeds},
                s.element.diagnostics(flows, oxidizer))
        return InjectorResult(float(p_face), results, sets, oxidizer)

    def holes(self) -> InjectorHoles:
        """Every hole, placed in the injector frame."""
        names, index, feeds, orifices = [], [], [], []
        exits, directions, targets, footprint = [], [], [], []
        for s in self.element_sets:
            local = s.element.holes()
            for k, (center, rotation) in enumerate(zip(s.centers, s.rotations())):
                for hole in local:
                    names.append(s.name)
                    index.append(k)
                    feeds.append(hole.feed)
                    orifices.append(hole.orifice)
                    exits.append(center + rotation @ hole.exit)
                    direction = rotation @ hole.direction
                    directions.append(direction)
                    targets.append(center + rotation @ hole.target)
                    footprint.append(0.5 * _diameter(hole.orifice) / direction[0])
        return InjectorHoles(tuple(names), np.asarray(index), tuple(feeds), tuple(orifices),
                             np.asarray(exits), np.asarray(directions), np.asarray(targets),
                             np.asarray(footprint))

    def checks(self, *, chamber_radius: float | None = None) -> InjectorChecks:
        """Layout checks: webs between holes, holes and targets against the wall, hole L/D.

        Hole outlines on the face are taken as circles of the
        :attr:`InjectorHoles.footprint` radius, which is conservative for an
        inclined hole's ellipse.
        """
        holes = self.holes()
        face = holes.exit[:, 1:]
        radius = holes.footprint
        min_web, pair, overlaps = math.inf, ("", ""), 0
        if len(holes) > 1:
            k = min(len(holes), 9)
            distances, neighbours = cKDTree(face).query(face, k=k)
            seen = set()
            for i in range(len(holes)):
                for d, j in zip(distances[i, 1:], neighbours[i, 1:]):
                    web = d - radius[i] - radius[j]
                    key = (min(i, j), max(i, j))
                    if web < 0.0 and key not in seen:
                        overlaps += 1
                    seen.add(key)
                    if web < min_web:
                        min_web, pair = web, (holes.set_name[i], holes.set_name[j])
        outer = float(np.max(np.linalg.norm(face, axis=1) + radius))
        target = float(np.max(np.linalg.norm(holes.target[:, 1:], axis=1)))
        wall = math.nan if chamber_radius is None else float(chamber_radius)
        l_over_d = min(o.length / _diameter(o) for o in holes.orifices)
        return InjectorChecks(float(min_web), pair, overlaps, outer, target, wall - outer,
                              wall - target, l_over_d)
