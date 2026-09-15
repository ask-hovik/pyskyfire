"""Steady single-phase heat exchange with conservative segment heat rates.

A circuit solution represents identical parallel copies. Static pressure obeys
Darcy friction plus rho*u*du; energy conserves h + u²/2. Exchanger heat rates
are iterated against segment midpoint states, including counterflow links.
"""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.optimize import brentq

from pyskyfire.common.friction import darcy_friction_factor, curvature_multiplier
from pyskyfire.common.momentum import acceleration_pressure_drop
from pyskyfire.common.physics import h_internal
from pyskyfire.thermo.state import FluidStateResolver
from .circuit import Inlet, ThermalLink, FlowConnection

__all__ = ["HeatExchanger", "HeatExchangerResult", "CircuitResult", "WallResult"]


def midpoint(a):
    a = np.asarray(a)
    return (a[:-1] + a[1:]) / 2


@dataclass
class CircuitResult:
    """One representative channel; heat_rate is TOTAL bank heat per segment."""

    circuit: object
    states: tuple
    mass_flow: float
    velocity: np.ndarray
    heat_rate: np.ndarray
    friction_loss: np.ndarray
    acceleration_drop: np.ndarray
    reynolds: np.ndarray
    darcy_factor: np.ndarray
    heat_transfer_coefficient: np.ndarray

    @property
    def T(self):
        return np.array([s.T for s in self.states])

    @property
    def p(self):
        return np.array([s.p for s in self.states])

    @property
    def h(self):
        return np.array([s.h for s in self.states])

    @property
    def s(self):
        return self.circuit.s

    @property
    def outlet(self):
        return self.states[-1]

    @property
    def total_heat(self):
        return float(np.sum(self.heat_rate))

    @property
    def energy_change(self):
        return self.mass_flow * (
            self.h[-1]
            - self.h[0]
            + 0.5 * (self.velocity[-1] ** 2 - self.velocity[0] ** 2)
        )


@dataclass
class WallResult:
    """One physical wall patch, with temperatures at both fluid interfaces.

    Heat rate [W] is positive from b (or prescribed gas) into a, summed across
    the bank. Area [m²] is per copied channel. Indices identify finite volumes.
    """

    circuit_a: str
    segment_a: int
    surface_a: str
    circuit_b: str | None
    segment_b: int | None
    surface_b: str | None
    temperature_a: float
    temperature_b: float
    heat_rate: float
    area: float


@dataclass
class HeatExchangerResult:
    circuits: dict[str, CircuitResult]
    walls: tuple[WallResult, ...]
    iterations: int
    relative_residual: float
    prescribed_heat: float
    displaced_volume: float

    @property
    def energy_balance_error(self):
        """Sum of all bank energy changes minus prescribed-gas heat [W]."""
        return (
            sum(r.energy_change for r in self.circuits.values()) - self.prescribed_heat
        )


def _property(transport, name, T, p):
    method = getattr(transport, "get_" + name, None) or getattr(transport, name, None)
    if not callable(method):
        raise TypeError(f"transport must provide get_{name}(T,p)")
    value = float(method(T, p))
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"invalid fluid {name}: {value}")
    return value


def _single_phase(resolver, h, p):
    state = resolver.from_hp(h, p)
    if state.is_two_phase:
        raise ValueError(
            "heat exchanger reached a two-phase state; only single-phase flow is supported"
        )
    return state


def _march(circuit, inlet, heat, resolver, friction_model):
    mdot = inlet.mass_flow / circuit.n_channels
    first = _single_phase(resolver, resolver.enthalpy(inlet.T, inlet.p), inlet.p)
    states = [first]
    velocity = [mdot / (first.rho * circuit.area[0])]
    losses, accelerations, re_values, fd_values, coefficients = [], [], [], [], []
    diameters, areas, roughness = (
        midpoint(circuit.diameter),
        midpoint(circuit.area),
        midpoint(circuit.roughness),
    )
    radii = circuit.path.curvature_radius
    for i, length in enumerate(circuit.path.lengths):
        up, u0 = states[-1], velocity[-1]
        h_total = up.h + u0 * u0 / 2 + heat[i] / mdot
        p1, h1 = up.p, h_total - u0 * u0 / 2
        for sweep in range(60):
            if p1 <= 0:
                raise RuntimeError(
                    f"{circuit.name}: pressure fell below zero at segment {i}"
                )
            down = _single_phase(resolver, h1, p1)
            u1 = mdot / (down.rho * circuit.area[i + 1])
            mid = _single_phase(resolver, 0.5 * (up.h + down.h), 0.5 * (up.p + down.p))
            u = mdot / (mid.rho * areas[i])
            mu = _property(circuit.transport, "mu", mid.T, mid.p)
            Re = mid.rho * u * diameters[i] / mu
            fd = darcy_friction_factor(
                Re, diameters[i], roughness[i], model=friction_model
            )
            factor = fd * curvature_multiplier(
                Re, diameters[i], radii[i], roughness=roughness[i], model=friction_model
            )
            loss = factor * length / diameters[i] * mid.rho * u * u / 2
            acceleration = acceleration_pressure_drop(up.rho, u0, down.rho, u1)
            target_p = up.p - loss - acceleration
            target_h = h_total - u1 * u1 / 2
            error = max(
                abs(target_p - p1) / max(up.p, 1),
                abs(target_h - h1) / max(abs(h_total), 1),
            )
            if error < 1e-10:
                break
            p1 = 0.5 * p1 + 0.5 * target_p
            h1 = 0.5 * h1 + 0.5 * target_h
        else:
            raise RuntimeError(
                f"{circuit.name}: pressure/energy iteration failed at segment {i}; reduce step or flow speed"
            )
        k = _property(circuit.transport, "k", mid.T, mid.p)
        Pr = _property(circuit.transport, "Pr", mid.T, mid.p)
        # Gnielinski receives the straight Darcy factor; coil friction is not
        # implicitly reused as an unvalidated heat-transfer enhancement.
        if circuit.heat_transfer is None:
            htc = h_internal(
                Re, Pr, k, diameters[i], fd, laminar_nusselt=circuit.laminar_nusselt
            )
        else:
            htc = float(
                circuit.heat_transfer(
                    Re=Re, Pr=Pr, k=k, Dh=diameters[i], darcy_factor=fd
                )
            )
        if not np.isfinite(htc) or htc <= 0:
            raise ValueError(
                "internal heat-transfer coefficient must be positive and finite"
            )
        states.append(down)
        velocity.append(u1)
        losses.append(loss)
        accelerations.append(acceleration)
        re_values.append(Re)
        fd_values.append(factor)
        coefficients.append(htc)
    return CircuitResult(
        circuit,
        tuple(states),
        inlet.mass_flow,
        np.asarray(velocity),
        np.asarray(heat) * circuit.n_channels,
        np.asarray(losses),
        np.asarray(accelerations),
        np.asarray(re_values),
        np.asarray(fd_values),
        np.asarray(coefficients),
    )


def _conductivity(wall, temperature):
    k = float(wall.material.get_k(temperature))
    if not np.isfinite(k) or k <= 0:
        raise ValueError("wall conductivity must be finite and positive")
    return k


def _wall_exchange(Ta, Tb, Ra, Rb, thickness, area, wall):
    if (
        not np.isfinite([Ra, Rb, thickness, area]).all()
        or min(Ra, Rb, area) <= 0
        or thickness < 0
    ):
        raise ValueError("invalid wall resistance, thickness, or contact area")
    if Ta == Tb:
        return 0.0, Ta, Tb

    # q is signed from b to a. Solve mean wall temperature with k(T).
    def temperatures(Tmean):
        resistance = thickness / (_conductivity(wall, Tmean) * area)
        q = (Tb - Ta) / (Ra + Rb + resistance)
        return q, Ta + q * Ra, Tb - q * Rb

    mean = brentq(
        lambda T: 0.5 * sum(temperatures(T)[1:]) - T,
        min(Ta, Tb),
        max(Ta, Tb),
        xtol=1e-8,
    )
    return temperatures(mean)


class HeatExchanger:
    """An assembly of banks, shared walls, and optional series headers.

    Boundary mappings and inlet mappings use circuit names. For series banks,
    supply an inlet only for the first bank. Thermal links may pair arbitrary
    segment indices; exposed-area budgets prevent overlapping thermal links.
    A prescribed boundary heats whatever exposed area is not assigned a link.
    """

    def __init__(self, circuits, *, thermal_links=(), flow_connections=()):
        circuits = tuple(circuits)
        self.circuits = {c.name: c for c in circuits}
        if not circuits or len(self.circuits) != len(circuits):
            raise ValueError("heat exchanger requires uniquely named circuits")
        self.thermal_links = tuple(thermal_links)
        self.flow_connections = tuple(flow_connections)
        self._parents = {}
        children = set()
        for edge in self.flow_connections:
            if (
                not isinstance(edge, FlowConnection)
                or edge.upstream not in self.circuits
                or edge.downstream not in self.circuits
            ):
                raise ValueError("flow connection names must identify circuits")
            if edge.downstream in self._parents or edge.upstream in children:
                raise ValueError(
                    "series flow permits one inlet and one outlet per bank"
                )
            a, b = self.circuits[edge.upstream], self.circuits[edge.downstream]
            if a.n_channels != b.n_channels or a.transport.fluid != b.transport.fluid:
                raise ValueError("series banks require matching copy counts and fluids")
            self._parents[edge.downstream] = edge.upstream
            children.add(edge.upstream)
        self._order = []
        pending = set(self.circuits)
        while pending:
            ready = [
                name
                for name in self.circuits
                if name in pending and self._parents.get(name) not in pending
            ]
            if not ready:
                raise ValueError("series flow connections contain a cycle")
            self._order.extend(ready)
            pending.difference_update(ready)
        self._contacts, self._remaining = self._build_contacts()

    def _build_contacts(self):
        remaining = {
            (name, surface): midpoint(p) * c.path.lengths
            for name, c in self.circuits.items()
            for surface, p in c.inner_perimeters.items()
        }
        contacts = []
        for link in self.thermal_links:
            if (
                not isinstance(link, ThermalLink)
                or link.circuit_a not in self.circuits
                or link.circuit_b not in self.circuits
                or link.circuit_a == link.circuit_b
            ):
                raise ValueError("a thermal link must identify two different circuits")
            a, b = self.circuits[link.circuit_a], self.circuits[link.circuit_b]
            if a.n_channels != b.n_channels:
                raise ValueError(
                    "thermal links require one-to-one copied channels with equal counts"
                )
            if (a.name, link.surface_a) not in remaining or (
                b.name,
                link.surface_b,
            ) not in remaining:
                raise ValueError("unknown thermal-link surface")
            pairs = tuple(link.pairs)
            areas = (
                None
                if link.areas is None
                else np.broadcast_to(np.asarray(link.areas, float), (len(pairs),))
            )
            for n, (i, j) in enumerate(pairs):
                if (
                    not isinstance(i, (int, np.integer))
                    or not isinstance(j, (int, np.integer))
                    or not 0 <= i < len(a.s) - 1
                    or not 0 <= j < len(b.s) - 1
                ):
                    raise ValueError("thermal-link indices must identify path segments")
                ra, rb = (
                    remaining[(a.name, link.surface_a)],
                    remaining[(b.name, link.surface_b)],
                )
                area = min(ra[i], rb[j]) if areas is None else float(areas[n])
                if (
                    not np.isfinite(area)
                    or area <= 0
                    or area > min(ra[i], rb[j]) + 1e-12
                ):
                    raise ValueError(
                        "thermal-link area exceeds available exposed surface or repeats an occupied patch"
                    )
                ra[i] -= area
                rb[j] -= area
                contacts.append((link, i, j, area))
        return contacts, remaining

    @property
    def displaced_volume(self):
        return sum(c.displaced_volume for c in self.circuits.values())

    def solve(
        self,
        inlets,
        *,
        boundaries=None,
        friction_model="fluids",
        max_iterations=250,
        tolerance=1e-7,
        relaxation=0.5,
    ):
        boundaries = {} if boundaries is None else dict(boundaries)
        if set(inlets) != set(self.circuits) - set(self._parents):
            raise ValueError(
                "supply an inlet for each independent circuit, and none for series-fed circuits"
            )
        if set(boundaries) - set(self.circuits):
            raise ValueError("boundary names must identify circuits")
        if (
            not 0 < relaxation <= 1
            or tolerance <= 0
            or not np.isfinite([relaxation, tolerance]).all()
            or not isinstance(max_iterations, int)
            or max_iterations < 1
        ):
            raise ValueError("invalid iteration settings")
        resolvers = {}
        for name in self._order:
            inlet = inlets.get(name)
            p = inlet.p if inlet is not None else resolvers[self._parents[name]][1]
            resolver = FluidStateResolver(self.circuits[name].transport, p)
            if not resolver.capability.enthalpy_march:
                raise ValueError(
                    f"{name}: an enthalpy-capable single-phase fluid backend is required"
                )
            resolvers[name] = (resolver, p)
        heat = {name: np.zeros(len(c.s) - 1) for name, c in self.circuits.items()}
        previous_walls = ()
        for iteration in range(1, max_iterations + 1):
            results = {}
            for name in self._order:
                c = self.circuits[name]
                if name in self._parents:
                    upstream = results[self._parents[name]]
                    inlet = Inlet(
                        upstream.outlet.T, upstream.outlet.p, upstream.mass_flow
                    )
                    # Ideal header has equal endpoint area, so carries total
                    # enthalpy without an unmodelled acceleration jump.
                    if not np.isclose(
                        upstream.circuit.area[-1], c.area[0], rtol=1e-6, atol=1e-12
                    ):
                        raise ValueError(
                            "series header endpoint areas must match; model a continuous area transition on the path"
                        )
                else:
                    inlet = inlets[name]
                results[name] = _march(
                    c, inlet, heat[name], resolvers[name][0], friction_model
                )
            candidate, walls = self._exchange(results, boundaries, previous_walls)
            scale = max(1.0, max(float(np.max(np.abs(q))) for q in candidate.values()))
            error = (
                max(float(np.max(np.abs(candidate[n] - heat[n]))) for n in heat) / scale
            )
            if walls and previous_walls:
                wall_error = max(
                    max(
                        abs(w.temperature_a - old.temperature_a),
                        abs(w.temperature_b - old.temperature_b),
                    )
                    / max(1.0, abs(w.temperature_a), abs(w.temperature_b))
                    for w, old in zip(walls, previous_walls)
                )
                error = max(error, wall_error)
            elif walls:
                error = max(error, 1.0)
            previous_walls = walls
            if error < tolerance:
                prescribed = sum(w.heat_rate for w in walls if w.circuit_b is None)
                return HeatExchangerResult(
                    results,
                    tuple(walls),
                    iteration,
                    error,
                    prescribed,
                    self.displaced_volume,
                )
            for name in heat:
                heat[name] += relaxation * (candidate[name] - heat[name])
        raise RuntimeError(
            f"heat exchanger did not converge in {max_iterations} iterations (relative coupling residual {error:.3g})"
        )

    def _exchange(self, results, boundaries, previous_walls=()):
        heat = {name: np.zeros(len(c.s) - 1) for name, c in self.circuits.items()}
        temperatures = {name: midpoint(r.T) for name, r in results.items()}
        walls = []
        # Fin enhancement is distributed over active base surfaces. Sidewalls
        # cannot also be active surfaces when symmetric ribs are enabled.
        enhancement = {}
        for name, c in self.circuits.items():
            if c.cross_section.enable_fin:
                h = np.interp(
                    c.s, midpoint(c.s), results[name].heat_transfer_coefficient
                )
                # The fin material uses the previous wall solution, not the
                # cryogenic bulk-fluid temperature. Iterate both heat and wall
                # temperatures to settle temperature-dependent conductivities.
                sums = np.zeros(len(c.s) - 1)
                areas = np.zeros_like(sums)
                for w in previous_walls:
                    i = (
                        w.segment_a
                        if w.circuit_a == name
                        else w.segment_b
                        if w.circuit_b == name
                        else None
                    )
                    if i is not None:
                        sums[i] += 0.5 * (w.temperature_a + w.temperature_b) * w.area
                        areas[i] += w.area
                wall_mid = np.divide(
                    sums, areas, out=midpoint(results[name].T), where=areas > 0
                )
                wall_T = np.interp(c.s, midpoint(c.s), wall_mid)
                k = np.array([_conductivity(c.wall, T) for T in wall_T])
                effective = 1 / (
                    c.cross_section.R_coolant_per_len(c.profiles, h, k) * h
                )
                base = c.cross_section.P_coolant(c.profiles)
                enhancement[name] = midpoint(
                    np.divide(effective, base, out=np.ones_like(base), where=base > 0)
                )
            else:
                enhancement[name] = np.ones(len(c.s) - 1)
        for link, i, j, area in self._contacts:
            a, b = self.circuits[link.circuit_a], self.circuits[link.circuit_b]
            ha = results[a.name].heat_transfer_coefficient[i] * enhancement[a.name][i]
            hb = results[b.name].heat_transfer_coefficient[j] * enhancement[b.name][j]
            thickness = float(link.wall.thickness(midpoint(a.s)[i]))
            q, Ta, Tb = _wall_exchange(
                temperatures[a.name][i],
                temperatures[b.name][j],
                1 / (ha * area),
                1 / (hb * area),
                thickness,
                area,
                link.wall,
            )
            heat[a.name][i] += q
            heat[b.name][j] -= q
            walls.append(
                WallResult(
                    a.name,
                    i,
                    link.surface_a,
                    b.name,
                    j,
                    link.surface_b,
                    Ta,
                    Tb,
                    q * a.n_channels,
                    area,
                )
            )
        for name, boundary in boundaries.items():
            c = self.circuits[name]
            positions = [midpoint(p.centerline) for p in c.placed_paths]
            for surface in c.cross_section.surface_names:
                inner = midpoint(c.inner_perimeters[surface]) * c.path.lengths
                outer = midpoint(c.outer_perimeters[surface]) * c.path.lengths
                for i, area in enumerate(self._remaining[(name, surface)]):
                    if area <= 1e-15:
                        continue
                    outer_area = outer[i] * area / inner[i]
                    if outer_area <= 0:
                        raise ValueError("exposed outer area must be positive")
                    hc = (
                        results[name].heat_transfer_coefficient[i]
                        * enhancement[name][i]
                    )
                    fluid_T = temperatures[name][i]
                    thickness = midpoint(c.thicknesses[surface])[i]

                    def exchange(wall_T):
                        samples = np.asarray(
                            [boundary.sample(points[i], wall_T) for points in positions]
                        )
                        if not np.allclose(samples, samples[0], rtol=1e-8, atol=1e-8):
                            raise ValueError(
                                f"{name}: copied channels see different gas conditions; use separate circuits"
                            )
                        gas_T, hg = samples[0]
                        return _wall_exchange(
                            fluid_T,
                            gas_T,
                            1 / (hc * area),
                            1 / (hg * outer_area),
                            thickness,
                            0.5 * (area + outer_area),
                            c.wall,
                        )

                    gas_T, _ = boundary.sample(positions[0][i], fluid_T)
                    if np.isclose(gas_T, fluid_T, rtol=0, atol=1e-10):
                        q, Ta, Tb = exchange(fluid_T)
                    else:
                        wall_T = brentq(
                            lambda T: exchange(T)[2] - T,
                            min(fluid_T, gas_T),
                            max(fluid_T, gas_T),
                            xtol=1e-7,
                        )
                        q, Ta, Tb = exchange(wall_T)
                    heat[name][i] += q
                    walls.append(
                        WallResult(
                            name,
                            i,
                            surface,
                            None,
                            None,
                            None,
                            Ta,
                            Tb,
                            q * c.n_channels,
                            area,
                        )
                    )
        return heat, walls
