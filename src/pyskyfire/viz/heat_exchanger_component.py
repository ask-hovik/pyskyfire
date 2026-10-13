"""Heat-exchanger geometry and results in the shared scene hierarchy."""

from dataclasses import dataclass
import numpy as np
from pyskyfire.heat_exchanger import (
    ChamberInstallation,
    HeatExchanger,
    HeatExchangerResult,
)
from pyskyfire.thermo.eos import fluid_name
from .components import _carrying, _surface
from .scene import DEFAULT_COLOR, Anchor, FieldData, SceneNode, Transform

__all__ = ["HeatExchangerComponent"]


def _fields(circuit, result, walls, samples):
    """Map station data and area-averaged wetted wall patches onto the mesh."""
    n = len(circuit.path.stations)
    fields = {}
    for name, values, units, label in (
        ("T_static", result.T, "K", "Static coolant temperature"),
        ("p_static", result.p, "Pa", "Static pressure"),
        ("velocity", result.velocity, "m/s", "Coolant velocity"),
        ("coolant_enthalpy", result.h, "J/kg", "Coolant enthalpy"),
    ):
        values = np.asarray(values, float)
        if values.shape != (n,):
            raise ValueError(f"{name} must have one value per heat-exchanger station")
        fields[name] = FieldData(np.repeat(values, 4 * samples), units, label)
    weight = np.zeros((n - 1, 4))
    inner = np.zeros_like(weight)
    names = circuit.cross_section.surface_names
    for wall in walls:
        if wall.circuit_a == circuit.name:
            i, surface, Ti = wall.segment_a, wall.surface_a, wall.temperature_a
        elif wall.circuit_b == circuit.name:
            i, surface, Ti = wall.segment_b, wall.surface_b, wall.temperature_b
        else:
            continue
        if (
            not isinstance(i, (int, np.integer))
            or not 0 <= i < n - 1
            or surface not in names
        ):
            raise ValueError(
                "wall result does not identify a heat-exchanger surface segment"
            )
        if not np.isfinite(wall.area) or wall.area <= 0:
            raise ValueError("wall result contact area must be positive and finite")
        j = names.index(surface)
        weight[i, j] += wall.area
        inner[i, j] += wall.area * Ti
    exposure = circuit.cross_section.exposure_fractions(circuit.profiles)
    edge_midpoints = (np.arange(samples) + 0.5) / samples
    values = np.divide(inner, weight, out=np.full_like(weight, np.nan), where=weight > 0)
    per_edge = np.repeat(values, samples, axis=1)
    for j, surface in enumerate(names):
        fraction = 0.5 * (exposure[surface][:-1] + exposure[surface][1:])
        active = np.abs(edge_midpoints[None, :] - 0.5) <= fraction[:, None] / 2
        per_edge[:, j * samples : (j + 1) * samples] = np.where(
            active, per_edge[:, j * samples : (j + 1) * samples], np.nan
        )
    fields["T_coolant_wall"] = FieldData(
        np.repeat(per_edge.ravel(), 2), "K", "Coolant wall temperature", "cell"
    )
    return fields


@dataclass
class HeatExchangerComponent:
    """Display an exchanger or chamber installation in the shared Scene.

    Each channel is drawn as its wetted passage, the surface the coolant
    touches, not the outer envelope of its walls. Raw exchangers use assembly
    coordinates. Installations carry their simulation transform, relative to
    the parent chamber; explicit Scene.add placement
    overrides it. Each channel has inlet/outlet anchors. Root anchors
    circuit_N.inlet/outlet refer to the first copy; origin is the assembly datum.

    Optional results provide T_static, p_static, velocity, coolant_enthalpy
    and T_coolant_wall, the wall temperature at the wetted interface. Wall
    temperatures are area averages per exposed surface and segment. Inactive
    portions remain NaN. The opposite (gas-side) interface is not drawn.
    With results, each circuit records its CoolProp fluid name as
    ``metadata["fluid"]``. ``circuit_colors`` paints each circuit apart;
    without it every channel takes the default colour.
    """

    heat_exchanger: object
    results: HeatExchangerResult | None = None
    samples: int = 8
    circuit_colors: object = None

    def build(self):
        if (
            isinstance(self.samples, bool)
            or not isinstance(self.samples, int)
            or self.samples < 2
        ):
            raise ValueError("samples must be an integer of at least two")
        installation = (
            self.heat_exchanger
            if isinstance(self.heat_exchanger, ChamberInstallation)
            else None
        )
        exchanger = (
            installation.exchanger if installation is not None else self.heat_exchanger
        )
        if not isinstance(exchanger, HeatExchanger):
            raise TypeError(
                "HeatExchangerComponent expects HeatExchanger or ChamberInstallation"
            )
        if self.results is not None:
            if not isinstance(self.results, HeatExchangerResult):
                raise TypeError("results must be a HeatExchangerResult")
            if set(self.results.circuits) != set(exchanger.circuits):
                raise ValueError("results must match the heat-exchanger circuits")
            for name, result in self.results.circuits.items():
                if result.circuit is not exchanger.circuits[name]:
                    raise ValueError(
                        f"results for {name!r} belong to a different circuit"
                    )
        if self.circuit_colors is not None and len(self.circuit_colors) < len(exchanger.circuits):
            raise ValueError("circuit_colors must provide a color for every circuit")
        root = SceneNode(
            "heat_exchanger",
            anchors={"origin": Anchor()},
            metadata={"kind": "heat_exchanger"},
        )
        if installation is not None:
            matrix = np.eye(4)
            matrix[:3, :3] = installation.gas.rotation
            matrix[:3, 3] = installation.gas.position
            root.transform = Transform(matrix=matrix)
        for ci, (name, circuit) in enumerate(exchanger.circuits.items()):
            group = SceneNode(
                f"circuit_{ci}",
                metadata=_carrying(
                    {
                        "kind": "heat_exchanger_circuit",
                        "circuit_index": ci,
                        "label": name,
                    },
                    None if self.results is None else fluid_name(getattr(circuit, "transport", None)),
                ),
            )
            root.children.append(group)
            fields = (
                {}
                if self.results is None
                else _fields(
                    circuit,
                    self.results.circuits[name],
                    self.results.walls,
                    self.samples,
                )
            )
            meshes = {}
            for i, placement in enumerate(circuit.placements):
                # A different roll on a curved template cannot in general be
                # reproduced by rigidly rotating the entire centerline mesh.
                local_up = placement.rotation.T @ circuit.frames[i][0, :, 0]
                key = tuple(np.round(local_up, 12))
                matrix = np.eye(4)
                matrix[:3, :3] = placement.rotation
                matrix[:3, 3] = (
                    np.asarray(placement.start_position)
                    - placement.rotation @ circuit.path.centerline[0]
                )
                transform = Transform(matrix=matrix)
                if key not in meshes:
                    sections = transform.inverse().apply(
                        circuit.section_points(i, exterior=False, samples=self.samples)
                    )
                    meshes[key] = _surface(sections, closed=True)
                anchors = {
                    "inlet": Anchor.from_direction(
                        circuit.path.centerline[0], -circuit.path.tangents[0]
                    ),
                    "outlet": Anchor.from_direction(
                        circuit.path.centerline[-1], circuit.path.tangents[-1]
                    ),
                }
                group.children.append(
                    SceneNode(
                        f"channel_{i}",
                        mesh=meshes[key],
                        fields=dict(fields),
                        anchors=anchors,
                        transform=transform,
                        color=(
                            DEFAULT_COLOR
                            if self.circuit_colors is None
                            else self.circuit_colors[ci]
                        ),
                        metadata={
                            "kind": "heat_exchanger_channel",
                            "circuit_index": ci,
                            "instance_index": i,
                        },
                    )
                )
                if i == 0:
                    group.anchors = {
                        k: Anchor(transform @ a.transform) for k, a in anchors.items()
                    }
                    root.anchors.update(
                        {f"circuit_{ci}.{k}": a for k, a in group.anchors.items()}
                    )
        return root
