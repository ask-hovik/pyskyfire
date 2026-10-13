"""Resolve nominal cooling circuits into patternable channel layouts.

The analytical resolver only applies constant angular phase shifts, so every
seed remains an exact copy of its circuit's nominal reference path.  The
numerical resolver retains the visualization packing solver and may deform the
angular coordinate along a seed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property, reduce
from math import gcd
from typing import Any

import numpy as np

from pyskyfire.regen.channel_helpers import (
    SectionPlanes,
    build_axis_normal_section_planes,
    build_section_planes,
    rotate_cylindrical_centerline,
)


__all__ = [
    "AnalyticalLayoutError",
    "ChannelLayout",
    "ResolvedChannel",
    "circular_pattern",
    "collect_channels",
    "collect_keep_outs",
    "gcd_list",
    "redistribute_channels",
    "resolve_analytical_layout",
    "resolve_numerical_layout",
    "tile_centerlines",
    "wedge_info",
]

TWO_PI = 2.0 * np.pi


class AnalyticalLayoutError(ValueError):
    """Raised when nominal circuits cannot form an exact analytical layout."""


@dataclass(frozen=True)
class ResolvedChannel:
    """One unique channel seed and the circular pattern that reproduces it."""

    name: str
    circuit_index: int
    seed_index: int
    reference_path: np.ndarray
    pattern_count: int
    pattern_angle_rad: float
    angular_offset_rad: float = 0.0

    def __post_init__(self) -> None:
        centerline = np.asarray(self.reference_path, dtype=float)
        if centerline.ndim != 2 or centerline.shape[1] != 3 or len(centerline) == 0:
            raise ValueError("Resolved channel centerlines must have shape (N, 3)")
        if not np.all(np.isfinite(centerline)):
            raise ValueError("Resolved channel centerlines must contain finite values")
        if int(self.pattern_count) < 1:
            raise ValueError("A channel pattern must contain at least one instance")
        if not np.isfinite(self.pattern_angle_rad):
            raise ValueError("A channel pattern angle must be finite")
        object.__setattr__(self, "reference_path", centerline)
        object.__setattr__(self, "pattern_count", int(self.pattern_count))

    @cached_property
    def section_planes(self) -> SectionPlanes:
        """Build and cache planes along the reference path."""

        return build_section_planes(self.reference_path)

    @cached_property
    def axis_normal_section_planes(self) -> SectionPlanes:
        """Transverse construction planes for squared channel sections."""
        return build_axis_normal_section_planes(self.reference_path)


@dataclass(frozen=True)
class ChannelLayout:
    """A consumer-neutral collection of unique, circularly patternable seeds."""

    method: str
    channels: tuple[ResolvedChannel, ...]
    metadata: dict[str, Any] = field(default_factory=dict)

    def for_circuits(self, circuit_indices) -> ChannelLayout:
        """Return a layout containing only the requested source circuits."""

        selected = {int(index) for index in circuit_indices}
        available = {channel.circuit_index for channel in self.channels}
        missing = selected - available
        if missing:
            missing_text = ", ".join(str(index) for index in sorted(missing))
            raise IndexError(
                f"Cooling-circuit indices not present in layout: {missing_text}"
            )
        return ChannelLayout(
            method=self.method,
            channels=tuple(
                channel
                for channel in self.channels
                if channel.circuit_index in selected
            ),
            metadata=dict(self.metadata),
        )

    def expanded_channels(self) -> tuple[ResolvedChannel, ...]:
        """Expand every seed into its complete physical circular pattern."""

        expanded = []
        for channel in self.channels:
            for instance_index in range(channel.pattern_count):
                rotation = instance_index * channel.pattern_angle_rad
                expanded.append(
                    ResolvedChannel(
                        name=f"{channel.name}_instance_{instance_index}",
                        circuit_index=channel.circuit_index,
                        seed_index=channel.seed_index,
                        reference_path=rotate_cylindrical_centerline(
                            channel.reference_path,
                            rotation,
                        ),
                        pattern_count=1,
                        pattern_angle_rad=TWO_PI,
                        angular_offset_rad=(channel.angular_offset_rad + rotation),
                    )
                )
        return tuple(expanded)


def gcd_list(integers) -> int:
    """Return the greatest common divisor of an integer sequence."""

    return reduce(gcd, integers) if len(integers) else 1


def wedge_info(thrust_chamber) -> tuple[int, float, list[int]]:
    """Return the numerical layout's symmetry and circuit channel counts.

    Keep-out sector counts join the symmetry, since every wedge has to hold
    the same share of them, but are not listed among the channel counts.
    """

    counts = [
        int(circuit.placement.n_channel_positions)
        for circuit in thrust_chamber.cooling_circuits
    ]
    keep_out_counts = [
        keep_out.count for keep_out in getattr(thrust_chamber, "keep_outs", ())
    ]
    symmetry_count = gcd_list(counts + keep_out_counts)
    wedge_angle = TWO_PI / symmetry_count
    return symmetry_count, wedge_angle, counts


def circular_pattern(centerline: np.ndarray, count: int) -> list[np.ndarray]:
    """Replicate a cylindrical centerline evenly around the engine axis."""

    return [
        rotate_cylindrical_centerline(centerline, angle)
        for angle in np.linspace(0.0, TWO_PI, int(count), endpoint=False)
    ]


def collect_channels(thrust_chamber, *, jitter_divisor: float = 100.0):
    """Collect the unique nominal channels in the numerical symmetry wedge."""

    symmetry_count, wedge_angle, counts = wedge_info(thrust_chamber)

    def jitter_step(angular_width: np.ndarray) -> float:
        angular_width = np.asarray(angular_width, float)
        positive = angular_width[angular_width > 0]
        angular_width_mean = np.nanmean(positive) if positive.size else 0.0
        return (
            angular_width_mean / jitter_divisor if angular_width_mean > 0.0 else 1.0e-6
        )

    centerlines_out = []
    angular_widths_out = []
    circuit_ids_out = []
    global_index = 0

    for circuit_index, circuit in enumerate(thrust_chamber.cooling_circuits):
        base = np.asarray(circuit.reference_path, float)
        count = int(circuit.placement.n_channel_positions)
        angular_width = np.asarray(circuit.channel_angular_width, float)
        step = jitter_step(angular_width)
        unique_count = count // symmetry_count if symmetry_count > 1 else count

        for centerline in circular_pattern(base, count)[:unique_count]:
            centerline = centerline.copy()
            centerline[:, 2] += global_index * step
            centerlines_out.append(centerline)
            angular_widths_out.append(angular_width.copy())
            circuit_ids_out.append(circuit_index)
            global_index += 1

    metadata = {
        "g": symmetry_count,
        "wedge_angle": wedge_angle,
        "n_list": counts,
    }
    return centerlines_out, angular_widths_out, circuit_ids_out, metadata


def collect_keep_outs(thrust_chamber, symmetry_count: int):
    """Collect the keep-out sectors in the numerical symmetry wedge.

    Returns their centerlines and angular widths, sampled at the corners of
    each width law so that linear interpolation reproduces it exactly.
    """

    centerlines_out = []
    angular_widths_out = []
    contour = thrust_chamber.contour
    for keep_out in getattr(thrust_chamber, "keep_outs", ()):
        base = keep_out.centerline(contour)
        widths = np.array(
            [keep_out.sector_angle(x, contour) for x in base[:, 0]], dtype=float
        )
        unique_count = keep_out.count // symmetry_count
        for centerline in circular_pattern(base, keep_out.count)[:unique_count]:
            centerlines_out.append(centerline)
            angular_widths_out.append(widths.copy())
    return centerlines_out, angular_widths_out


def tile_centerlines(
    centerlines_subset,
    angular_widths_subset,
    circuit_ids_subset,
    symmetry_count: int,
    wedge_angle: float,
):
    """Replicate a numerical wedge subset around the circumference."""

    full_centerlines = []
    full_angular_widths = []
    full_circuit_ids = []
    for rotation_index in range(symmetry_count):
        rotation = rotation_index * wedge_angle
        for centerline, angular_width, circuit_id in zip(
            centerlines_subset,
            angular_widths_subset,
            circuit_ids_subset,
            strict=True,
        ):
            full_centerlines.append(rotate_cylindrical_centerline(centerline, rotation))
            full_angular_widths.append(angular_width.copy())
            full_circuit_ids.append(circuit_id)
    return full_centerlines, full_angular_widths, full_circuit_ids


def redistribute_channels(
    centerlines,
    angular_widths,
    *,
    k_pack: float = 1.0,
    k_bend: float = 1.0,
    k_slope: float = 3.0,
    desired_slope=None,
    dt: float = 0.25,
    tol: float = 1.0e-3,
    max_iters: int = 800,
    gamma_damp: float = 0.4,
    wedge_angle: float | None = None,
    recorder: list | None = None,
    endpoint_slope0=None,
    endpoint_slopeN=None,
    fixed=None,
    shape_offsets=None,
):
    """Relax channel angles so they pack without overlap.

    This is the existing visualization solver, now shared by visualization and
    external exporters through :func:`resolve_numerical_layout`. Entries
    flagged in ``fixed`` push on their neighbours but never move; keep-outs
    enter this way. ``shape_offsets`` holds, per entry, an angular detour
    the bending and slope terms treat as straight, so a channel already
    bent around a keep-out is not pulled back into it.
    """

    channel_count = len(centerlines)
    if fixed is None:
        fixed = [False] * channel_count
    if len(fixed) != channel_count:
        raise ValueError("fixed must match the number of centerlines")
    if shape_offsets is None:
        shape_offsets = [None] * channel_count
    if len(shape_offsets) != channel_count:
        raise ValueError("shape_offsets must match the number of centerlines")
    shape_offsets = [
        np.zeros(len(centerline)) if offset is None else np.asarray(offset, float)
        for centerline, offset in zip(centerlines, shape_offsets, strict=True)
    ]
    x_values = [centerline[:, 0] for centerline in centerlines]
    radii = [centerline[:, 1] for centerline in centerlines]
    angles = [centerline[:, 2].copy() for centerline in centerlines]
    widths = [angular_width.copy() for angular_width in angular_widths]
    velocities = [np.zeros_like(angle) for angle in angles]

    if endpoint_slope0 is None:
        endpoint_slope0 = [None] * channel_count
    if endpoint_slopeN is None:
        endpoint_slopeN = [None] * channel_count
    if len(endpoint_slope0) != channel_count or len(endpoint_slopeN) != channel_count:
        raise ValueError(
            "endpoint_slope0/endpoint_slopeN must match the number of centerlines"
        )

    if desired_slope is None:
        desired_slope = [None] * channel_count
    if len(desired_slope) != channel_count:
        raise ValueError("desired_slope must match the number of centerlines")

    def record() -> None:
        if recorder is not None:
            recorder.append(
                [
                    np.column_stack([x_values[i], radii[i], angles[i]]).astype(
                        float,
                        copy=True,
                    )
                    for i in range(channel_count)
                ]
            )

    def apply_hard_endpoint_slopes() -> None:
        for index in range(channel_count):
            if fixed[index]:
                continue
            slope_start = endpoint_slope0[index]
            slope_end = endpoint_slopeN[index]
            if slope_start is None and slope_end is None:
                continue

            angle = angles[index]
            offset = shape_offsets[index]
            axial = x_values[index]
            if angle.size < 2:
                continue
            delta_start = float(axial[1] - axial[0])
            delta_end = float(axial[-1] - axial[-2])
            if abs(delta_start) < 1.0e-14 or abs(delta_end) < 1.0e-14:
                continue
            if slope_start is not None:
                angle[0] = (angle[1] - offset[1] + offset[0]
                            - float(slope_start) * delta_start)
            if slope_end is not None:
                angle[-1] = (angle[-2] - offset[-2] + offset[-1]
                             + float(slope_end) * delta_end)

    master_x = np.unique(np.concatenate(x_values))
    master_count = master_x.size
    active_masks = []
    for axial in x_values:
        minimum, maximum = sorted((axial[0], axial[-1]))
        active_masks.append((master_x >= minimum) & (master_x <= maximum))

    record()
    apply_hard_endpoint_slopes()

    for _ in range(int(max_iters)):
        previous_angles = np.concatenate(angles)
        angle_master = np.zeros((channel_count, master_count), float)
        width_master = np.zeros((channel_count, master_count), float)
        for index in range(channel_count):
            minimum, maximum = sorted((x_values[index][0], x_values[index][-1]))
            interpolated_angle = np.interp(
                master_x,
                x_values[index],
                angles[index],
                left=angles[index][0],
                right=angles[index][-1],
            )
            interpolated_width = np.interp(
                master_x,
                x_values[index],
                widths[index],
                left=widths[index][0],
                right=widths[index][-1],
            )
            active = (master_x >= minimum) & (master_x <= maximum)
            angle_master[index, active] = interpolated_angle[active]
            width_master[index, active] = interpolated_width[active]

        forces_master = np.zeros_like(angle_master)
        for station_index in range(master_count):
            active = [
                index
                for index in range(channel_count)
                if active_masks[index][station_index]
            ]
            if len(active) < 2:
                continue

            if wedge_angle is None:
                station_angles = np.array(
                    [angle_master[index, station_index] for index in active]
                )
                ordered = [active[position] for position in np.argsort(station_angles)]
                for position, first in enumerate(ordered):
                    second = ordered[(position + 1) % len(ordered)]
                    overlap = (
                        angle_master[first, station_index]
                        + 0.5 * width_master[first, station_index]
                        - angle_master[second, station_index]
                        + 0.5 * width_master[second, station_index]
                    )
                    overlap = (overlap + np.pi) % TWO_PI - np.pi
                    forces_master[first, station_index] -= k_pack * overlap
                    forces_master[second, station_index] += k_pack * overlap
            else:
                wedge = float(wedge_angle)
                phases = np.mod(
                    np.array([angle_master[index, station_index] for index in active]),
                    wedge,
                )
                station_widths = np.array(
                    [width_master[index, station_index] for index in active]
                )
                all_phases = np.concatenate([phases - wedge, phases, phases + wedge])
                all_widths = np.concatenate([station_widths] * 3)
                all_ids = np.concatenate([active, active, active])
                order = np.argsort(all_phases)
                middle_positions = [
                    position
                    for position in order
                    if 0.0 <= all_phases[position] < wedge
                ]

                def wrap_wedge(value):
                    return (value + 0.5 * wedge) % wedge - 0.5 * wedge

                for position in middle_positions:
                    first = int(all_ids[position])
                    ordered_index = np.where(order == position)[0][0]
                    next_position = order[(ordered_index + 1) % len(order)]
                    second = int(all_ids[next_position])
                    overlap = wrap_wedge(
                        all_phases[position]
                        + 0.5 * all_widths[position]
                        - all_phases[next_position]
                        + 0.5 * all_widths[next_position]
                    )
                    forces_master[first, station_index] -= k_pack * overlap
                    forces_master[second, station_index] += k_pack * overlap

        for index in range(channel_count):
            if fixed[index]:
                continue
            force = np.interp(
                x_values[index],
                master_x,
                forces_master[index],
                left=0.0,
                right=0.0,
            )
            angle = angles[index]
            shape = angle - shape_offsets[index]
            laplacian = np.zeros_like(angle)
            if angle.size > 2:
                laplacian[1:-1] = shape[:-2] - 2.0 * shape[1:-1] + shape[2:]
            if angle.size > 1:
                laplacian[0] = 0.0
                laplacian[-1] = 0.0
            force += k_bend * laplacian

            if k_slope > 0.0 and desired_slope[index] is not None and angle.size > 2:
                desired = np.asarray(desired_slope[index], dtype=float)
                if desired.shape[0] != angle.shape[0]:
                    raise ValueError(
                        f"desired_slope[{index}] must match its centerline length"
                    )
                segment_dx = np.diff(x_values[index])
                slope = np.diff(shape) / (segment_dx + 1.0e-16)
                desired_segment = 0.5 * (desired[:-1] + desired[1:])
                error = slope - desired_segment
                gradient = np.zeros_like(angle)
                weights = segment_dx / (np.mean(segment_dx) + 1.0e-16)
                gradient[1:-1] = weights[:-1] * error[:-1] - weights[1:] * error[1:]
                gradient[0] = 0.0
                gradient[-1] = 0.0
                force += -k_slope * gradient

            velocity = (1.0 - gamma_damp) * velocities[index] + dt * force
            angles[index] = angle + dt * velocity
            velocities[index] = velocity

        apply_hard_endpoint_slopes()
        record()
        if np.max(np.abs(np.concatenate(angles) - previous_angles)) < tol:
            break

    resolved = []
    for index, centerline in enumerate(centerlines):
        result = centerline.copy()
        result[:, 2] = angles[index]
        resolved.append(result)
    return resolved


def _seat_around_keep_outs(
    centerlines,
    circuit_ids,
    obstacles,
    obstacle_widths,
    counts,
    period: float,
):
    """Open room at each keep-out before the packer relaxes the rest.

    Relaxation alone cannot do this: a sector several pitches wide appears on
    top of a fully packed row, and the packer would push channels through one
    another. So each circuit crossing a keep-out is first turned by a constant
    phase, less than one pitch, that puts the keep-out centre between two of
    its channels rather than on one. Then, at every station, each stretch of
    channels between two keep-out centres is scaled linearly onto the free
    arc between their edges. Where a keep-out has no width the map is the
    identity, so channels leave the nominal pattern smoothly along a taper.
    """

    if not obstacles:
        return centerlines
    spans = [(obstacle[0, 0], obstacle[-1, 0]) for obstacle in obstacles]
    seated = [centerline.copy() for centerline in centerlines]

    turned = {}
    for circuit_index in sorted(set(circuit_ids)):
        seeds = [
            index for index, owner in enumerate(circuit_ids) if owner == circuit_index
        ]
        axial = seated[seeds[0]][:, 0]
        low, high = min(axial[0], axial[-1]), max(axial[0], axial[-1])
        crossing = [
            index
            for index, (start, end) in enumerate(spans)
            if start < high and end > low
        ]
        if not crossing:
            continue
        obstacle = obstacles[crossing[0]]
        middle = float(np.clip(0.5 * sum(spans[crossing[0]]), low, high))
        order = np.argsort(axial)
        seed_angle = float(np.interp(middle, axial[order], seated[seeds[0]][order, 2]))
        pitch = TWO_PI / counts[circuit_index]
        turned[circuit_index] = (obstacle[0, 2] - 0.5 * pitch - seed_angle) % pitch
        # The seeds are laid on the exact pitch from the first one. The
        # tie-breaking jitter grows with each seed's index, and would carry
        # the last of them into the keep-out; the first keeps its share, which
        # still separates this circuit from any other.
        first = seated[seeds[0]][:, 2].copy()
        for position, index in enumerate(seeds):
            seated[index][:, 2] = first + position * pitch + turned[circuit_index]

    for centerline in seated:
        for station, x in enumerate(centerline[:, 0]):
            active = [
                (obstacle[0, 2] % period, 0.5 * float(np.interp(x, obstacle[:, 0], width)))
                for obstacle, width, (start, end) in zip(
                    obstacles, obstacle_widths, spans, strict=True
                )
                if start <= x <= end
            ]
            if not active:
                continue
            active.sort()
            centres = np.array([centre for centre, _ in active])
            halves = np.array([half for _, half in active])
            angle = centerline[station, 2]
            phase = angle % period
            after = int(np.searchsorted(centres, phase, side="right")) - 1
            start_centre = centres[after]
            start_half = halves[after]
            end_centre = centres[(after + 1) % len(centres)]
            end_half = halves[(after + 1) % len(centres)]
            nominal = (end_centre - start_centre) % period or period
            free = nominal - start_half - end_half
            offset = (phase - start_centre) % period
            mapped = start_half + offset * free / nominal
            centerline[station, 2] = angle + (mapped - offset)
    return seated


def _placement_slope_data(thrust_chamber, centerlines, circuit_ids):
    desired_slopes = []
    slope_starts = []
    slope_ends = []
    contour = thrust_chamber.contour
    for centerline, circuit_index in zip(centerlines, circuit_ids, strict=True):
        circuit = thrust_chamber.cooling_circuits[circuit_index]
        placement = circuit.placement
        axial = centerline[:, 0]
        slope = np.asarray(
            [placement.dtheta_dx(float(x), contour, circuit) for x in axial],
            dtype=float,
        )
        desired_slopes.append(slope)
        slope_starts.append(float(slope[0]))
        slope_ends.append(float(slope[-1]))
    return desired_slopes, slope_starts, slope_ends


def resolve_numerical_layout(thrust_chamber, **solver_options) -> ChannelLayout:
    """Resolve the historical numerical packing layout.

    The result stores only the unique seeds in the greatest-common-divisor
    wedge.  Consumers may call :meth:`ChannelLayout.expanded_channels` when
    they need all physical instances.
    """

    centerlines, widths, circuit_ids, metadata = collect_channels(thrust_chamber)
    if not centerlines:
        return ChannelLayout("numerical", (), metadata)
    desired_slopes, slope_starts, slope_ends = _placement_slope_data(
        thrust_chamber,
        centerlines,
        circuit_ids,
    )
    # Keep-outs ride along as pinned obstacles and are dropped afterwards.
    obstacles, obstacle_widths = collect_keep_outs(thrust_chamber, metadata["g"])
    channel_count = len(centerlines)
    padding = [None] * len(obstacles)
    metadata["keep_out_sectors"] = len(obstacles)
    period = metadata["wedge_angle"] if metadata["g"] > 1 else TWO_PI
    seated = _seat_around_keep_outs(
        centerlines,
        circuit_ids,
        obstacles,
        obstacle_widths,
        metadata["n_list"],
        period,
    )
    detours = [
        after[:, 2] - before[:, 2]
        for before, after in zip(centerlines, seated, strict=True)
    ]
    centerlines = seated
    options = {
        "k_pack": 2.0,
        "k_bend": 0.1,
        "k_slope": 0.05,
        "dt": 0.15,
        "tol": 1.0e-5,
        "max_iters": 800,
        "gamma_damp": 0.4,
        "wedge_angle": metadata["wedge_angle"] if metadata["g"] > 1 else None,
        "endpoint_slope0": slope_starts + padding,
        "endpoint_slopeN": slope_ends + padding,
        "desired_slope": desired_slopes + padding,
        "fixed": [False] * channel_count + [True] * len(obstacles),
        "shape_offsets": detours + padding,
    }
    options.update(solver_options)
    resolved = redistribute_channels(
        centerlines + obstacles,
        widths + obstacle_widths,
        **options,
    )[:channel_count]

    seed_numbers: dict[int, int] = {}
    channels = []
    for centerline, circuit_index in zip(resolved, circuit_ids, strict=True):
        seed_index = seed_numbers.get(circuit_index, 0)
        seed_numbers[circuit_index] = seed_index + 1
        circuit = thrust_chamber.cooling_circuits[circuit_index]
        channels.append(
            ResolvedChannel(
                name=f"{circuit.name}_seed_{seed_index}",
                circuit_index=circuit_index,
                seed_index=seed_index,
                reference_path=centerline,
                pattern_count=metadata["g"],
                pattern_angle_rad=metadata["wedge_angle"],
            )
        )
    return ChannelLayout("numerical", tuple(channels), metadata)


def _axial_span(circuit) -> tuple[float, float]:
    centerline = np.asarray(circuit.reference_path, dtype=float)
    return tuple(sorted((float(centerline[0, 0]), float(centerline[-1, 0]))))


def _positive_overlap(first_span, second_span, tolerance: float) -> bool:
    return (
        min(first_span[1], second_span[1]) - max(first_span[0], second_span[0])
        > tolerance
    )


def _overlap_components(circuits, *, tolerance: float) -> list[list[int]]:
    occluding = [
        index
        for index, circuit in enumerate(circuits)
        if bool(getattr(circuit.placement, "occludes", True))
    ]
    spans = {index: _axial_span(circuits[index]) for index in occluding}
    adjacency = {index: set() for index in occluding}
    for position, first in enumerate(occluding):
        for second in occluding[position + 1 :]:
            if _positive_overlap(spans[first], spans[second], tolerance):
                adjacency[first].add(second)
                adjacency[second].add(first)

    components = []
    unseen = set(occluding)
    while unseen:
        start = min(unseen)
        stack = [start]
        component = []
        unseen.remove(start)
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbor in sorted(adjacency[current], reverse=True):
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    stack.append(neighbor)
        components.append(sorted(component))

    for index in range(len(circuits)):
        if index not in occluding:
            components.append([index])
    return sorted(components, key=lambda component: component[0])


def _validate_matching_angular_evolution(
    first,
    second,
    *,
    span_tolerance: float,
    angular_rate_tolerance: float,
) -> float:
    """Validate the two nominal winding laws and return their phase gap."""
    first_centerline = np.asarray(first.reference_path, dtype=float)
    second_centerline = np.asarray(second.reference_path, dtype=float)
    lower = max(_axial_span(first)[0], _axial_span(second)[0])
    upper = min(_axial_span(first)[1], _axial_span(second)[1])
    if upper - lower <= span_tolerance:
        return 0.0

    samples = np.unique(
        np.concatenate(
            (
                first_centerline[
                    (first_centerline[:, 0] >= lower)
                    & (first_centerline[:, 0] <= upper),
                    0,
                ],
                second_centerline[
                    (second_centerline[:, 0] >= lower)
                    & (second_centerline[:, 0] <= upper),
                    0,
                ],
                np.array([lower, upper]),
            )
        )
    )
    first_rate = np.asarray(
        [
            first.placement.dtheta_dx(float(x), getattr(first, "contour", None), first)
            for x in samples
        ]
    )
    second_rate = np.asarray(
        [
            second.placement.dtheta_dx(float(x), getattr(second, "contour", None), second)
            for x in samples
        ]
    )
    if not np.allclose(
        first_rate,
        second_rate,
        atol=angular_rate_tolerance,
        rtol=0.0,
    ):
        raise AnalyticalLayoutError(
            f"Circuits {first.name!r} and {second.name!r} do not follow the "
            "same angular/helix evolution over their overlap"
        )

    first_origin = getattr(first, "reference_angle_origin_rad", None)
    second_origin = getattr(second, "reference_angle_origin_rad", None)
    if first_origin is not None and second_origin is not None:
        return float(first_origin) - float(second_origin)

    def reference_angle_at(centerline, x):
        order = np.argsort(centerline[:, 0])
        axial = centerline[order, 0]
        angles = np.unwrap(centerline[order, 2])
        return float(np.interp(x, axial, angles))

    first_angle = reference_angle_at(first_centerline, lower)
    second_angle = reference_angle_at(second_centerline, lower)
    return second_angle - first_angle


def resolve_analytical_layout(
    thrust_chamber,
    *,
    span_tolerance: float = 1.0e-9,
    angular_rate_tolerance: float = 1.0e-7,
) -> ChannelLayout:
    """Resolve an exact nominal layout using only constant phase rotations.

    Overlapping, occluding circuits are inferred to interleave.  Such circuits
    must have equal channel counts and matching reference-path angular
    rates.  The coolant-centroid path is intentionally ignored.  Two may
    have different axial endpoints; groups of three or more must have identical
    spans.  The first circuit has zero added phase and later circuits divide a
    single channel pitch evenly.
    """

    circuits = thrust_chamber.cooling_circuits
    if not circuits:
        return ChannelLayout("analytical", (), {"groups": []})
    if getattr(thrust_chamber, "keep_outs", ()):
        raise AnalyticalLayoutError(
            "Channels can only bend around keep-outs in the numerical layout"
        )
    for circuit in circuits:
        if circuit.reference_path is None or len(circuit.reference_path) == 0:
            raise AnalyticalLayoutError(
                f"Circuit {circuit.name!r} has no reference path"
            )

    groups = _overlap_components(circuits, tolerance=span_tolerance)
    phases: dict[int, float] = {}
    group_metadata = []
    for group in groups:
        counts = [int(circuits[index].placement.n_channel_positions) for index in group]
        if len(set(counts)) != 1:
            names = ", ".join(repr(circuits[index].name) for index in group)
            raise AnalyticalLayoutError(
                f"Interleaving circuits must have equal channel counts: {names}"
            )

        if len(group) >= 3:
            reference_span = _axial_span(circuits[group[0]])
            for index in group[1:]:
                if not np.allclose(
                    _axial_span(circuits[index]),
                    reference_span,
                    atol=span_tolerance,
                    rtol=0.0,
                ):
                    names = ", ".join(repr(circuits[item].name) for item in group)
                    raise AnalyticalLayoutError(
                        "Three or more interleaving circuits must have identical "
                        f"axial spans: {names}"
                    )

        nominal_phase_from_first = {group[0]: 0.0}
        for position, first_index in enumerate(group):
            for second_index in group[position + 1 :]:
                phase_difference = _validate_matching_angular_evolution(
                    circuits[first_index],
                    circuits[second_index],
                    span_tolerance=span_tolerance,
                    angular_rate_tolerance=angular_rate_tolerance,
                )
                if first_index == group[0]:
                    nominal_phase_from_first[second_index] = phase_difference

        count = counts[0]
        combined_pitch = TWO_PI / (count * len(group))
        for phase_index, circuit_index in enumerate(group):
            phases[circuit_index] = (
                phase_index * combined_pitch
                - nominal_phase_from_first.get(circuit_index, 0.0)
            )
        group_metadata.append(
            {
                "circuit_indices": list(group),
                "channel_count": count,
                "phase_step_rad": combined_pitch,
            }
        )

    channels = []
    for circuit_index, circuit in enumerate(circuits):
        count = int(circuit.placement.n_channel_positions)
        phase = phases[circuit_index]
        centerline = rotate_cylindrical_centerline(circuit.reference_path, phase)
        channels.append(
            ResolvedChannel(
                name=f"{circuit.name}_seed_0",
                circuit_index=circuit_index,
                seed_index=0,
                reference_path=centerline,
                pattern_count=count,
                pattern_angle_rad=TWO_PI / count,
                angular_offset_rad=phase,
            )
        )

    return ChannelLayout(
        method="analytical",
        channels=tuple(channels),
        metadata={"groups": group_metadata},
    )
