"""Build and solve the RL10 turbine-discharge duct from Figure G1.

The duct runs from the turbine housing discharge to the main fuel shut-off
valve. Figure G1 is the tube-bending table for the part and Table G1 is the
pressure-loss study built on it; between them, and the fabrication drawing
views, the tube is very nearly fully specified.

Points 1 and 6 are the tube ends. Points 2 to 5 are the bend locations, each
one the intersection vertex of the two straight control legs that meet there.

Bore diameter
-------------
The drawing gives tube OD 2.510/2.490 and wall 0.054/0.044 (0.054/0.035 at
bends), so the nominal bore is 2.500 - 2 x 0.049 = 2.402 in. Section A-A
states 2.408/2.396 AVG DIA independently, and the two agree exactly. That
average bore is the flow diameter used here.

Bend planes
-----------
Figure G1's tube-rotation column gives the dihedral angle between consecutive
bend planes: 71 deg 06', 65 deg 39', 116 deg 35'. The magnitudes are
tabulated; the sign of each rotation is not.

The drawing narrows it. Its end view draws the two tube ends at a clear angle
to each other -- one reads as a circle, the other as an ellipse -- so any
layout leaving the end legs near-collinear is wrong. That eliminates two of
the four sign patterns outright:

    (+, +, +)   end legs  20.4 deg apart   -- rejected, near-collinear
    (+, -, -)   end legs 159.6 deg apart   -- rejected, near-collinear
    (+, +, -)   end legs  53.8 deg apart   -- possible
    (+, -, +)   end legs 126.5 deg apart   -- possible

The two survivors are near-supplementary, so they foreshorten the far end by
the same factor and produce an equally convincing ellipse. The end view cannot
tell them apart that way.

What does separate them, weakly, is the drawn silhouette. Measured off the
scan, the end view inks roughly 25.2 x 16.9 in. Searching in-plane rotations
of each candidate's projection along leg 1-2 for that box gives a total
mismatch of 1.5 in for ``(+, +, -)`` against 3.65 in for ``(+, -, +)``, so
``(+, +, -)`` is the default here. That measurement is worth maybe one or two
inches either way, so treat this as a leaning rather than a determination:
the report's "Bend planes" tab draws both candidates in the same views for
comparison against the drawing, and ``ROLL_SIGNS`` switches between them.

This ambiguity is cosmetic. Roll angles set the three-dimensional shape, but
the loss model sees only lengths, angles, radii, and diameter, so the pressure
drop is identical for every sign pattern.

Assumed
-------
Only the operating point: representative RL10A-3-3A turbine-discharge hydrogen
conditions. The geometry and the 4.6e-5 in roughness height are from the
source.

Two connection models run over the identical centerline:

1. :class:`~pyskyfire.ducts.RoundedBend` -- one lumped fitting loss per bend
   from ``fluids.fittings.bend_rounded``. This matches what Table G1 is doing,
   and is the model to use for discrete piping.
2. :class:`~pyskyfire.ducts.Continuous` -- the cooling-channel helical-coil
   curvature multiplier applied to a discretized arc. Shown for contrast; it
   models curvature as enhanced wall friction and is *not* a discrete-bend
   model, so it runs outside its documented range here and says so.
"""

import argparse
import math
from pathlib import Path
import warnings

import numpy as np
import plotly.graph_objects as go

import pyskyfire as psf


INCH = 0.0254


def dms(degrees: float, minutes: float) -> float:
    """Degrees and arcminutes, as the drawing states angles, in decimal degrees."""
    return degrees + minutes / 60.0


#: Figure G1 "true length of straight" [in], tangent point to tangent point,
#: ordered 1-2 through 5-6.
STRAIGHTS = (0.185, 2.825, 1.519, 3.572, 6.452)

#: Figure G1 bends: centerline radius [in] and turn angle [deg], at points 2-5.
BENDS = (
    (3.750, dms(97, 31)),
    (12.000, dms(67, 30)),
    (6.250, dms(104, 26)),
    (10.000, dms(29, 54)),
)

#: Figure G1 tube rotations [deg]: the dihedral angle from one bend plane to
#: the next. Zero would leave two bends coplanar and turning the same way.
ROLL_MAGNITUDES = (dms(71, 6), dms(65, 39), dms(116, 35))

#: Sign patterns not excluded by the drawing, keyed by the resulting angle
#: between the two end legs. See the module docstring.
CANDIDATE_ROLL_SIGNS = {"obtuse ends": (1, -1, 1), "acute ends": (1, 1, -1)}
ROLL_SIGNS = CANDIDATE_ROLL_SIGNS["acute ends"]

#: Silhouette of the drawing's end view, long then short side [in], measured
#: off the scan. Good to an inch or two, which is why it only leans.
DRAWING_END_VIEW_EXTENT = (25.2, 16.9)

#: Figure G1's tangent-distance column [in]. Independently equal to
#: ``R tan(angle / 2)``, which is what makes it a check on the columns left.
TABULATED_TANGENT_DISTANCES = (4.277, 8.018, 8.062, 2.670)

#: Table G1 equivalent lengths [in]: manual estimate, then software estimate.
REFERENCE_EQUIVALENT_LENGTHS = ((35.21, 34.8), (33.36, 28.80), (30.63, 32.40), (5.212, 13.20))
REFERENCE_TOTAL_EQUIVALENT_LENGTH = (119.0, 123.8)

#: Table G1's stated actual duct length [in].
REFERENCE_DUCT_LENGTH = 51.7

#: Drawing bore: OD 2.510/2.490 less twice the 0.054/0.044 wall, and equal to
#: the 2.408/2.396 AVG DIA of Section A-A [in].
BORE_DIAMETER = 2.402
OUTER_DIAMETER = 2.500
WALL_THICKNESS = 0.049

#: Roughness height stated with Table G1 [m].
ROUGHNESS = 4.6e-5 * INCH

#: Representative RL10A-3-3A turbine-discharge hydrogen conditions.
INLET_TEMPERATURE = 130.0
INLET_PRESSURE = 4.1e6
MASS_FLOW = 2.77

TUBE_COLOR = "#333333"
AXIS_COLOR = "#888888"
COLORS = ("#0072B2", "#D55E00")


# --- geometry ----------------------------------------------------------------


def rotate(vector, axis, angle):
    """Rotate ``vector`` about ``axis`` by ``angle`` [rad] (Rodrigues)."""
    axis = axis / np.linalg.norm(axis)
    return (
        vector * math.cos(angle)
        + np.cross(axis, vector) * math.sin(angle)
        + axis * np.dot(axis, vector) * (1.0 - math.cos(angle))
    )


def rounded_bend(radius):
    """One lumped ``fluids`` rounded-bend loss between the tangent stations."""
    return psf.ducts.RoundedBend(radius=radius)


def build_route(diameter, bend, signs=ROLL_SIGNS):
    """Turn the Figure G1 bending table into ``PathNode`` control vertices.

    The table gives tangent-to-tangent straight lengths, but a ``PathNode``
    sits at the *intersection* of two straight control legs. Each corner is
    therefore offset from its tangent points by the tangent distance
    :math:`R \\tan(\\theta/2)`, which is added to both neighbouring legs so the
    resolved straights come back out at their tabulated lengths.

    The bend plane is carried by the bend's rotation axis, the binormal. Each
    tabulated tube rotation rolls that binormal about the intervening straight.

    ``bend`` is called with a radius [m] and returns the connection object
    attached to that corner.
    """
    position = np.zeros(3)
    tangent = np.array([1.0, 0.0, 0.0])
    binormal = np.array([0.0, 0.0, 1.0])

    nodes = [psf.ducts.PathNode(tuple(position), diameter)]
    leg = STRAIGHTS[0] * INCH

    for index, (radius_inches, angle_degrees) in enumerate(BENDS):
        if index > 0:
            roll = math.radians(signs[index - 1] * ROLL_MAGNITUDES[index - 1])
            binormal = rotate(binormal, tangent, roll)
            binormal -= np.dot(binormal, tangent) * tangent
            binormal /= np.linalg.norm(binormal)

        angle = math.radians(angle_degrees)
        radius = radius_inches * INCH
        tangent_distance = radius * math.tan(0.5 * angle)

        position = position + (leg + tangent_distance) * tangent
        nodes.append(psf.ducts.PathNode(tuple(position), diameter, bend(radius)))

        tangent = rotate(tangent, binormal, angle)
        leg = tangent_distance + STRAIGHTS[index + 1] * INCH

    nodes.append(psf.ducts.PathNode(tuple(position + leg * tangent), diameter))
    return nodes


def build_pipe(diameter, bend, signs=ROLL_SIGNS):
    """Resolve the route into a ``Pipe`` carrying the tabulated roughness."""
    stations, connections = psf.ducts.generate_pipe_stations(build_route(diameter, bend, signs))
    return psf.ducts.Pipe(stations, connections, roughness=ROUGHNESS)


def route_vertices(diameter, signs=ROLL_SIGNS):
    """Drawing points 1-6 in inches."""
    nodes = build_route(diameter, rounded_bend, signs)
    return np.asarray([node.position for node in nodes]) / INCH


def dense_centerline(diameter, signs=ROLL_SIGNS, segments=60):
    """A finely resolved centerline in inches, for drawing rather than solving."""
    pipe = build_pipe(diameter, lambda radius: psf.ducts.Continuous(radius, segments), signs)
    return pipe.centerline / INCH


def end_leg_angle(diameter, signs=ROLL_SIGNS) -> float:
    """Angle between the two end legs [deg]; the drawing shows it is not small."""
    vertices = route_vertices(diameter, signs)
    inlet = vertices[1] - vertices[0]
    outlet = vertices[5] - vertices[4]
    cosine = float(np.dot(inlet, outlet) / (np.linalg.norm(inlet) * np.linalg.norm(outlet)))
    return math.degrees(math.acos(np.clip(cosine, -1.0, 1.0)))


def part_frame(vertices):
    """A reproducible part-intrinsic frame, close to how the drawing sits.

    ``x`` runs along the long end leg 5-6, ``z`` is the perpendicular component
    of the tall side of the hook, and ``y`` completes the right-handed set.
    """
    x = vertices[5] - vertices[4]
    x /= np.linalg.norm(x)
    z = vertices[1] - vertices[4]
    z -= np.dot(z, x) * x
    z /= np.linalg.norm(z)
    return x, np.cross(z, x), z


# --- orthographic drawing ----------------------------------------------------


def tube_outline(centerline, radius, horizontal, vertical):
    """Silhouette of a circular tube in orthographic projection.

    A circular tube of radius ``r`` has a silhouette exactly ``2r`` wide
    measured across its projected centerline, whatever the viewing angle, so
    each edge is the projected centerline offset by ``r`` along the projected
    normal.
    """
    projected = np.stack([centerline @ horizontal, centerline @ vertical], axis=1)
    tangent = np.gradient(projected, axis=0)
    length = np.linalg.norm(tangent, axis=1, keepdims=True)
    tangent = tangent / np.where(length > 0.0, length, 1.0)
    normal = np.stack([-tangent[:, 1], tangent[:, 0]], axis=1)
    return projected + radius * normal, projected - radius * normal


def end_face(position, tangent, radius, horizontal, vertical, points=72):
    """The circular end face, which projects to a circle or an ellipse."""
    reference = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(reference, tangent)) > 0.9:
        reference = np.array([1.0, 0.0, 0.0])
    first = np.cross(tangent, reference)
    first /= np.linalg.norm(first)
    second = np.cross(tangent, first)
    theta = np.linspace(0.0, 2.0 * math.pi, points)
    circle = position + radius * (np.outer(np.cos(theta), first) + np.outer(np.sin(theta), second))
    return circle @ horizontal, circle @ vertical


def make_orthographic_figure(centerline, vertices, radius, frame, view, title):
    """One orthographic view, drawn as tube walls, centerline, and points 1-6."""
    horizontal, vertical, labels = view
    figure = go.Figure()

    for edge in tube_outline(centerline, radius, horizontal, vertical):
        figure.add_trace(
            go.Scatter(
                x=edge[:, 0],
                y=edge[:, 1],
                mode="lines",
                line=dict(color=TUBE_COLOR, width=2),
                showlegend=False,
                hoverinfo="skip",
            )
        )

    for index in (0, -1):
        tangent = vertices[1] - vertices[0] if index == 0 else vertices[5] - vertices[4]
        tangent = tangent / np.linalg.norm(tangent)
        x, y = end_face(vertices[index], tangent, radius, horizontal, vertical)
        figure.add_trace(
            go.Scatter(
                x=x,
                y=y,
                mode="lines",
                line=dict(color=TUBE_COLOR, width=2),
                showlegend=False,
                hoverinfo="skip",
            )
        )

    figure.add_trace(
        go.Scatter(
            x=centerline @ horizontal,
            y=centerline @ vertical,
            mode="lines",
            line=dict(color=AXIS_COLOR, width=1.5, dash="dashdot"),
            name="Centerline",
            hoverinfo="skip",
        )
    )
    figure.add_trace(
        go.Scatter(
            x=vertices @ horizontal,
            y=vertices @ vertical,
            mode="markers+text",
            marker=dict(color=COLORS[1], size=8, symbol="circle-open", line=dict(width=2)),
            text=[str(index + 1) for index in range(len(vertices))],
            textposition="top right",
            textfont=dict(color=COLORS[1], size=13),
            name="Points 1-6",
            hovertemplate="point %{text}<br>%{x:.3f}, %{y:.3f} in<extra></extra>",
        )
    )

    figure.update_xaxes(title_text=f"{labels[0]} (in)", constrain="domain")
    figure.update_yaxes(title_text=f"{labels[1]} (in)", scaleanchor="x", scaleratio=1.0)
    figure.update_layout(
        title=title,
        template="plotly_white",
        height=560,
        margin=dict(l=60, r=30, t=70, b=60),
        legend=dict(orientation="h", yanchor="bottom", y=1.0),
    )
    return figure


def end_on_view(vertices):
    """Looking straight down end leg 1-2, the view the drawing's end view uses."""
    axis = vertices[1] - vertices[0]
    axis = axis / np.linalg.norm(axis)
    reference = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(axis, reference)) > 0.9:
        reference = np.array([1.0, 0.0, 0.0])
    horizontal = np.cross(axis, reference)
    horizontal /= np.linalg.norm(horizontal)
    return "End view, looking along leg 1-2", (horizontal, np.cross(axis, horizontal), ("u", "w"))


def drawn_extent(centerline, radius, view):
    """Best-fit silhouette bounding box [in] over in-plane rotations.

    The drawing's orientation within its own view plane is arbitrary, so the
    comparable quantity is the box the silhouette can occupy, minimised
    against the extent measured off the scan.
    """
    horizontal, vertical, _ = view
    points = np.concatenate(tube_outline(centerline, radius, horizontal, vertical))
    best = None
    for degrees in np.arange(0.0, 180.0, 0.5):
        angle = math.radians(degrees)
        rotation = np.array(
            [[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]]
        )
        box = np.sort(np.ptp(points @ rotation.T, axis=0))[::-1]
        mismatch = float(
            abs(box[0] - DRAWING_END_VIEW_EXTENT[0]) + abs(box[1] - DRAWING_END_VIEW_EXTENT[1])
        )
        if best is None or mismatch < best[0]:
            best = (mismatch, box)
    return best[1], best[0]


def orthographic_views(frame):
    """The three standard views as ``(horizontal, vertical, axis labels)``."""
    x, y, z = frame
    return (
        ("Front view, looking along -y", (x, z, ("x", "z"))),
        ("Side view, looking along +x", (y, z, ("y", "z"))),
        ("Top view, looking along -z", (x, y, ("x", "y"))),
    )


def make_isometric_figure(centerline, vertices, title):
    """The resolved centerline in 3D."""
    figure = go.Figure()
    figure.add_trace(
        go.Scatter3d(
            x=centerline[:, 0],
            y=centerline[:, 1],
            z=centerline[:, 2],
            mode="lines",
            line=dict(color=COLORS[0], width=7),
            name="Centerline",
        )
    )
    figure.add_trace(
        go.Scatter3d(
            x=vertices[:, 0],
            y=vertices[:, 1],
            z=vertices[:, 2],
            mode="markers+text",
            marker=dict(color=COLORS[1], size=5),
            text=[str(index + 1) for index in range(len(vertices))],
            textposition="top center",
            name="Points 1-6",
        )
    )
    figure.update_layout(
        title=title,
        template="plotly_white",
        height=650,
        margin=dict(l=30, r=30, t=70, b=30),
        scene=dict(
            xaxis_title="x (in)",
            yaxis_title="y (in)",
            zaxis_title="z (in)",
            aspectmode="data",
        ),
    )
    return figure


def make_scene(pipe, result):
    """The solved duct as a shaded 3D tube, coloured by its own march.

    ``PipeComponent`` takes the geometry straight off the ``Pipe``. Each lumped
    ``RoundedBend`` is one connection between tangent points, so the solver
    centerline cuts its corner as a chord; the component reopens those into
    arcs for drawing using the ``center`` and ``normal`` the generator
    resolved. Only the drawn shape changes -- every length, angle, and loss
    still comes from the stations and connections.

    Scene coordinates are metres, unlike the inch-based drawing views.
    """
    scene = psf.viz.Scene(title="RL10 turbine-discharge duct", color_by="p_static")
    scene.add(
        psf.viz.PipeComponent(pipe, result, n_theta=48, segments=24),
        name="duct",
    )
    return scene


def make_pressure_figure(result, arc_result):
    """Static pressure along the duct for both connection models."""
    figure = go.Figure()
    for color, label, target in zip(
        COLORS,
        ("Lumped rounded bends", "Continuous curvature"),
        (result, arc_result),
    ):
        figure.add_trace(
            go.Scatter(
                x=target.s / INCH,
                y=target.p / 1.0e5,
                mode="lines",
                name=label,
                line=dict(color=color, width=2.5),
            )
        )
    figure.update_xaxes(title_text="Arc length (in)")
    figure.update_yaxes(title_text="Static pressure (bar)")
    figure.update_layout(
        title="Static pressure along the duct",
        template="plotly_white",
        height=520,
        margin=dict(l=70, r=30, t=70, b=60),
        legend=dict(orientation="h", yanchor="bottom", y=1.0),
    )
    return figure


def make_loss_figure(pipe, result):
    """Irreversible loss contributed by each connection."""
    labels = []
    bend_index = 0
    for index, connection in enumerate(pipe.connections):
        if isinstance(connection, psf.ducts.Straight):
            labels.append(f"straight {np.diff(result.s)[index] / INCH:.3f} in")
        else:
            radius, angle = BENDS[bend_index]
            labels.append(f"bend {angle:.2f} deg, R={radius:.3f} in")
            bend_index += 1

    figure = go.Figure(
        go.Bar(
            x=labels,
            y=result.pressure_loss / 1.0e3,
            marker_color=[
                COLORS[0] if isinstance(item, psf.ducts.Straight) else COLORS[1]
                for item in pipe.connections
            ],
            hovertemplate="%{x}<br>%{y:.2f} kPa<extra></extra>",
        )
    )
    figure.update_yaxes(title_text="Pressure loss (kPa)")
    figure.update_layout(
        title="Loss by connection, lumped rounded bends",
        template="plotly_white",
        height=520,
        margin=dict(l=70, r=30, t=70, b=140),
    )
    return figure


# --- solving and reporting ---------------------------------------------------


def solve(pipe, transport, inlet_state):
    """Run one loss march and return its result plus any correlation warnings."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = pipe.solve(
            inlet_state=inlet_state,
            mass_flow=MASS_FLOW,
            transport=transport,
            friction_model="fluids",
        )
    messages = sorted({str(entry.message) for entry in caught})
    return result, messages


def geometry_table(pipe, result, diameter):
    """Resolved geometry against the drawing and Table G1."""
    lengths = np.diff(result.s)
    total = lengths.sum() / INCH
    straight = sum(
        length
        for length, connection in zip(lengths, pipe.connections)
        if isinstance(connection, psf.ducts.Straight)
    ) / INCH
    return {
        "Bore diameter (in)": diameter / INCH,
        "Outer diameter, drawing (in)": OUTER_DIAMETER,
        "Wall thickness, drawing (in)": WALL_THICKNESS,
        "Roughness height (micro-in)": ROUGHNESS / INCH * 1.0e6,
        "Straight sections (in)": straight,
        "Straight sections, Figure G1 (in)": sum(STRAIGHTS),
        "Bend arcs (in)": total - straight,
        "Total centerline length (in)": total,
        "Actual duct length, Table G1 (in)": REFERENCE_DUCT_LENGTH,
        "Angle between end legs (deg)": end_leg_angle(diameter),
    }


def equivalent_lengths(pipe, result, diameter):
    """Each connection as an equivalent length of straight pipe [in].

    Table G1 reports each bend as the length of straight duct that would cause
    the same loss, so the comparable pyskyfire quantity is :math:`K D / f`
    using the straight-pipe Darcy factor at the same conditions.
    """
    friction = float(np.nanmean(result.darcy_factor))
    lengths = np.diff(result.s)
    rows = []
    total = 0.0
    bend_index = 0
    for index, connection in enumerate(pipe.connections):
        if isinstance(connection, psf.ducts.Straight):
            label = f"{lengths[index] / INCH:.3f} in straight"
            equivalent = lengths[index] / INCH
            references = (equivalent, equivalent)
        else:
            radius, angle = BENDS[bend_index]
            label = f"{angle:.2f} deg bend (R={radius:.3f})"
            equivalent = result.loss_coefficient[index] * diameter / friction / INCH
            references = REFERENCE_EQUIVALENT_LENGTHS[bend_index]
            bend_index += 1
        total += equivalent
        rows.append((label, equivalent, references[0], references[1]))
    rows.append(("total", total, *REFERENCE_TOTAL_EQUIVALENT_LENGTH))
    return rows, friction


def print_console_summary(pipe, result, messages, arc_result, arc_messages, diameter) -> None:
    """Keep the terminal output of the example self-contained."""
    print("RL10 turbine-discharge duct, Figure G1\n")
    print("Bore diameter")
    print(f"  OD 2.510/2.490                {OUTER_DIAMETER:8.4f} in")
    print(f"  wall .054/.044                {WALL_THICKNESS:8.4f} in")
    print(f"  OD less twice the wall        {OUTER_DIAMETER - 2 * WALL_THICKNESS:8.4f} in")
    print(f"  section A-A 2.408/2.396 AVG   {BORE_DIAMETER:8.4f} in")

    print("\nTangent distances, in")
    print(f"  {'bend':22s} {'R tan(angle/2)':>15s} {'figure G1':>10s}")
    for (radius, angle), tabulated in zip(BENDS, TABULATED_TANGENT_DISTANCES):
        resolved = radius * math.tan(math.radians(angle) / 2.0)
        print(f"  {f'{angle:.4f} deg, R={radius:.3f}':22s} {resolved:15.4f} {tabulated:10.3f}")

    print("\nResolved geometry")
    for key, value in geometry_table(pipe, result, diameter).items():
        print(f"  {key:38s} {value:10.4f}")

    print("\nBend-plane sign patterns (drawing rejects near-collinear end legs)")
    print(f"  drawing end view inks about {DRAWING_END_VIEW_EXTENT[0]} x "
          f"{DRAWING_END_VIEW_EXTENT[1]} in, good to an inch or two")
    for label, signs in CANDIDATE_ROLL_SIGNS.items():
        marker = "  <- used" if signs == ROLL_SIGNS else ""
        candidate_vertices = route_vertices(diameter, signs)
        _, view = end_on_view(candidate_vertices)
        box, mismatch = drawn_extent(dense_centerline(diameter, signs), diameter / INCH / 2.0, view)
        print(f"  {label:14s} {str(signs):12s} end legs {end_leg_angle(diameter, signs):7.2f} deg, "
              f"silhouette {box[0]:6.2f} x {box[1]:6.2f} in, mismatch {mismatch:5.2f} in{marker}")

    print(
        f"\nInlet hydrogen: {INLET_TEMPERATURE:.1f} K, {INLET_PRESSURE / 1.0e5:.1f} bar, "
        f"{MASS_FLOW:.2f} kg/s (assumed operating point)"
    )
    for label, target, notes in (
        ("Lumped rounded bends", result, messages),
        ("Continuous curvature, for contrast", arc_result, arc_messages),
    ):
        print(f"\n{label}")
        print(f"  total pressure loss      {target.total_pressure_loss / 1.0e3:8.2f} kPa"
              f"  ({target.total_pressure_loss / 6894.76:6.3f} psi)")
        print(f"  outlet pressure          {target.p[-1] / 1.0e5:8.2f} bar")
        print(f"  inlet velocity           {target.velocity[0]:8.2f} m/s")
        print(f"  Reynolds number          {target.reynolds.mean():8.3e}")
        for message in notes:
            print(f"  correlation warning: {message}")

    rows, friction = equivalent_lengths(pipe, result, diameter)
    print("\nEquivalent lengths, in")
    print(f"  {'section':28s} {'pyskyfire':>10s} {'manual':>8s} {'software':>9s}")
    for label, equivalent, manual, software in rows:
        print(f"  {label:28s} {equivalent:10.3f} {manual:8.3f} {software:9.3f}")
    print(f"  straight-pipe Darcy factor used: {friction:.5f}")


def build_report(pipe, result, messages, arc_result, arc_messages, diameter):
    """Assemble the duct report, one figure per block."""
    radius = diameter / INCH / 2.0
    vertices = route_vertices(diameter)
    centerline = dense_centerline(diameter)
    frame = part_frame(vertices)

    report = psf.viz.Report("RL10 turbine-discharge duct")

    tab = report.add_tab("Geometry")
    tab.add_markdown(
        "Figure G1 is the tube-bending table for the duct between the turbine "
        "housing discharge and the main fuel shut-off valve. Points 1 and 6 are "
        "the tube ends; points 2 to 5 are the bend vertices, each the "
        "intersection of the two straight control legs meeting there.\n\n"
        "The bore comes from the drawing rather than the loss study: OD "
        "2.510/2.490 less twice the 0.054/0.044 wall gives 2.402 in, which "
        "equals the 2.408/2.396 AVG DIA of Section A-A exactly."
    )
    tab.add_table(
        geometry_table(pipe, result, diameter),
        caption="Resolved geometry against Figure G1 and Table G1",
        key_title="Quantity",
        value_title="Value",
        precision=4,
    )
    tab.add_table(
        {
            f"{angle:.4f} deg bend, R={radius_inches:.3f} in": (
                f"R tan(angle/2) = {radius_inches * math.tan(math.radians(angle) / 2):.4f} in, "
                f"Figure G1 gives {tabulated:.3f} in"
            )
            for (radius_inches, angle), tabulated in zip(BENDS, TABULATED_TANGENT_DISTANCES)
        },
        caption="Figure G1's tangent-distance column reproduced from the angles and radii",
        key_title="Bend",
        value_title="Tangent distance",
    )
    tab.add_figure(
        make_isometric_figure(centerline, vertices, "Resolved duct centerline"),
        caption="Drag to rotate. Points 1-6 are the drawing's numbered vertices.",
    )

    tab = report.add_tab("3D view")
    tab.add_markdown(
        "The same `Pipe` object the solver marches, handed straight to "
        "`psf.viz.Scene`. The tube is swept along the resolved stations, with "
        "each lumped `RoundedBend` reopened into its arc for drawing; the loss "
        "model still sees only the two tangent stations and one fitting.\n\n"
        "Colour is the marched static pressure, so the four bends read as the "
        "steps they are. Scene coordinates are metres."
    )
    tab.add_iframe(
        make_scene(pipe, result).data_url,
        caption=(
            "Drag to orbit, scroll to zoom. Pressure falls from "
            f"{result.p[0] / 1e5:.1f} to {result.p[-1] / 1e5:.1f} bar along the duct."
        ),
    )

    tab = report.add_tab("Orthographic views")
    tab.add_markdown(
        "Three standard views of the resolved tube, drawn as silhouette walls "
        "with the centerline and the numbered points, for comparison against "
        "the fabrication drawing. The frame is part-intrinsic: **x** runs along "
        "the long end leg 5-6, **z** is the perpendicular component of the tall "
        "side of the hook, and **y** completes the set.\n\n"
        "The end faces are drawn as true projected circles, so the degree to "
        "which each one flattens into an ellipse is the same cue the drawing "
        "uses."
    )
    for title, view in (*orthographic_views(frame), end_on_view(vertices)):
        tab.add_figure(
            make_orthographic_figure(centerline, vertices, radius, frame, view, title)
        )

    tab = report.add_tab("Bend planes")
    tab.add_markdown(
        "Figure G1 tabulates the magnitude of each tube rotation "
        f"({ROLL_MAGNITUDES[0]:.3f}, {ROLL_MAGNITUDES[1]:.3f}, "
        f"{ROLL_MAGNITUDES[2]:.3f} deg) but not its sign. Of the four sign "
        "patterns, two leave the end legs near-collinear (20.4 and 159.6 deg "
        "apart). The drawing's end view shows the ends at a clear angle to each "
        "other -- one reads as a circle, the other as an ellipse -- so those two "
        "are rejected.\n\n"
        "The two survivors below are near-supplementary, so they foreshorten the "
        "far end by the same factor and give an equally convincing ellipse. The "
        "ellipse cannot separate them. The drawn silhouette leans toward "
        f"`{ROLL_SIGNS}`, but only by about two inches of fit against a scan "
        "measurement worth about that much, so this is **not settled**. Compare "
        "the views below against the drawing and set `ROLL_SIGNS` accordingly.\n\n"
        "**The choice does not change any pressure result** -- roll angles never "
        "enter the loss model."
    )
    extents = {}
    for label, signs in CANDIDATE_ROLL_SIGNS.items():
        candidate_vertices = route_vertices(diameter, signs)
        candidate_centerline = dense_centerline(diameter, signs)
        candidate_frame = part_frame(candidate_vertices)
        angle = end_leg_angle(diameter, signs)
        chosen = " (selected)" if signs == ROLL_SIGNS else ""
        view_title, view = end_on_view(candidate_vertices)
        box, mismatch = drawn_extent(candidate_centerline, radius, view)
        extents[f"{label} {signs}"] = (
            f"best silhouette {box[0]:.2f} x {box[1]:.2f} in, "
            f"mismatch {mismatch:.2f} in, end legs {angle:.1f} deg apart"
        )
        for title, candidate_view in ((view_title, view), *orthographic_views(candidate_frame)[:2]):
            tab.add_figure(
                make_orthographic_figure(
                    candidate_centerline,
                    candidate_vertices,
                    radius,
                    candidate_frame,
                    candidate_view,
                    f"{label}, signs {signs}, end legs {angle:.1f} deg apart{chosen}"
                    f"<br><sub>{title}</sub>",
                )
            )
    tab.add_table(
        extents,
        caption=(
            f"Silhouette fit against the drawing's end view, measured at roughly "
            f"{DRAWING_END_VIEW_EXTENT[0]} x {DRAWING_END_VIEW_EXTENT[1]} in"
        ),
        key_title="Candidate",
        value_title="Best in-plane fit",
    )

    tab = report.add_tab("Solution")
    tab.add_markdown(
        f"Representative RL10A-3-3A turbine-discharge hydrogen at "
        f"{INLET_TEMPERATURE:.0f} K, {INLET_PRESSURE / 1e5:.0f} bar, "
        f"{MASS_FLOW:.2f} kg/s. The operating point is the one assumed quantity "
        "here; the geometry and the roughness height are from the source.\n\n"
        "The lumped `RoundedBend` model is the correct one for discrete piping "
        "and is what Table G1 is doing. The `Continuous` trace applies the "
        "cooling-channel helical-coil curvature multiplier to the same arcs; it "
        "is shown only for contrast and runs outside its documented range."
    )
    tab.add_figure(make_pressure_figure(result, arc_result))
    tab.add_figure(make_loss_figure(pipe, result))

    summary = {
        "Total pressure loss (kPa)": result.total_pressure_loss / 1.0e3,
        "Total pressure loss (psi)": result.total_pressure_loss / 6894.76,
        "Outlet pressure (bar)": result.p[-1] / 1.0e5,
        "Inlet velocity (m/s)": result.velocity[0],
        "Outlet velocity (m/s)": result.velocity[-1],
        "Mean Reynolds number": result.reynolds.mean(),
        "Continuous-curvature loss (kPa)": arc_result.total_pressure_loss / 1.0e3,
    }
    tab.add_table(summary, caption="Loss march summary", key_title="Quantity", value_title="Value", precision=4)
    for message in messages + arc_messages:
        tab.add_markdown(f"Correlation warning: {message}")

    tab = report.add_tab("Equivalent lengths")
    tab.add_markdown(
        "Table G1 reports each section as the length of straight duct causing "
        "the same loss. The comparable pyskyfire quantity is `K D / f`, using "
        "the straight-pipe Darcy factor at the same conditions."
    )
    rows, friction = equivalent_lengths(pipe, result, diameter)
    tab.add_table(
        {
            label: f"pyskyfire {value:8.3f}   manual {manual:7.3f}   software {software:7.3f}"
            for label, value, manual, software in rows
        },
        caption=f"Equivalent lengths in inches, straight-pipe Darcy factor {friction:.5f}",
        key_title="Section",
        value_title="Equivalent length",
    )
    return report


def main(output_dir: Path | None = None, diameter_inches: float = BORE_DIAMETER) -> None:
    if output_dir is None:
        output_dir = Path(__file__).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    diameter = diameter_inches * INCH

    hydrogen = psf.thermo.Fluid(type="fuel", propellants=["Hydrogen"], fractions=[1.0])
    transport = psf.thermo.FluidTransport(hydrogen)
    resolver = psf.thermo.FluidStateResolver(transport, p_ref=INLET_PRESSURE)
    inlet_state = resolver.from_tp(INLET_TEMPERATURE, INLET_PRESSURE)

    pipe = build_pipe(diameter, rounded_bend)
    arc_pipe = build_pipe(diameter, lambda radius: psf.ducts.Continuous(radius, 40))

    result, messages = solve(pipe, transport, inlet_state)
    arc_result, arc_messages = solve(arc_pipe, transport, inlet_state)

    print_console_summary(pipe, result, messages, arc_result, arc_messages, diameter)

    report = build_report(pipe, result, messages, arc_result, arc_messages, diameter)
    report_path = output_dir / "rl10-duct-report.html"
    report.save_html(str(report_path))
    print(f"\nReport saved to {report_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, help="Directory for the generated report.")
    parser.add_argument(
        "--diameter",
        type=float,
        default=BORE_DIAMETER,
        help="Duct bore diameter in inches. Defaults to the drawing's 2.402 in bore.",
    )
    arguments = parser.parse_args()
    main(output_dir=arguments.output_dir, diameter_inches=arguments.diameter)
