"""Turbine blade sections from pyTurbo 1.3.8, wrapped onto the annulus.

pyTurbo's ``Airfoil2D`` draws a cascade section flowing top to bottom (-y),
with the trailing edge at the origin and angles from axial. It turns the flow
toward +x. This package's angles are positive in the rotation direction, and
a row can turn the flow either way (a stator towards the rotation, a rotor
back against it), so each section is drawn turning the pyTurbo way and
mirrored where the row turns the other way.

Each section lies in the cascade plane unrolled at its own radius: axial
distance, and pitchwise distance ``r*theta``. The builder stacks the sections
on their centroids, then wraps each onto the annulus. A point's radius comes
from its axial position and span fraction, and its angle is its pitchwise
distance over the section radius. The throat and the narrowest gap are
measured in that unrolled plane, before wrapping.
"""

from __future__ import annotations

import math
import warnings
from importlib.metadata import PackageNotFoundError, version

import numpy as np

from pyskyfire.common.blade_surface import BladeSurface

from .geometry import RowGeometry, TurbineGeometry, _max_thickness, cascade_gaps

#: Leading-edge control offset per unit leading-edge diameter.
_LEADING_EDGE = 2/3
#: Section outline points used to measure the throat and thickness.
_MEASURE_POINTS = 400


def _require_pyturbo():
    try:
        installed = version("pyturbo-aero")
    except PackageNotFoundError as error:
        raise ImportError("Turbine blade generation uses pyTurbo: pip install "
                          "'pyskyfire[turbine]'") from error
    if installed != "1.3.8":
        raise ImportError(f"The turbine adapter supports pyturbo-aero==1.3.8; found {installed}. "
                          "Install 'pyskyfire[turbine]'.")
    return installed


def _section(row, inlet, outlet, stagger, n):
    """One section's suction and pressure sides as ``(axial, pitchwise)`` from its LE.

    Returns the two sides at ``n`` points and the outline at the measuring
    resolution.
    """
    from pyturbo.aero import Airfoil2D
    sign = 1.0 if outlet >= inlet else -1.0
    a = Airfoil2D(-sign*math.degrees(inlet), sign*math.degrees(outlet),
                  row.axial_chord, sign*math.degrees(stagger))
    chord = a.chord
    t = row.max_thickness/row.chord
    d = row.thickness_distribution
    a.add_le_thickness(_LEADING_EDGE*row.leading_edge_diameter/chord)
    a.add_ps_thickness([t*v for v in d.pressure], expansion_ratio=d.spacing,
                       camberPercent=d.pressure_extent)
    a.add_ss_thickness([t*v for v in d.suction], expansion_ratio=d.spacing,
                       camberPercent=d.suction_extent)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a.match_le_thickness()
    a.te_create(row.trailing_edge_thickness/2, *(math.degrees(w) for w in
                                                  row.trailing_edge_wedge_angles))
    y_le = a.cambBezierY[0]

    def plane(points):
        return np.column_stack((y_le - points[:, 1], sign*points[:, 0]))

    ss, ps = (plane(p) for p in a.get_points(n))
    ss_m, ps_m = (plane(p) for p in a.get_points(_MEASURE_POINTS))
    return ss, ps, np.vstack((ss_m, ps_m))


def _row(definition, row, span_points, chord_points):
    passage = definition.passage
    if row.tip_clearance > 0:
        passage = passage.with_clearance(row.tip_clearance, wall="shroud")
    if row.hub_clearance > 0:
        passage = passage.with_clearance(row.hub_clearance, wall="hub")
    hub_x, hub_r = passage.hub[:, 0], passage.hub[:, 1]
    shroud_x, shroud_r = passage.shroud[:, 0], passage.shroud[:, 1]

    def radius(x, f):
        h = np.interp(x, hub_x, hub_r)
        return h + f*(np.interp(x, shroud_x, shroud_r) - h)

    x_mid = 0.5*(row.x_le + row.x_te)
    spans = np.linspace(0.0, 1.0, span_points)
    sides = {"suction": [], "pressure": []}
    radii, inlets, outlets, throats, narrowest, where, thicknesses = [], [], [], [], [], [], []
    for f in spans:
        r_ref = float(radius(x_mid, f))
        inlet, outlet, stagger = row.section_angles(f, r_ref)
        ss, ps, outline = _section(row, inlet, outlet, stagger, chord_points)
        centroid = 0.5*(ss[:, 1].mean() + ps[:, 1].mean())
        pitch = 2*np.pi*r_ref/row.blade_count
        # The trailing edge: every outline point within a trailing-edge
        # thickness of its tip, which takes in the whole edge circle.
        tip = outline[_MEASURE_POINTS - 1]
        edge = outline[np.linalg.norm(outline - tip, axis=1) <= 1.5*row.trailing_edge_thickness]
        throat, gap, x_gap = cascade_gaps(outline, edge, pitch)
        throats.append(throat)
        narrowest.append(gap)
        where.append(min(max(x_gap/row.axial_chord, 0.0), 1.0))
        thicknesses.append(_max_thickness(outline[:_MEASURE_POINTS], outline[_MEASURE_POINTS:]))
        radii.append(r_ref)
        inlets.append(inlet)
        outlets.append(outlet)
        for name, side in (("suction", ss), ("pressure", ps)):
            x = row.x_le + side[:, 0]
            theta = definition.rotation_sign*(side[:, 1] - centroid)/r_ref
            sides[name].append(np.column_stack((x, radius(x, f), theta)))
    blade = BladeSurface(np.array(sides["suction"]), np.array(sides["pressure"]))
    return RowGeometry(row, blade, spans, radii, inlets, outlets, throats, thicknesses, narrowest,
                       where)


def build_geometry(definition, *, span_points, chord_points):
    for label, count, minimum in (("span_points", span_points, 2), ("chord_points", chord_points, 33)):
        if isinstance(count, bool) or not isinstance(count, (int, np.integer)) or count < minimum:
            raise ValueError(f"{label} must be an integer >= {minimum}")
    installed = _require_pyturbo()
    rows = tuple(_row(definition, row, span_points, chord_points) for row in definition.rows)
    return TurbineGeometry(definition.passage, rows, definition.rotation, installed)
