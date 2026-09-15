"""Turn a fabrication bending table into an explicit ``Pipe`` centerline.

A tube on a fabrication drawing is specified the way the shop bends it: a run
of straight lengths measured tangent point to tangent point, a bend radius and
turn angle at each corner, and the dihedral angle from one bend plane to the
next. This module converts that table into ``psf.ducts.PathNode`` control
vertices, which is what ``generate_path_stations`` wants.

Two details make the conversion non-obvious, and both are handled here so each
component file only has to carry its own numbers:

Tangent distances
    A ``PathNode`` sits at the *intersection* of two straight control legs, but
    the table measures straights between tangent points. Each corner is
    therefore offset from its tangent points by :math:`R \\tan(\\theta/2)`,
    which is added to both neighbouring legs so the resolved straights come
    back out at their tabulated lengths. Drawings usually tabulate that tangent
    distance too, which makes it a free check on the radius and angle columns.

Bend planes
    The bend plane is carried by the bend's rotation axis, its binormal. Each
    tabulated tube rotation rolls that binormal about the intervening straight.
    Drawings give the magnitude of each rotation but not its sign; negating
    every sign produces the mirror image, which is a genuinely different part
    that a fixed pair of drawing views can distinguish even though it marches
    to the same pressure drop. Signs therefore belong with the component,
    resolved against its drawing views.
"""

import math

import numpy as np

import pyskyfire as psf


#: Metres per inch. Bending tables are almost always in inches.
INCH = 0.0254


def dms(degrees: float, minutes: float) -> float:
    """Degrees and arcminutes, as drawings state angles, in decimal degrees."""
    return degrees + minutes / 60.0


def rotate(vector, axis, angle):
    """Rotate ``vector`` about ``axis`` by ``angle`` [rad] (Rodrigues)."""
    axis = axis / np.linalg.norm(axis)
    return (
        vector * math.cos(angle)
        + np.cross(axis, vector) * math.sin(angle)
        + axis * np.dot(axis, vector) * (1.0 - math.cos(angle))
    )


def tangent_distances(bends):
    """``R tan(angle / 2)`` for each bend, in the radius' own units.

    Compare against the drawing's tangent-distance column to check the radius
    and angle columns that were read off it.
    """
    return tuple(
        radius * math.tan(math.radians(angle) / 2.0) for radius, angle in bends
    )


def route_from_bending_table(straights, bends, rolls, bend):
    """Build ``PathNode`` control vertices from a bending table.

    Parameters
    ----------
    straights : sequence of float
        Tangent-to-tangent straight lengths [m], one more than there are
        bends, ordered from the inlet end.
    bends : sequence of (float, float)
        Centerline bend radius [m] and turn angle [deg] at each corner.
    rolls : sequence of float
        Signed dihedral angle [deg] from one bend plane to the next, so one
        fewer than there are bends. Zero leaves two bends coplanar and turning
        the same way.
    bend : callable
        Called with a radius [m], returns the connection object attached to
        that corner, e.g. ``lambda r: psf.ducts.RoundedBend(radius=r)``.

    Returns
    -------
    list of psf.ducts.PathNode
        The drawing's numbered points, ends included.
    """
    if len(straights) != len(bends) + 1:
        raise ValueError("expected one more straight than there are bends")
    if len(rolls) != max(len(bends) - 1, 0):
        raise ValueError("expected one fewer roll than there are bends")

    position = np.zeros(3)
    tangent = np.array([1.0, 0.0, 0.0])
    binormal = np.array([0.0, 0.0, 1.0])

    nodes = [psf.ducts.PathNode(tuple(position))]
    leg = straights[0]

    for index, (radius, angle_degrees) in enumerate(bends):
        if index > 0:
            binormal = rotate(binormal, tangent, math.radians(rolls[index - 1]))
            binormal -= np.dot(binormal, tangent) * tangent
            binormal /= np.linalg.norm(binormal)

        angle = math.radians(angle_degrees)
        tangent_distance = radius * math.tan(0.5 * angle)

        position = position + (leg + tangent_distance) * tangent
        nodes.append(psf.ducts.PathNode(tuple(position), bend(radius)))

        tangent = rotate(tangent, binormal, angle)
        leg = tangent_distance + straights[index + 1]

    nodes.append(psf.ducts.PathNode(tuple(position + leg * tangent)))
    return nodes


def pipe_from_bending_table(straights, bends, rolls, diameter, bend, roughness):
    """Resolve a bending table straight into a ``Pipe``.

    Arguments are those of :func:`route_from_bending_table`, plus the
    roughness height [m] the pipe should carry.
    """
    nodes = route_from_bending_table(straights, bends, rolls, bend)
    stations, connections = psf.ducts.generate_path_stations(nodes)
    return psf.ducts.Pipe(stations, connections, diameter=diameter, roughness=roughness)
