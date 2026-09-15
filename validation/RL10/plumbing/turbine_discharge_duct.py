"""The RL10 turbine-discharge duct, from Figure G1.

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
Figure G1 tabulates the dihedral angle between consecutive bend planes --
71 deg 06', 65 deg 39', 116 deg 35' -- but not the sign of each rotation.
Everything else about the shape is read directly off the table: the five
straight lengths, the four bend radii, and the four bend angles, cross-checked
against the table's own tangent-distance column.

``ROLL_SIGNS`` carries the signs, resolved by drawing each candidate layout and
comparing it against the fabrication drawing. Negating all three gives the
mirror image, which is a genuinely different part: a mirrored tube cannot be
rotated to match a fixed pair of drawing views, even though it marches to the
same pressure drop. Roll angles set the shape and nothing else -- the loss
model sees only lengths, angles, radii, and diameter.

Where this sits in the cycle
----------------------------
This is one section of the ``duct_turbine_injector`` leg in ``sizer_sim.py``,
not the whole leg: the main fuel shut-off valve and the injector inlet are
downstream of point 6 and are not modelled yet. The pipe is built and
available here, but nothing wires it into the engine network. Substituting it
for that leg's ``SimpleDuctBlock`` would trade a calibrated pressure ratio for
a physical model of part of the leg, which would move the cycle balance without
improving it. Retire the ratio once the rest of the leg exists.
"""

import pyskyfire as psf

from .bending_table import INCH, dms, pipe_from_bending_table, route_from_bending_table


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

#: Sign of each tabulated tube rotation, matched against the fabrication
#: drawing. See the module docstring.
ROLL_SIGNS = (1, -1, -1)

#: Signed rotations as the route builder wants them [deg].
ROLLS = tuple(sign * magnitude for sign, magnitude in zip(ROLL_SIGNS, ROLL_MAGNITUDES))

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


def rounded_bend(radius):
    """One lumped ``fluids`` rounded-bend loss between the tangent stations.

    This is the model to use for discrete piping, and it is what Table G1 is
    doing. The alternative, ``psf.ducts.Continuous``, is the cooling-channel
    helical-coil curvature multiplier over a discretized arc; it models
    curvature as enhanced wall friction rather than a discrete fitting, so it
    runs outside its documented range on a tube like this.
    """
    return psf.ducts.RoundedBend(radius=radius)


def _table(diameter):
    """The Figure G1 table converted to SI, ready for the route builders.

    Inches stop here. Every constant above transcribes a drawing, which is in
    inches; everything this module hands out is metres.
    """
    return (
        tuple(length * INCH for length in STRAIGHTS),
        tuple((radius * INCH, angle) for radius, angle in BENDS),
        ROLLS,
        diameter,
    )


def build_pipe(bend=rounded_bend, diameter=BORE_DIAMETER * INCH):
    """The duct as a solvable ``Pipe`` carrying the tabulated roughness.

    ``diameter`` is a bore in metres, defaulting to the drawing's 2.402 in.
    """
    return pipe_from_bending_table(*_table(diameter), bend, ROUGHNESS)


def build_route(bend=rounded_bend):
    """The duct's ``PathNode`` control vertices, drawing points 1 to 6 [m]."""
    return route_from_bending_table(*_table(BORE_DIAMETER * INCH)[:3], bend)
