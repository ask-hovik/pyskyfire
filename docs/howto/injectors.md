# Design and solve an injector

An injector in `pyskyfire.injector` is a list of element sets. Each set is one
element type repeated at face positions you supply. Holes are
`ducts.Orifice` objects, so sizing, discharge coefficients, choking and
cavitation work as in {doc}`Size and solve orifices <orifices>`.

The injector frame is the scene frame: the face is the plane x = 0, the
chamber lies towards +x, and a position `(r, θ)` on the face is
`y = r cos θ`, `z = r sin θ`.

## Impinging elements

An `ImpingingElement` is a set of `Jet`s aimed at one point a height above the
face. Each jet names its feed, its hole, its inclination from the chamber axis
and its azimuth around the element centre:

```python
import numpy as np
from pyskyfire.ducts import Orifice
from pyskyfire.injector import Jet, ImpingingElement

fuel_hole = Orifice.circular(1.35e-3, cd=0.75, length=6e-3)
ox_hole = Orifice.circular(1.80e-3, cd=0.75, length=7e-3)
film_hole = Orifice.circular(0.55e-3, cd=0.72, length=4e-3)

pair = ImpingingElement([
    Jet("main fuel", fuel_hole, inclination=np.radians(30), azimuth=0.0),
    Jet("oxidizer",  ox_hole,   inclination=np.radians(12), azimuth=np.pi),
], impingement_height=6e-3)

# One jet, aimed at the wall: its impingement point is where it strikes.
film = ImpingingElement([
    Jet("film ring", film_hole, inclination=np.radians(15), azimuth=np.pi),
], impingement_height=26e-3)
```

Like and unlike doublets, triplets and pentads are all the same class, with
different jets. Inclinations need not match: leaning the jet with more
momentum less sends the resultant down the axis.

A single jet is an element too. For film cooling, aim it shallowly at the
chamber wall, where it strikes and spreads along the wall as a film: centre
the element on the wall, and put the hole half a turn round (`azimuth=π`)
so it sits inboard and leans outward. The hole is then
`impingement_height·tan(inclination)` inside the wall, 7 mm here.

## Feeds

A feed is whatever supplies a hole. It can be feed hardware, such as a
`ducts.Manifold`, or a name where the hardware is not modelled yet. Feeds
compare as dictionary keys, so hardware objects compare by identity. Holes
on the same feed share one plenum; holes on different feeds, such as a film
ring supplied separately from the main plenum, never combine.

## Placing elements

Positions are plain arrays, from any source. Rings take a few lines of numpy:

```python
from pyskyfire.injector import ElementSet, Injector

counts, radii = [6, 12, 18, 24], [0.012, 0.024, 0.036, 0.048]
r = np.concatenate([np.full(n, x) for n, x in zip(counts, radii)])
theta = np.concatenate([2*np.pi*np.arange(n)/n for n in counts])

injector = Injector([
    ElementSet(pair, r, theta, orientation=0.0, name="main"),
    # Film elements centred on the 65 mm wall, where each jet strikes it.
    ElementSet(film, np.full(30, 0.065), 2*np.pi*np.arange(30)/30,
               orientation=0.0, name="film"),
])
```

`orientation` turns each element's azimuth zero from the local radial
direction towards the local tangential one.

## Sizing holes

Size each hole for its share of its feed's flow at the feed's pressure drop:

```python
fuel_hole = Orifice.sized(mass_flow=1.71/60, inlet_state=plenum, dp=9e5,
                          fluid="Methane", cd=0.75, length=6e-3)
```

## Solving

All the holes on one feed share its plenum and discharge to the same face, so
they run at one mass flux. A feed is exactly one orifice with the summed area
and `Cd·A` (`injector.equivalent_orifice(feed)`), and each hole takes the
share of the flow its own `Cd·A` sets.

`solve` takes one `FeedCondition` per feed. A condition with a pressure finds
the flow; one with a mass flow finds the plenum pressure. The two can be
mixed:

```python
from pyskyfire.injector import FeedCondition

result = injector.solve({
    "main fuel": FeedCondition("Methane", 120.0, pressure=69e5),
    "film ring": FeedCondition("Methane", 120.0, pressure=78e5),
    "oxidizer":  FeedCondition("Oxygen",  95.0,  mass_flow=6.0),
}, p_face=60e5, oxidizer=["oxidizer"])

result.mixture_ratio
result.feeds["film ring"].mass_flow
result.feeds["oxidizer"].plenum_pressure
result.sets["main"].diagnostics.rupe
```

`oxidizer` names the feeds carrying oxidizer; mixture ratios need it. Each
element type reports its own mixing quantities on `result.sets[name].diagnostics`.
For an impinging element these are the resultant momentum and its angle from
the axis, the element mixture ratio, and, for an element with one fuel and one
oxidizer jet, Rupe's ratio and the momentum-flux ratio. They are reported,
not enforced.

## Checks

```python
checks = injector.checks(chamber_radius=0.065)
checks.min_web, checks.overlaps, checks.wall_clearance, checks.target_clearance
```

These cover what the injector owns: the narrowest web between hole outlines
on the face, overlapping holes, the room between the holes and the wall
(`wall_clearance`), how far the impingement points sit inside it
(`target_clearance`, zero for jets aimed at the wall on purpose), and the
shortest hole relative to its diameter. Whether an
element sits where it should (clear of a baffle, say) is for whatever placed it.

## In an engine network

Each feed is one block, from its plenum station to a face station:

```python
blocks.append(psf.network.InjectorBlock(
    name="inj_film", st_in="film_ring", st_out="film_face",
    injector=injector, feed="film ring", fluid="Methane",
))
```

It emits `dp_<name>`, so pump `overcome` lists pick up the loss. The current
engine network prescribes flows, so the block answers the sizing question: the
drop each feed's flow causes. The split between two feeds from one source is
set by the network's splitter, not solved from the hydraulics.

## Viewing it

```python
from pyskyfire.viz import InjectorComponent, Scene, plot_injector_face

scene = Scene(title="Injector")
scene.add(InjectorComponent(injector, result=result, face_radius=0.065,
                            show_impingement=True), name="injector")
scene.save_html("injector.html")

plot_injector_face(injector, face_radius=0.065).show()
```

The scene draws wetted surface only: each hole's bore, one colour per feed,
ending where it meets the face plate. Hole terminations such as counterbores
or chamfers are not modelled yet. `show_impingement=True` adds grey arrows
from each hole to its impingement point; they are annotation, not hardware.
With a `result`, each bore carries `hole_mass_flow` and `jet_velocity` fields
to colour by.

`examples/injector/injector_report.py` puts this together: a LOX/methane
injector with a separately fed film ring, sized, rounded to drills, and
rated both ways.
