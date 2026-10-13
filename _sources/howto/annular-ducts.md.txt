# Define and solve an annular duct

An `AnnularDuct` carries flow between a hub and a shroud that are both walls of
revolution about the machine's shaft. Its walls are curves in the meridional
`(x, r)` plane, so it can turn between axial and radial: a turbine's inlet or
exhaust annulus, or the crossover that takes a pump stage's radial discharge
round to the next stage's eye.

An annulus about a centreline that bends anywhere in space, such as an exhaust
with a centrebody that turns to one side, or coaxial lines, is a `Pipe` with a
`core_diameter`; see [Give a pipe a core](pipes.md#give-a-pipe-a-core). Both
are solved by the same march and return a `DuctResult`.

## A duct about the shaft

Give the mean line as `(x, r)` stations and the passage height across it. This
crossover runs radially outward, turns 180 degrees, and runs back in:

```python
import numpy as np
import pyskyfire as psf

R, r_in, r_top = 0.015, 0.06, 0.10
turn = np.linspace(0.0, np.pi, 25)
out = np.column_stack((np.zeros(10), np.linspace(r_in, r_top, 10)))
arc = np.column_stack((R - R*np.cos(turn), r_top + R*np.sin(turn)))
back = np.column_stack((np.full(10, 2*R), np.linspace(r_top, r_in, 10)))

duct = psf.ducts.AnnularDuct(np.vstack((out[:-1], arc, back[1:])), height=0.01)
```

Each span line is drawn through a mean-line station, normal to the mean line.
The flow area is the surface that line sweeps round the axis,
`pi*(r_hub + r_shroud)*height`. That is the area the through-flow crosses. It
stays exact with swirl, because the swirl runs along that surface and never
through it. The stations are also the march stations, so give enough of them
to resolve every turn.

`height` can be one value, one value per station, or a function of arc length
along the mean line. `duct.hub` and `duct.shroud` are the resulting walls, and
`duct.passage` is the `MeridionalPassage` they form. `duct.end("outlet")` gives
the `AnnulusEnd` that whatever follows has to match.

### Starting from walls

If you already have the walls, for example a turbine passage's, re-pair them:

```python
duct = psf.ducts.AnnularDuct.from_walls(hub, shroud, stations=41)
```

Any pairing between the points of the two walls is discarded. The mean line is
found by iteration: each interior span line is the chord cut between the walls
by the normal to the mean line, and the mean line is the chords' midpoints.
The two end span lines are the walls' own ends, left as given, because they are
where the duct joins its neighbours.

### Limits

The duct may not run backward along the shaft (toward -x). Its walls are named
hub and shroud by radius, and a backward stretch would swap them. A turn's
radius must also exceed half the passage height, or the inner wall would fold.

## Solve

It takes the same arguments as `Pipe.solve`, without the bend options:

```python
result = duct.solve(inlet_state=inlet, mass_flow=5.0, transport=transport,
                    inlet_swirl=12.0)
print(result.total_pressure_loss)
print(result.swirl_velocity[-1])
```

`inlet_swirl` is the circumferential velocity about the shaft at the inlet mean
radius [m/s].

The march is the one `Pipe` uses:

- It is adiabatic and marches stagnation pressure, less the irreversible
  losses.
- Density is taken at the stagnation state.
- `p_static` subtracts the local dynamic head, swirl included, as a
  first-order correction.

What it charges:

- **Wall friction** on the hydraulic diameter. The wall shear lies along the
  total velocity, so a swirling flow rubs along a helix that is longer than
  the duct by `c/c_m`. The loss therefore grows as `(c/c_m)**3` on the
  through-flow head.
- **Laminar friction** uses the exact concentric-annulus `f*Re` of the annulus
  the section resembles: the hub's for an axial duct, parallel plates for a
  radial one (`annulus_laminar_factor`). Turbulent friction is left
  uncorrected; on the hydraulic diameter it differs by a few percent.
- **Swirl** is carried as angular momentum `r*c_u`. A change of radius leaves
  it unchanged, because that is reversible: a crossover's return channel spins
  the flow up as it moves inward. Wall friction is what takes it away.
- **Meridional turns** are charged as friction times the curvature multiplier,
  like a `Pipe`'s `Continuous` connection. The turn presses the flow toward one
  wall equally all the way round, so no secondary flow runs round the
  circumference. Tighter than one hydraulic diameter, the turn separates and
  this underestimates its loss; the solver warns there with a
  `PipeCorrelationWarning`.
- **Changes of area** lose as a `Pipe`'s do (see
  [Losses where the area changes](pipes.md#losses-where-the-area-changes)): each run over which
  the area keeps growing is a diffuser and each over which it keeps shrinking
  a contraction, charged Rennels' conical correlation for the cone of the same
  area. A radial leg of constant height widens or narrows its area with its
  radius, so it is charged too. This is rougher for an annulus than for a
  pipe: how soon an annular diffuser separates depends on both walls' angles,
  not only on its area, and a radial passage behaves more like a vaneless
  diffuser than a cone. The swirl's own slowing as the radius grows is
  reversible and is not charged. `duct.cone_angle` gives the local angle of
  the equivalent cone, and `duct.area_change_loss_coefficient` the loss.

Every `loss_coefficient` is referred to the through-flow head `rho*c_m**2/2`.

## Draw

`viz.AnnularDuctComponent` draws the duct, optionally coloured by a result:

```python
scene = psf.viz.Scene()
scene.add(psf.viz.AnnularDuctComponent(duct, result), name="crossover")
```

The hub and shroud are revolved about the local +x shaft, the shroud
translucent so the hub shows.
