# Volute design study

Designs three volutes, compares their losses, and renders all three in one 3D
scene.

```
python examples/volute/volute_report.py --output-dir examples/volute/output
```

Add `--screenshot` to also write a PyVista PNG preview. Requires
`pyskyfire[viz]`; no browser or display is needed for the HTML output.

## What it builds

| Case | Inlet | Notes |
|---|---|---|
| LOX pump discharge | radial | Dense propellant, shallow inlet angle, friction-led loss |
| Methane pump discharge | radial | Same duty class at a third of the density, steeper inlet |
| Coolant inlet manifold | inclined 45° | Large radius, modest swirl, steepest inlet of the three |

All three are sized on conservation of angular momentum at a sizing parameter
of 1, which is the condition that zeroes the tangential-velocity dump loss.

## Outputs

* `volute-report.html` — tabs for the scenes, design inputs, geometry and losses
* `volute-scene-<case>.html` — one standalone 3D scene per volute
* `volute-preview-<case>.png` — one per volute, with `--screenshot`

Each volute gets its own scene rather than sharing one, so each keeps its own
colour scale; they differ in through-flow velocity by about a factor of three.
