# Heat-exchanger examples

Run from the repository root:

```sh
python examples/heat_exchanger/baffles.py
python examples/heat_exchanger/two_fluid.py
```

Both accept `--nodes` (segments per straight path, default 24) and
`--output-dir` (default this directory's `output/`).

- **baffles.py:** two slightly inclined concentric rings forming a shallow V,
  plus six radial spokes with paired straight downflow/return banks. All
  sections are squared, with arced ring floors/roofs and rectangular spokes.
  CEA supplies prescribed chamber-gas recovery temperatures. The exterior
  coefficient is fixed at 1000 W/m²/K for a reproducible convergent example.
- **two_fluid.py:** hydrogen heats oxygen through explicitly paired segments
  in counterflow. Both streams use CoolProp and have independent inlet states.

Reports include the swept geometry, fluid states, heat budgets, wall interface
temperatures, and displaced volume. Return headers are ideal flow connections;
their turn geometry and fitting losses are not included. These examples test
the framework and convergence rather than establish a qualified baffle design.

See [the API guide](../../docs/howto/heat-exchangers.md) for path placement,
section dimensions, thermal links, pressure models, and current limitations.
