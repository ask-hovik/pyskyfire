# Axial turbine design and report

Install `pyskyfire[impeller,viz]`, then run:

```sh
python examples/turbine/turbine_report.py
```

Open `output/turbine-report.html`. The turbine is a single-stage impulse
oxygen turbine: 26 kg/s of oxygen at 700 K expanding from 220 bar to 140 bar
at 30 000 rpm. The nozzle takes the whole static enthalpy drop, and the
shrouded rotor only turns the flow. The numbers are the right size for an
oxygen-rich LOX turbopump turbine, not any one engine's data. The
efficiency is supplied, not predicted.

The report has three tabs:

- **Turbine**: the interactive 3D scene, with stator vanes grey and rotor
  blades red inside the translucent hub and casing. It also exists on its own
  as `output/turbine-scene.html`.
- **Design**: the duty, the meanline, the flow at every row edge and the
  velocity triangles.
- **Blade rows**: each row's hardware, with the throat and gauging angle
  measured on the built blades; the meridional passage with its rotating
  hub stretches; the mid-span sections unrolled; and metal angle, gauging
  angle, throat and pitch from hub to tip.

Change the duty and design choices in `build_turbine()`. Each kind of row
has its own `SectionProfile`. The nozzle turns the flow to nearly 80°, so
its chord leans 80 % of the way to the outlet angle (`stagger_fraction=0.8`)
before the built vanes' throats leave the flow where the triangles ask. The
rotor uses `SectionProfile.impulse`, which is thick at mid-chord and thin at
both edges, so its passage never narrows below its steep inlet. The Blade
rows tab checks, in one dimension, that neither row chokes.

Use `--screenshot` to also create a PNG through PyVista. `--output-dir PATH`
selects another folder.
