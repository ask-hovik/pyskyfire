# Impeller scene and report

Install `pyskyfire[impeller,viz]`, then run:

```sh
python examples/impeller/impeller_report.py
```

Open `output/impeller-report.html`. The first tab embeds an interactive Scene;
the second shows dimensions, supplied efficiencies, the meridional passage and
the resolved inlet/outlet blade angles and wrap. The scene also exists as
`output/impeller-scene.html`. Sampled cylindrical arrays are saved in
`output/impeller-geometry.npz`.

The duty is a low specific speed one, chosen so splitters earn their place:
four full blades keep the small eye open, and a splitter in each passage
doubles the blade count over the outer half of the wheel.

Change the duty and design inputs in `build_impeller()`. `blade_thickness` can
be a float in metres or a function requesting any of `chord_fraction`,
`span_fraction`, `x`, `r`. `chord_fraction` is normalized conformal meridional
distance (m-prime) from the main-blade leading edge to trailing edge. Span runs
from the hub to the blade tip. Splitters follow the main camber and use the
same coordinate, beginning at their specified fraction.

Use `--screenshot` to also create a PNG through PyVista. `--output-dir PATH`
selects another folder.

The blade body thickness is fitted through pyTurbo controls. Edge closures taper
to zero; they are excluded from the body-thickness check. The generated surfaces
are suitable for inspection and mesh export. Hub and shroud are flow boundaries,
not a solid wheel with a shaft bore and manufacturing fillets.
