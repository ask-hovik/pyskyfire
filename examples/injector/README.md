# Impinging injector design and report

Install `pyskyfire[viz]`, then run:

```sh
python examples/injector/injector_report.py
```

Open `output/injector-report.html`. The injector is for LOX/methane at a
60 bar face: 6.0 kg/s of oxygen and 1.9 kg/s of methane. There are 60 unlike
doublets in four rings, fed from a main fuel plenum and an oxidizer plenum,
and an outer ring of 30 single fuel jets aimed shallowly at the chamber wall,
fed from a separate film ring at a higher pressure. Each film jet leans 15°
outward and strikes the wall about 26 mm downstream of the face, where it
spreads into a film. A film jet is an impinging element with one jet, whose
impingement point is on the wall. The numbers are the right size for an engine of roughly
20–25 kN, not any one engine's data.

Holes are sized for the design point with `Orifice.sized`, then rounded to the
nearest 0.05 mm drill. The report has four tabs:

- **Injector**: the interactive 3D scene of the wetted surface: each hole's
  bore, coloured by feed, ending at the face plate, with grey arrows from each
  hole to its impingement point. Colour it by *Hole mass flow* or
  *Jet velocity*. The face view below it shows each hole's outline and the
  jets projected onto the face. The scene also exists on its own as
  `output/injector-scene.html`.
- **Design**: the choices, the designed and drilled diameters, the mixing
  diagnostics and the layout checks.
- **As drilled**: the drilled injector answering two questions. What flows
  does it pass at the design plenum pressures? What plenum pressures does it
  need for the design flows?
- **Film ring**: the film share of the fuel as the film ring pressure varies.

Change the design point, stiffness and layout at the top of the script. The
unlike doublets lean their oxidizer jet less than their fuel jet, so the
resultant momentum leaves along the element axis.

Use `--screenshot` to also create a PNG through PyVista. `--output-dir PATH`
selects another folder.
