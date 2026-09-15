# Small Onshape hot-wall export

This geometry-only example uses real `Contour`, `CoolingCircuit`, and
`ThrustChamber` objects without a cooling solve or saved simulation. Edit the
`params` dictionary in [export_hot_wall.py](export_hot_wall.py), then run:

```bash
python examples/hot_wall_export/export_hot_wall.py
```

Defaults are 50 squared channels, a 15° helix, seven axial stations, and a
180 mm contour with a throat and bell. Geometry lengths are in metres and
angle settings are in degrees. Rounded cases use 2 mm **base height**, so their
full height includes additional lobes. These are CAD debugging dimensions.
Command-line arguments override only the settings explicitly supplied.

## Construction modes

- **Squared only:** revolve the contour over the covered circuit spans.
  Straight and helical channels use the same wall. In Onshape, enable
  **Cover entire hot wall** to replace these spans with one full-contour
  revolution, including uncovered ends and gaps. The checkbox defaults off.
- **Rounded or mixed:** loft the rounded/mixed spans curve by curve with
  tangent circumferential blends and preserved helical correspondence.
  Squared-only spans still use revolves. The full-contour checkbox is hidden.
- Default construction never bridges an uncovered axial gap. With rounded
  channels present, circuit entries/exits into shared spans leave an open
  zone centered on each boundary. `transition_fraction` (CLI
  `--transition-fraction`) sets its total width as a fraction of the full
  contour length. The default 0.02 makes a 2%-wide gap, with 1% on each side.
  The shared span between boundaries keeps its wall. Blend across the open
  cuts manually in CAD.

`trim_fraction` / `--trim-fraction` controls the tangent blends around a
profile: the fraction of each arc parameter trimmed per blended side, between
zero and one half (exclusive). It does not control axial gaps.

## Starter files

| File | Channels | Hot-wall construction |
| --- | --- | --- |
| [hot-wall-squared.json](hot-wall-squared.json) | 100 squared, straight | One revolve |
| [hot-wall-rounded.json](hot-wall-rounded.json) | 45 rounded, straight | 90 strip lofts |
| [hot-wall-rounded-50.json](hot-wall-rounded-50.json) | 50 rounded, straight | 100 strip lofts |
| [hot-wall-rounded-helical.json](hot-wall-rounded-helical.json) | 45 rounded, 15° helix | 90 strip lofts |
| [hot-wall-mixed.json](hot-wall-mixed.json) | 10 squared + 10 rounded, straight | 40 strip lofts |

These circuits cover the full contour and have no internal interleave sites.
Regenerate them using:

```bash
python examples/hot_wall_export/export_hot_wall.py --section squared --channels 100 --helix-angle-deg 0
python examples/hot_wall_export/export_hot_wall.py --section rounded --channels 45 --helix-angle-deg 0
python examples/hot_wall_export/export_hot_wall.py --section mixed --channels 20 --helix-angle-deg 0
python examples/hot_wall_export/export_hot_wall.py --section rounded --channels 50 \
    --helix-angle-deg 0 --output-name hot-wall-rounded-50.json
python examples/hot_wall_export/export_hot_wall.py --section rounded --channels 45 \
    --helix-angle-deg 15 --output-name hot-wall-rounded-helical.json
```

Use `--stations 3` for fewer coolant/loft profiles, `--output-dir PATH` for
another directory, or `--output-name NAME.json` for a new variant. Otherwise,
rerunning a section overwrites its corresponding JSON file. To exercise
uncovered ends, shorten `circuit_span` in `params`; the squared-only full-wall
option should still revolve the entire contour.

## Try in Onshape

Use the shared [FeatureScript](../advanced/onshape_cooling_channels.fs) and
version 8 JSON from the shared [exporter](../advanced/export_onshape.py).
**Older JSON must be regenerated**, and existing Feature Studios must be
updated to the current source.

The large checked-in advanced and RL10 exports still use older schemas. Their
local saved results reference the removed `pyskyfire.common.fluids` module;
generate fresh results before re-exporting those cases. The starter files
above are current and need no saved results.

1. Import a starter JSON into a fresh Part Studio and select it as **Channel
   data** in the PySkyFire cooling channels feature.
2. Enable **Construct hot wall**. The result is an inner surface open at its
   axial ends, not a wall solid with thickness.
3. For squared data, check that **Cover entire hot wall** appears and defaults
   off. For rounded and mixed data, check that it is hidden, including when
   changing the selected JSON from squared to rounded.
4. For lofts, **Knit each hot-wall span** defaults on. Turn it off to inspect
   individual strips. The helical case should follow the coolant channels.

**Show section planes** displays coolant section planes. FeatureScript checks
that loft sketches enclose one region and each exported curve resolves to one
edge. It reports failing curves and kernel errors in FeatureScript notices.
Successful Python export does not verify the Onshape CAD kernel.

See [the construction documentation](../../docs/explanations/channel-construction.md#hot-wall-export)
for the schema, coverage rules, and advanced exporter controls.
