# Small Onshape hot-wall export

A geometry-only debug engine with editable settings in the `params` dictionary
near the top of [export_hot_wall.py](export_hot_wall.py). It uses real
`Contour`, `CoolingCircuit`, and `ThrustChamber` objects but runs no CEA
calculation, cooling analysis, or visualiser. No saved simulation is needed.

The current settings use **45 rounded channels, 7 axial stations, and a 180 mm
contour**, with a constant 30 mm radius so every station carries the same
profile. That isolates the CAD loft from contour variation. Switch
`contour_rs` to the commented-out varying array for a throat and a bell.

Channel counts here are for CAD debugging, not a realistic thermal design.
Rounded tubes are deliberately large; their 2 mm **base height** excludes the
lobes, so their full height is larger than the squared channels' 2 mm height.

Run from the repository root with your PySkyFire Python environment:

```bash
python examples/hot_wall_export/export_hot_wall.py
```

Edit `params` directly, then run the file normally (or use your IDE's Run
button):

```python
params = dict(
    section="rounded",   # Or "squared" / "mixed"
    n_channels=45,
    n_stations=7,
    helix_angle_deg=0.0,
    # ... keep the other entries in the existing dictionary
)
```

The same dictionary contains channel height, rib width, wall thickness,
material, circuit span, contour coordinates, hot-wall trim/transition settings,
and output location. Geometry lengths are in metres; the angle entries are in
degrees. Command-line options only override values explicitly supplied, so they
will not silently replace your dictionary edits.

Hot-wall data is always included. The generated starter files are:

| File | Channels | Coolant loft seeds | Hot-wall curves per profile | Strip lofts |
| --- | --- | --- | --- | --- |
| [hot-wall-squared.json](hot-wall-squared.json) | 100 squared | 1 | 100 | 100 |
| [hot-wall-rounded.json](hot-wall-rounded.json) | 45 rounded | 1 | 90 | 90 |
| [hot-wall-rounded-50.json](hot-wall-rounded-50.json) | 50 rounded | 1 | 100 | 100 |
| [hot-wall-rounded-helical.json](hot-wall-rounded-helical.json) | 45 rounded, 15° helix | 1 | 90 | 90 |
| [hot-wall-mixed.json](hot-wall-mixed.json) | 10 squared + 10 rounded | 2 | 40 | 40 |

All circuits run the full length. There are no circuit hand-offs, axial gaps,
liners, or hot-side ribs in these cases, so each file is a single run of
constant profile topology and needs no transition bridge. The squared case has
no Bézier blends; the rounded and mixed cases add the tangent connectors.

## How the hot wall is lofted

A hot-wall profile is one curve per squared channel, or an arc plus a tangent
blend per rounded channel. Those curves are alike and evenly spaced around the
axis, which is the worst possible input for a loft that pairs profiles by
proximity: it fails outright past roughly ninety curves, and a helix that
advances more than half a channel pitch between stations pairs each curve with
its neighbour, twisting the wall against the channels.

The export therefore carries `hotWall.runs`, spans of consecutive sections
whose curves are the same channels in the same order. FeatureScript lofts each
curve along its run on its own, so every loft profile is a single edge and
there is nothing left to pair by guesswork. The strips are then knitted into
one surface body; if the kernel declines to knit them the feature keeps the
strips and says so in a notice rather than failing.

One strip is one loft, so these files are 40–100 operations each and take a few
seconds to regenerate. A full engine export is far larger — see
`--circuit` in [the advanced example](../advanced/export_onshape.py).

## Try in Onshape

Use the current shared [FeatureScript](../advanced/onshape_cooling_channels.fs).
This example calls the same [exporter](../advanced/export_onshape.py) as the
advanced example and writes the same version 6 JSON schema.

1. Start with `hot-wall-squared.json` in a fresh Part Studio.
2. Select that JSON as **Channel data** in the PySkyFire cooling channels
   feature. Leave **Construct hot wall** off initially to check the coolant
   loft alone.
3. Enable **Construct hot wall**. The expected addition is one connected inner
   surface, open at the two axial ends. It is not a wall solid with thickness.
   **Knit hot wall into one surface** appears with it; turn it off to inspect
   the individual strips.
4. Test `hot-wall-rounded.json`, then `hot-wall-rounded-50.json`, then
   `hot-wall-mixed.json` in separate fresh Part Studios to add the ellipse and
   Bézier geometry and to pass the count that used to fail.
5. Test `hot-wall-rounded-helical.json` last. The hot wall must twist the same
   way as the coolant channels.

**Show section planes** displays coolant section planes. If a case fails, note
which JSON was used, whether the hot-wall checkbox was enabled, and the exact
Onshape error. The feature checks that every sketch encloses exactly one region
and that every exported curve became exactly one sketch edge, and it names the
failing curve and run when a strip loft fails. Kernel errors are also printed
to the FeatureScript notices. Successful Python export does not verify the
Onshape loft kernel.

If updating from an earlier FeatureScript, replace the Feature Studio source
with the current shared file and update the custom feature reference if it uses
a version. **Version 4 and 5 JSON files must be regenerated**: version 5 added
the `hotWall.runs` block the strip lofts read, and version 6 the
`hotWall.bridges` block that pairs the profiles either side of a hand-off.

These single-circuit files have no hand-off, so none of them exercise a bridge.
The [advanced example](../advanced/onshape-cooling-channels.json) does: it hands
off from 140 rounded channels to 50, then picks up a squared return pass.

## Change one thing at a time

Regenerate the starter variants:

```bash
python examples/hot_wall_export/export_hot_wall.py --section squared --channels 100
python examples/hot_wall_export/export_hot_wall.py --section rounded --channels 45
python examples/hot_wall_export/export_hot_wall.py --section mixed --channels 20
python examples/hot_wall_export/export_hot_wall.py --section rounded --channels 50 \
    --output-name hot-wall-rounded-50.json
python examples/hot_wall_export/export_hot_wall.py --section rounded --channels 45 \
    --helix-angle-deg 15 --output-name hot-wall-rounded-helical.json
```

Reduce the baseline further to just three axial profiles:

```bash
python examples/hot_wall_export/export_hot_wall.py --stations 3 \
    --output-name hot-wall-three-stations.json
```

Other controls:

- `--channels 45`: total count, at least four; mixed cases require an even count.
- `--stations 7`: profiles per circuit, at least three. The default includes the
  throat; use an odd number to retain it on this symmetric axial grid.
- `--trim-fraction 0.1`: fraction of each arc's parameter span trimmed from each
  blended side, between zero and one half (exclusive).
- `--output-dir PATH`: directory for JSON output.
- `--output-name NAME.json`: filename override. Otherwise the filename depends
  only on `--section`, so rerunning a section overwrites its previous export.

The script prints the output size and sketch counts. All example settings are in
`params` near the top of the script.
