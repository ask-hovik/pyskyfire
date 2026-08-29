# Rao Nozzle Contour Angles

A Rao-style bell nozzle replaces a longer conical nozzle with a curved contour.
The contour begins at the end of the throat blend with wall angle $\theta_n$
and reaches the nozzle exit with wall angle $\theta_e$. Both angles depend on
the nozzle expansion ratio and the selected nozzle length.

```{raw} html
<div class="psf-demo-frame psf-plot-demo">
  <iframe
    src="../_static/explanation-artifacts/rao-nozzle-angles.html"
    title="Interactive Rao nozzle contour-angle chart"
    loading="lazy"
    sandbox="allow-scripts allow-same-origin">
  </iframe>
</div>
```

The expansion ratio is

$$
\varepsilon = \frac{A_e}{A_t},
$$

where $A_e$ is the nozzle-exit area and $A_t$ is the throat area. The percentage
shown beside each curve is the bell-nozzle length expressed as a fraction of
the length of a 15-degree conical nozzle with the same expansion ratio. Shorter
nozzles save length and mass, but require different initial and exit angles to
form the bell contour.

Pyskyfire uses $\theta_n$ and $\theta_e$ as the endpoint slopes of a parabolic
bell. For values between the plotted curves, it first interpolates each angle
in expansion ratio and then interpolates between the two surrounding nozzle
length fractions.

## How the chart is generated

`tools/generate_engineering_charts.py` creates the chart with
`pyskyfire.viz.PlotThetaVsEpsilon`. The plotting class reads the digitised
$\theta_n$ and $\theta_e$ curves in
`src/pyskyfire/regen/data/theta_n.json` and
`src/pyskyfire/regen/data/theta_e.json`, then plots the data for nozzle length
fractions from 60% to 100%. These are the same data tables used by Pyskyfire's
Rao contour generator.

## Source

The angle curves were digitised from page 15 of *Liquid Rocket Engine Nozzles*,
NASA Space Vehicle Design Criteria (Chemical Propulsion), NASA SP-8120, July
1976. The original report is available from the [NASA Technical Reports
Server](https://ntrs.nasa.gov/citations/19770009165).
