# The Moody Diagram

The Moody diagram relates the Darcy friction factor, $f$, to Reynolds number
and relative surface roughness. Pyskyfire uses this friction factor when it
calculates pressure loss in a cooling channel.

```{raw} html
<div class="psf-demo-frame psf-plot-demo">
  <iframe
    src="../_static/explanation-artifacts/moody-diagram.html"
    title="Interactive Moody diagram"
    loading="lazy"
    sandbox="allow-scripts allow-same-origin">
  </iframe>
</div>
```

## Quantities in the diagram

The Reynolds number based on hydraulic diameter is

$$
Re_{D_h} = \frac{\rho u D_h}{\mu},
$$

where $\rho$ is the fluid density, $u$ is the mean channel velocity, $D_h$ is
the hydraulic diameter, and $\mu$ is the dynamic viscosity. Each curve in the
diagram represents a different relative roughness, $\epsilon/D_h$, where
$\epsilon$ is the absolute surface roughness.

For fully developed flow through a channel of length $L$, the friction factor
appears in the Darcy--Weisbach pressure-loss equation:

$$
\Delta p = f\frac{L}{D_h}\frac{\rho u^2}{2}.
$$

## Friction-factor model

Pyskyfire treats flow below $Re_{D_h}=2300$ as laminar and uses

$$
f_{lam} = \frac{64}{Re_{D_h}}.
$$

For turbulent flow with a specified roughness, it solves the Colebrook--White
equation iteratively:

$$
\frac{1}{\sqrt{f_{turb}}}
+ 2\log_{10}\left(
\frac{\epsilon}{3.71D_h}
+ \frac{2.51}{Re_{D_h}\sqrt{f_{turb}}}
\right) = 0.
$$

When no roughness is specified, the smooth-wall turbulent expression is

$$
f_{turb} = \left(0.79\ln Re_{D_h} - 1.64\right)^{-2}.
$$

Between $Re_{D_h}=2300$ and $Re_{D_h}=3500$, Pyskyfire uses a linear blend
rather than an abrupt switch. Defining

$$
\alpha = \frac{Re_{D_h}-2300}{3500-2300},
$$

the blended friction factor is

$$
f = (1-\alpha)f_{lam} + \alpha f_{turb}.
$$

Above $Re_{D_h}=3500$, the turbulent result is used directly.

## How the chart is generated

`tools/generate_engineering_charts.py` creates the diagram with
`pyskyfire.viz.PlotMoodyDiagram`. The plot evaluates the same
`pyskyfire.regen.f_darcy` function used by the cooling solver at 400
logarithmically spaced Reynolds numbers from $7\times10^2$ to $10^8$. It sets
$D_h=1$ so the absolute roughness value passed to the calculation is also the
relative roughness $\epsilon/D_h$. Each roughness curve is therefore a direct
visualisation of the implemented friction-factor model.
