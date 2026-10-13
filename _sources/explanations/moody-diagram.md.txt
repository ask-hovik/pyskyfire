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

## Selecting the friction backend

This chart illustrates the **`fluids` backend for straight ducts**, including
pyskyfire's smooth transition policy. It is a representative Moody chart, not a
validation envelope or a plot of the native transition interpolation.
**`native` remains the default** for existing simulations and validation cases;
selecting `fluids` is explicit.

For regenerative cooling, pass `friction_model="fluids"`:

```python
import pyskyfire as psf

# Using an existing thrust chamber and coolant boundary conditions:
results = psf.regen.coupled_steady_heating_analysis(
    thrust_chamber,
    boundary_conditions,
    friction_model="fluids",
    pressure_curvature_correction=True,
)
```

`RegenBlock` also accepts `friction_model="fluids"` for network simulations.
At the scalar duct-physics level, the keyword is instead `model`:

```python
from pyskyfire.ducts import darcy_friction_factor

f = darcy_friction_factor(
    ReDh=1e5, Dh=0.004, roughness=4e-6, model="fluids",
)
```

Here roughness is an absolute length in metres, not a callable or a relative
roughness. Regen evaluates the circuit's local roughness before passing it to
the duct backend.

## Straight ducts: direct branches and smooth transition

The `fluids` backend calls the upstream branch functions directly. The laminar
branch is

$$
f_{lam} = \frac{64}{Re_{D_h}},
$$

and the turbulent branch is the value returned by Clamond for relative
roughness $\epsilon/D_h$. Clamond is a near-machine-precision numerical solver
for the Colebrook--White equation, including for smooth walls:

$$
\frac{1}{\sqrt{f_{turb}}}
+ 2\log_{10}\left(
\frac{\epsilon}{3.7D_h}
+ \frac{2.51}{Re_{D_h}\sqrt{f_{turb}}}
\right) = 0.
$$

Zero or `None` roughness means a smooth wall. Clamond is a Colebrook **solver**,
not an empirical explicit approximation such as Haaland or Swamee--Jain.

Pyskyfire uses the laminar branch through the straight-duct critical Reynolds
number $Re_{crit}=2040$. It preserves the relative width of its established
2300--3500 transition band by setting

$$
Re_{end}=Re_{crit}\frac{3500}{2300}
=3104.3478\ldots.
$$

Within that interval, define

$$
x=\frac{Re_{D_h}-Re_{crit}}{Re_{end}-Re_{crit}},
\qquad w=3x^2-2x^3,
$$

and blend the two direct branches as

$$
f=(1-w)f_{lam}+wf_{turb}.
$$

Thus $f=f_{lam}$ at and below $Re_{crit}$, $f=f_{turb}$ at and above
$Re_{end}$, and both value and blend weight slope join continuously at the
bounds. This transition width and cubic smoothstep are an explicit **pyskyfire
numerical/model policy**. They are not a transition width prescribed by the
upstream `fluids` library. The critical value likewise does not claim that
physical transition always occurs at exactly $Re=2040$.

The wide numerical accuracy range reported for Clamond is not an empirical
validation range for straight-pipe flow or for a rocket cooling channel. No
additional straight-pipe empirical limits are inferred from it.

## Curved ducts: Schmidt correlations

The straight-duct Moody chart does not include curvature. With
`pressure_curvature_correction=True`, the `fluids` backend calls the forced
Schmidt branch functions directly: `helical_laminar_fd_Schmidt` and
`helical_turbulent_fd_Schmidt`. It obtains the curved critical Reynolds number
from `helical_transition_Re_Schmidt`. Here $D_i=D_h$ is the
hydraulic-diameter substitution for the channel and $D_c=2|R|$ is the coil
diameter measured between tube centerlines. The sign of the centerline radius
does not affect friction. Schmidt's transition criterion is

$$
Re_{crit} = 2300\left[1+8.6\left(\frac{D_h}{2|R|}\right)^{0.45}\right].
$$

The same smoothstep policy is scaled from this geometry-dependent value:

$$
Re_{end}=Re_{crit}\frac{3500}{2300},
\qquad
x=\frac{Re-Re_{crit}}{Re_{end}-Re_{crit}},
\qquad w=3x^2-2x^3.
$$

The curved factor is the direct Schmidt laminar branch through $Re_{crit}$, the
direct Schmidt turbulent branch from $Re_{end}$ onward, and
$(1-w)f_{lam,curved}+wf_{turb,curved}$ between them. The $3500/2300$ scaling and
smoothstep are pyskyfire policy, not a transition width prescribed by Schmidt
or `fluids`. Consequently, the result matches the selected upstream Schmidt
branches exactly outside the transition interval; the blended result inside it
is intentionally pyskyfire's composite model, not the exact output of the
upstream generic curved-flow dispatcher.

The curvature multiplier is defined consistently as

$$
K_{curv} = \frac{f_{curved,blended}(Re,D_h,2|R|,\epsilon)}
{f_{straight,blended}(Re,D_h,\epsilon)},
\qquad f = f_{straight,blended}K_{curv} = f_{curved,blended}.
$$

**The same absolute roughness, and therefore the same relative roughness for a
given $D_h$, is used in the curved numerator and straight denominator.** This
avoids combining a rough straight factor with a smooth-wall curved ratio.
Upstream Schmidt laminar ignores roughness, so the final curved laminar factor
does too. Schmidt turbulent accepts the same absolute roughness and is used
unchanged in the fully turbulent region.

### Correlation warnings and applicability

`pyskyfire.ducts.FrictionCorrelationWarning` is the exported warning category
for finite curved evaluations outside these documented conditions:

| Warning condition | Meaning |
| --- | --- |
| $Re\leq100$ | Below the Schmidt laminar range, documented as $100<Re<Re_{crit}$. |
| $Re>150000$ | Above the documented Schmidt turbulent upper Reynolds number. |
| $D_h/(2|R|)\geq0.14$ | Outside the recommendation for the **Schmidt transition criterion**, $D_i/D_c<0.14$. |

The last condition is **not a known geometry range for the entire Schmidt
friction fit**. The calculation continues unmodified after a warning: there is
no clipping, fallback correlation, or suppression of curvature. These warnings
do not impose empirical limits on straight-duct Clamond evaluations.

These are continuously curved/helical-tube correlations, not discrete elbow
loss coefficients. Applying them locally to a noncircular cooling channel
through hydraulic diameter is a modelling assumption; numerical agreement with
the upstream library alone does not validate that assumption.

## Native backend retained for validation

The native backend keeps its existing equations: $64/Re$ below 2300, a linear
laminar/turbulent blend from 2300 to 3500, and a turbulent result above 3500.
For smooth walls it uses $(0.79\ln Re-1.64)^{-2}$; for nonzero roughness it
solves Colebrook with **3.71**, rather than the `fluids` value of 3.7. The native
rough-wall numerical solve uses bracketed `brentq` in $z=1/\sqrt{f}$ instead of
unconstrained `fsolve` in $f$:

$$
z + 2\log_{10}\left(\frac{\epsilon}{3.71D_h}
+ \frac{2.51z}{Re}\right)=0,\qquad f=z^{-2}.
$$

This changes the numerical solution method, not the native correlation or its
constants. Native curvature remains the RTE/Ito multiplier based on
$Re[(D_h/4)/|R|]^2$. Keeping native as the default preserves the backend choice
used by existing validation cases; the chart's use of `fluids` does not switch
those cases.

## How the chart is generated

`pyskyfire.viz.PlotMoodyDiagram` calls
`pyskyfire.ducts.darcy_friction_factor(..., model="fluids")`. By default it
samples 400 logarithmically spaced Reynolds numbers from $7\times10^2$ to
$10^8$ and explicitly includes both straight transition bounds, 2040 and
$2040(3500/2300)\approx3104.35$. There is no NaN plotting gap: each roughness
curve continuously displays the implemented smooth blend. The chart sets
$D_h=1$ so the absolute roughness value equals the relative roughness
$\epsilon/D_h$. The displayed Reynolds-number span is a plotting choice, not a
validation range.

Generate **only** the Moody HTML artifact from the repository root:

```sh
python tools/generate_engineering_charts.py --chart moody
```

This writes `docs/_static/explanation-artifacts/moody-diagram.html` without
regenerating the Rao nozzle chart. The HTML uses CDN-hosted Plotly and MathJax,
so viewing it requires network access.

See the upstream [`fluids` friction documentation](https://fluids.readthedocs.io/fluids.friction.html)
for Clamond, Schmidt's laminar and turbulent correlations, and the separate
Schmidt transition recommendation. See also the
[coolant pressure-marching explanation](regenerative-cooling-explanation.md)
for how the selected factor enters the cooling solver.
