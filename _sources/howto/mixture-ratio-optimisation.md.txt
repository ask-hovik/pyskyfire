# Optimise mixture ratio

This how-to finds the oxidizer-to-fuel mixture ratio that gives the highest vacuum specific impulse for a fixed engine design point.

The example uses:

| Parameter | Value |
|---|---:|
| Fuel | Methane |
| Oxidizer | Oxygen |
| Chamber pressure | 100 bar |
| Thrust | 100 kN |
| Nozzle area ratio | 60 |
| Characteristic length | 1.2 m |

## Objective

The mixture ratio is defined as:

$$
MR = \frac{\dot{m}_{ox}}{\dot{m}_{fu}}
$$

For each trial mixture ratio, Pyskyfire builds a new `Aerothermodynamics` object and reads the calculated vacuum specific impulse.

The optimisation problem is:

$$
\max_{MR} I_{sp,vac}(MR)
$$

The script solves this by minimizing the negative specific impulse:

$$
\min_{MR} -I_{sp,vac}(MR)
$$

## Run the optimisation

Run the full script at: [`examples/MR_optimisation/MR_opt.py`](https://github.com/ask-hovik/pyskyfire/blob/main/examples/MR_optimisation/MR_opt.py)

The central optimisation step is:

```{literalinclude} ../../examples/MR_optimisation/MR_opt.py
:language: python
:start-after: tutorial:start:optimisation
:end-before: tutorial:end:optimisation
:dedent: 4
```

The script uses `scipy.optimize.minimize_scalar` with bounded optimisation over the interval:

$$
1.5 \le MR \le 6.0
$$

## Optimisation curve

The curve below shows the vacuum specific impulse calculated over the mixture-ratio range. The marker shows the optimum found by the optimiser.

```{raw} html
<div class="psf-wide-frame">
  <iframe
    src="../_static/howto-artifacts/mixture-ratio-optimization/mixture-ratio-optimisation.html"
    title="Mixture-ratio optimisation curve"
    loading="lazy"
    sandbox="allow-scripts allow-same-origin">
  </iframe>
</div>
```

## Optimised contour

As a curiosity, after finding the optimum mixture ratio, we can easily find and plot a corresponding contour. 

```{literalinclude} ../../examples/MR_optimisation/MR_opt.py
:language: python
:start-after: tutorial:start:contour
:end-before: tutorial:end:contour
:dedent: 4
```

Yielding: 

```{raw} html
<div class="psf-wide-frame">
  <iframe
    src="../_static/howto-artifacts/mixture-ratio-optimization/optimized-contour.html"
    title="Optimised thrust-chamber contour"
    loading="lazy"
    sandbox="allow-scripts allow-same-origin">
  </iframe>
</div>
```

## Interpretation

This optimisation only varies mixture ratio. Chamber pressure, thrust, area ratio, characteristic length, propellant choice, and inlet temperatures are fixed. The result is therefore the mixture ratio that maximises ideal vacuum performance for this design point.