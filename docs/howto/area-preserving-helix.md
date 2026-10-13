# Hold the coolant area constant with a helix

A squared channel of fixed height, rib width and count gets wider wherever the
wall it wraps gets wider. On an engine with area ratio 8, an axial channel
built this way has more than three times the throat area at the nozzle exit. To
keep the area constant, wind the channel into a helix that steepens as the
radius grows. The extra circumference then goes into the winding instead of the
channel width.

This guide finds that helix by root-finding the helix angle against the area
Pyskyfire computes. The result below has 40 channels, each holding the same
coolant area to within one part in a million from the injector to the nozzle
exit. Drag to rotate, scroll to zoom and right-drag to pan.

```{raw} html
<div class="psf-demo-frame">
  <iframe
    src="../_static/howto-artifacts/area-preserving-helix/engine-3d.html"
    title="Thrust chamber with area-preserving helical channels"
    loading="lazy"
    sandbox="allow-scripts allow-same-origin">
  </iframe>
</div>
```

The complete script is `examples/area_preserving_helix/area_preserving_helix.py`.

## Why the helix opens out with radius

Consider $n$ squared channels of height $h$ whose floor sits at radius $r_b$.
Each channel and its share of the rib fill an annular sector. On a helix of
angle $\gamma$, measured from the meridian, that sector's angle measured
normal to the channel is $2\pi\cos\gamma/n$. Subtracting a rib of width
$2s$, where $s$ is the `sidewall_half_thickness`, leaves approximately

$$
A \approx h\left[\frac{\pi\cos\gamma}{n}\left(2r_b + h\right) - 2s\right].
$$

With $h$, $s$ and $n$ fixed, $A$ is constant when $\cos\gamma\,(2r_b + h)$ is.
Taking a station where the floor radius is $r_{b,0}$ and the helix angle is
$\gamma_0$,

$$
\cos\gamma = \cos\gamma_0\,\frac{2r_{b,0} + h}{2r_b + h}.
$$

As the radius grows, $\cos\gamma$ falls and the helix opens out. Since
$\cos\gamma \le 1$, no station can have a smaller radius than the reference.
The reference is therefore the narrowest point, which is usually the throat.

This formula treats the wall as a cylinder. It misses two effects that grow
with the helix angle. Both rib faces are helicoids, so the rib is wider below
the roof than its nominal width (see {doc}`../explanations/channel-construction`).
On a sloping wall, the wall stack stands off along the normal, so $r_b$ is
larger than $r + t_\text{wall}$. Used directly as a helix law on the engine
above, the formula undershoots the target area by up to 10 % on the convergent,
where the helix has reached 55°.

The root finder below matches the area Pyskyfire actually computes, so it
includes both effects without modelling them.

## Set up the geometry

Write a function that builds the chamber for any helix law, and a helix law
that interpolates $\cos\gamma$ from a table on the solver stations. A helix law
is a design input that requests quantities by parameter name (see
{doc}`../explanations/design-input-functions`), so its argument must be called
`x`.

```{literalinclude} ../../examples/area_preserving_helix/area_preserving_helix.py
:language: python
:start-after: tutorial:start:geometry
:end-before: tutorial:end:geometry
:dedent: 4
```

The transports are `None` because only the geometry is needed. In a real case,
pass the transports you already have.

## Root-find the helix angle

The area at each station depends only on the helix angle at that station. One
build therefore evaluates the residual at every station, and a secant step
applied to the whole array solves them all together:

```{literalinclude} ../../examples/area_preserving_helix/area_preserving_helix.py
:language: python
:start-after: tutorial:start:root-find
:end-before: tutorial:end:root-find
:dedent: 4
```

```text
max area error: 2.8e-02
max area error: 6.2e-04
max area error: 5.4e-07
```

Five builds in total, each well under a second.

- **The target** is the largest area every station can reach, which is the
  axial channel's area at its narrowest point. You can also set `A_target`
  directly. Stations whose axial area falls short of the target remain axial,
  because $\cos\gamma$ is clipped at 1.
- **Iterate on $\cos\gamma$, not $\gamma$.** The area is nearly linear in
  $\cos\gamma$, so the secant step lands close to the root. Stepping in
  $\gamma$ can overshoot to an angle where the rib fills the whole pitch, and
  the build then raises `sidewall_half_thickness leaves no room for the channel`.

## Build the 3D view

Build the chamber once more with the converged law and save the scene. The
chamber uses 400 stations and `stride=1` draws every one of them, so the helix
renders cleanly where it is steepest:

```{literalinclude} ../../examples/area_preserving_helix/area_preserving_helix.py
:language: python
:start-after: tutorial:start:scene
:end-before: tutorial:end:scene
:dedent: 4
```

```text
target area 6.420 mm², helix angle 0.0° to 68.0°
```

## Limits

- The tabulated law depends only on `x`, so several circuits can share the
  placement. It is tied to the contour, wall stack, channel height, channel
  count and `n_nodes` it was solved against. Solve it again whenever any of
  these change.
- A change in channel count, such as under an overlapping return pass, changes
  the pitch abruptly. No helix law can hold the area across that step.
- The independence of stations was checked for a single circuit on a
  `SurfacePlacement`. When several circuits are interleaved, keep the
  convergence printout: if a station's area also depends on its neighbours,
  convergence slows down rather than failing silently.
- Large helix angles lengthen the passage by roughly $1/\cos\gamma$, which
  increases the pressure drop. At 68° the channel is about 2.7 times longer
  than an axial one. Check the pressure drop in the coupled analysis, not just
  the area.
