# Specific Impulse and Thrust Coefficient

An `Aerothermodynamics` design point reports specific impulse under four
different names: `Isp_vac`, `Isp_optimum`, `Isp_amb`, and `Isp_SL`. They can
differ by more than an order of magnitude for the same engine. This page
explains what each one assumes, and when the ambient values stop being
meaningful.

## Where the differences come from

The thrust of a rocket nozzle has two parts:

$$
F = \dot{m}v_e + \left(p_e - p_{amb}\right)A_e,
$$

where $\dot{m}$ is the propellant mass flow, $v_e$ is the exhaust velocity at
the nozzle exit, $p_e$ is the static pressure at the exit plane, $A_e$ is the
exit area, and $p_{amb}$ is the pressure of the surrounding atmosphere.

The first term is fixed by the combustion and expansion process. The second
term is pure geometry and pressure: it depends on where the engine is flying,
not on how well it burns. Every specific impulse below is the same first term
combined with a different assumption about the second.

Specific impulse follows from thrust as

$$
I_{sp} = \frac{F}{\dot{m}g_0},
$$

with $g_0 = 9.81\ \mathrm{m\,s^{-2}}$. It is convenient to separate the
contributions of the combustion chamber and the nozzle by writing thrust as

$$
F = C_F p_c A_t,
$$

where $p_c$ is the chamber pressure and $A_t$ the throat area. The
*characteristic velocity* $c^*$ measures the chamber, and the *thrust
coefficient* $C_F$ measures the nozzle:

$$
c^* = \frac{p_c A_t}{\dot{m}},
\qquad
I_{sp} = \frac{C_F c^*}{g_0}.
$$

Pyskyfire computes one $c^*$ and one vacuum thrust coefficient, then obtains
every other operating condition by subtracting the ambient pressure term:

$$
C_{F,amb} = C_{F,vac} - \frac{p_{amb}}{p_c}\varepsilon,
\qquad
\varepsilon = \frac{A_e}{A_t}.
$$

## The four quantities

`Isp_vac` is the vacuum specific impulse, $p_{amb} = 0$. The full pressure
term is added, so this is the largest of the four. It is the value to use for
upper stages and in-space engines, and it is what the mixture-ratio optimiser
in {doc}`Mixture ratio optimisation <../howto/mixture-ratio-optimisation>`
maximises.

`Isp_optimum` is the specific impulse at *optimum expansion*, meaning the
condition $p_e = p_{amb}$ where the nozzle exit pressure exactly matches the
surrounding atmosphere. The pressure term vanishes and only $\dot{m}v_e$
remains, so this quantity is simply $v_e/g_0$. It is the number NASA CEA
reports as `Isp`.

Two properties of `Isp_optimum` are easy to misread. It does not depend on the
`p_amb` given to the constructor — CEA never receives that value. And it is not
tied to a fixed altitude: it floats with whatever exit pressure the expansion
ratio happens to produce, so it describes a different altitude for every
nozzle.

`Isp_amb` is the specific impulse at the ambient pressure supplied as `p_amb`
when the design point is constructed. `Isp_SL` is the same calculation at a
fixed sea-level pressure of 101325 Pa. Both subtract the pressure term above,
and both carry the validity limit described in the next section.

The thrust coefficients `CF_vac`, `CF_amb`, and `CF_SL` correspond one-to-one
with these, and are related to them by $I_{sp} = C_F c^*/g_0$.

## A worked example

The RL10A-3-3A validation case in `validation/RL10/` is a hydrogen/oxygen
upper-stage engine with $\varepsilon = 61$, $p_c = 33.2$ bar, and
$c^* = 2375\ \mathrm{m\,s^{-1}}$, constructed with `p_amb=1e5`:

| Quantity | $C_F$ | $I_{sp}$ | Assumed $p_{amb}$ |
| --- | --- | --- | --- |
| Vacuum | 1.937 | 469.1 s | 0 |
| Optimum expansion | — | 452.0 s | $p_e = 3.85$ kPa |
| Ambient | 0.102 | 24.7 s | 100 kPa |
| Sea level | 0.078 | 18.8 s | 101.325 kPa |

The gap between the vacuum and optimum values, 17.1 s, is exactly the pressure
thrust $p_e A_e/(\dot{m}g_0)$. The exit pressure of 3.85 kPa corresponds to an
altitude near 22 km, which is the altitude `Isp_optimum` describes for this
particular nozzle.

The sea-level figure is a different matter. It is not a small correction to
469 s; it is 4% of it.

## When the ambient values stop being valid

The expression for $C_{F,amb}$ assumes the nozzle flows full, with the exhaust
attached to the wall all the way to the exit plane. A nozzle sized for vacuum
violates that assumption at sea level.

In the example above the exit pressure is 3.85 kPa against an ambient of
101 kPa, an overexpansion of a factor of 26. Real flow does not tolerate this.
The exhaust separates from the wall well upstream of the exit, at roughly the
station where the wall pressure falls to $0.3$–$0.4$ times ambient, which for
this gas is somewhere near $\varepsilon = 10$ rather than 61. The engine then
behaves like a much smaller nozzle, delivers substantially more than the
tabulated 18.8 s, and generates severe lateral loads on the nozzle extension.

Engine practice reflects this. Ground testing of the RL10 is done with altitude
simulation — either in a true altitude chamber, such as the NASA Plum Brook B-2
facility that fired more than a hundred RL10s during Centaur development, or on
a stand at sea level where a diffuser and steam-ejector system holds the
pressure at the nozzle exit down to a few psia. And when a sea-level-capable
RL10 was needed for the DC-X vehicle, the RL10A-5 variant was built with the
expansion ratio cut from 57–61 to 4. The tabulated sea-level figure for the
$\varepsilon = 61$ engine therefore does not correspond to a condition the
engine is operated in.

Note that the expansion ratio alone does not decide this; the ratio $p_e/p_{amb}$
does. The Space Shuttle Main Engine has a comparable $\varepsilon \approx 69$,
but its chamber pressure of roughly 200 bar puts the exit pressure near 20 kPa
rather than 3.85 kPa. It is overexpanded by a factor of about five at sea level
instead of 26, and is routinely fired at sea level.

There is a numerical problem alongside the physical one. At $\varepsilon = 61$
the pressure term is $1.860$ against a vacuum thrust coefficient of $1.937$, so
`CF_SL` is the small difference of two nearly equal numbers. A 1% error in
$C_{F,vac}$ moves `Isp_SL` by about 25%. Even as an idealisation, the value
carries almost no significant figures.

Pyskyfire reports `Isp_amb` and `Isp_SL` unconditionally and does not check for
separation. Treat them as valid only when the nozzle is close to matched — as a
rule of thumb, when $p_e$ is within a factor of two or three of $p_{amb}$. For
a high-expansion vacuum nozzle, use `Isp_vac`.

## From ideal to delivered performance

All four values are ideal: they come from a one-dimensional equilibrium
expansion and assume complete combustion. Real engines lose performance to
two-dimensional divergence, the wall boundary layer, finite-rate chemistry, and
imperfect injection.

## If you are interested in reading more

* [Binder, M., Tomsik, T., and Veres, J. P., "RL10A-3-3A Rocket Engine Modeling
  Project", NASA TM-107318,
  1997](https://ntrs.nasa.gov/api/citations/19970010379/downloads/19970010379.pdf)
  The source of the RL10A-3-3A design point and cooling-jacket data used in the
  `validation/RL10/` case. Appendix E.2 states the efficiency relation, and
  table E3 tabulates $\eta_{cf}$.

* [Binder, M., "An RL10A-3-3A Rocket Engine Model Using the Rocket Engine
  Transient Simulator (ROCETS) Software", NASA CR-190786,
  1993](https://ntrs.nasa.gov/api/citations/19950017370/downloads/19950017370.pdf)
  The earlier model of the same engine by the same author, and the predecessor
  to the report above.

* [RL10 variants and expansion ratios
  (Wikipedia)](https://en.wikipedia.org/wiki/RL10)
  Convenient comparison of expansion ratios across the RL10 family, from the
  4:1 sea-level RL10A-5 built for the DC-X to the 280:1 RL10B-2.

* [The RL10, a Usenet discussion by Bruce Dunn, Gary Hudson and Henry
  Spencer](https://yarchive.net/space/rocket/rl10.html)
  Informal but well-informed discussion of why the RL10 cannot be operated at
  sea level, including an estimate of the performance cost of shortening the
  nozzle enough to avoid separation.

* [NASA Plum Brook Station B-2 Spacecraft Propulsion
  Facility](https://azoony.com/hattrips/b2.html)
  Description of the altitude chamber in which much RL10 and Centaur testing
  was carried out.

* [RL10 rocket engine literature index
  (science.gov)](https://www.science.gov/topicpages/r/rl10+rocket+engine.html)
  Aggregated abstracts, including descriptions of the steam-ejector and
  diffuser altitude-simulation systems used for static firing at sea level.
