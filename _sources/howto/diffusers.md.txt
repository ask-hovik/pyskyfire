# Design a vaned diffuser

A vaned diffuser sits between an impeller and its volute. Its vanes take the
swirl out of the impeller discharge and turn the flow outward through widening
channels, so the volute behind it collects slow, low-swirl flow instead of flow
at the impeller's own speed.

Design one from the impeller it follows:

```python
from pyskyfire.pump import Impeller, VanedDiffuser

impeller = Impeller(Q=0.03, H=400., n=6000., shaft_diameter=0.03,
                    inlet_hub_diameter=0.04, blade_count=6,
                    blade_thickness=0.003, hydraulic_efficiency=0.95, rho=1000.)
diffuser = VanedDiffuser.from_impeller(impeller, vane_count=8)
print(diffuser)
print(diffuser.checks)
```

The inlet swirl is the impeller's slip prediction on its mean streamline, and
the head its outlet head. `rho` is needed for the gap rule and defaults to the
impeller's.

## How it is sized

The design follows Gülich, *Centrifugal Pumps*, chapter 7.9:

- **Gap to the impeller**, `gap_ratio = d3/d2`. It defaults to the Table 10.2
  minimum for low pressure pulsations,
  `1.015 + 0.08*(rho*H/(1000*1000) - 0.1)**0.8` below `nq = 40`. The gap grows
  with the energy density `rho*H`, so for a high-head pump it is ten percent of
  `d2` or more.
- **Inlet width**, `inlet_width_ratio = b3/b2`, 1.05 to 1.3. Default 1.15.
- **Approach flow** at `d3` (T7.7.3 to T7.7.5). The swirl falls across the gap
  as a free vortex. The vane inlet angle follows from the throat (see below),
  and its incidence on this flow is reported.
- **Throat**, from a chosen deceleration `deceleration = c3q/c2` (T7.7.11), 0.7
  to 0.85, default 0.8; or from `throat_factor`, Gülich's `f_a3` times the
  angular-momentum throat (T7.7.10), normally 1.1 to 1.3.
- **Channel length**, `length_ratio = L/a3`, 2.5 to 6, default 4. That sets the
  outer diameter `d4`. Give `outer_diameter_ratio` instead to fix `d4` and let
  the length follow. Gülich's own `d4` rule assumes a gap of a few percent, and
  the gap a high-head pump needs eats most of it, so the channel is sized
  directly.
- **Area ratio.** Each channel opens to the area ratio that recovers the most
  pressure for its length, `1.05 + 0.184*L/R_eq` (T1.7.13), so its recovery is
  Gülich's optimum, `0.36*(L/R_eq)**0.26` (T1.7.14), corrected for the
  channel's curvature (Eq. 1.46). His closed forms only cover the optimum, so
  the area ratio is not an input.
- **Flare**, `outlet_width_ratio = b4/b3`. The walls open linearly with radius.
  The flare opens the channel too, and counts toward the same opening-angle
  limit. Default 1, a plane diffuser.

## The vanes

The vanes are thin: `vane_thickness_ratio` of `d2`, 0.012 by default, within
Gülich's 0.01 to 0.015 (7.9.2). They are held at that thickness and rounded to
a point at both ends. Each vane is a camber line in two parts:

- **The leading spiral.** From the leading edge on `d3`, the vane runs as a log
  spiral, the shape of a free-vortex streamline. Its angle is the one that
  makes the throat to the next vane's leading edge exactly `a3`. The
  difference from the approach flow angle is the incidence, reported on
  `checks.incidence_deg`. It follows from the throat instead of being chosen,
  and lands within Gülich's ±3° for sensible throats.
- **The channel.** From the throat the vanes turn outward just enough that the
  channel between two of them widens at a steady rate, from `a3` to the exit
  width the optimum area ratio asks for. The gap between two vanes is about the
  pitch times the radius times the sine of their angle, less a vane thickness.
  That estimate is refined until the exit width measured on the drawn vanes is
  the design exit (`checks.drawn_exit_ratio`).

The channels therefore curve. `checks.deflection_deg` reports how far they turn
the flow. Gülich applies straight-diffuser data to about 25° of deflection. The
pressure recovery carries his curved-diffuser correction (Eq. 1.46)
regardless, which is conservative below that.

Straight channels would leave the throat much steeper than the flow reaching
it. The vane between two of them would then be a wedge filling most of the
ring, which is why they are not used here.

`build_geometry()` resolves one vane as a `BladeSurface` and the two walls,
symmetric about the impeller's outlet mid-plane. Both vane faces share their
camber stations, so the end caps join matching points. `viz.DiffuserComponent`
draws them in the impeller's own frame, so it takes the impeller's rotation and
a translation to that mid-plane.

## Checking a design

`diffuser.checks` holds Gülich's criteria for a finished diffuser (7.9.1 item
10). They are reported, never enforced; only a geometry that cannot be drawn
raises.

| Check | Gülich's range |
|---|---|
| `gap_ratio` against `minimum_gap_ratio` | at least the minimum (Table 10.2) |
| `inlet_width_ratio` | 1.05 to 1.3 |
| `throat_factor` | 1.1 to 1.3; below 1 at low `nq` is sometimes needed |
| `deceleration`, `c3q/c2` | 0.7 to 0.85, not much below 0.8 (Table 10.2) |
| `length_ratio`, `L/a3` | 2.5 to 6 |
| `aspect_ratio`, `b3/a3` | 0.8 to 2 |
| `exit_energy`, `c4²/(2gH)` | 0.02 to 0.04 for multistage pumps |
| `opening_angle_deg` | below `allowed_opening_angle_deg`, `16.5°*sqrt(R_eq/L)` |
| `vane_count_conflicts` | empty: no `|nu3*z_Le - nu2*z_La|` of 0 or 1 for `nu` up to 3 |
| `incidence_deg` | within ±3° |
| `deflection_deg` | up to about 25° for straight-diffuser data |
| `drawn_throat_ratio`, `drawn_exit_ratio` | 1: the drawn vanes meet the design |

## Losses

`diffuser.solve(inlet_state=..., mass_flow=..., transport=...)` charges
Gülich's Table 3.8 (2), each term referred to `u2²/2`:

- friction across the vaneless gap (T3.8.23);
- friction in the triangular inlet region (T3.8.18);
- the shock of decelerating `c2` onto the throat, the 0.3 term of T3.8.19;
- the diffusing channel, `1 - cp - 1/AR²` times the throat head;
- the channel jets mixing out over the outlet ring. They keep their swirl but
  lose their excess radial velocity (Borda-Carnot on that component).

The energy left in the mixed flow is not charged: the volute collects it and
charges its own loss. `diffuser.outlet` is what the volute should be designed
for: the outlet radius and width, the mean radial velocity over the ring, and
the mass-averaged swirl.

In an engine network, `network.DiffuserBlock` goes between the pump block,
whose outlet is the impeller outlet, and the volute block, which then takes the
diffuser's outlet station.

## Does it pay?

Not always. The channels recover only the optimum `cp` for their length,
typically 0.55 to 0.65, of a throat head that is still most of the impeller's
dynamic head. A vaneless volute behind a small, fast impeller can lose less.
In the dual-expander project, Gülich's diffuser and an Aungier volute behind it
lose about twice what the vaneless volute alone does. A vaned diffuser still
evens out the radial thrust, which a single volute does not. Compare the two
with the same loss models before choosing.
