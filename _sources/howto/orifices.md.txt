# Size and solve orifices

`pyskyfire.ducts` has two restrictions described by a flow area and a
discharge coefficient:

- `Orifice`: a hole discharging into a large volume, such as an injector
  orifice into the chamber or a bleed into a cavity. The jet's kinetic energy is
  not recovered, so the whole drop is lost.
- `OrificePlate`: a thin round plate in a straight pipe, following ISO 5167.
  Part of the drop recovers behind it.

Both take the fluid as anything `thermo.as_equation_of_state` reads: a
CoolProp name, a `thermo.Fluid`, or an `EquationOfState` such as
`Incompressible(rho)`.

## Define an orifice

```python
from pyskyfire.ducts import Orifice

hole = Orifice.circular(1.2e-3, cd=0.75, length=6e-3)
gap = Orifice.annular(3.0e-3, 3.6e-3, cd=0.8)       # around a coaxial post
slot = Orifice.slot(4e-3, 0.5e-3, cd=0.7)
quoted = Orifice(area=1.48e-3, cd=1.0)             # an effective area given as Cd·A
```

The discharge coefficient is always supplied. Length is geometry only; it
does not change `cd`.

## Pressure drop at a flow, or flow at a pressure drop

```python
from pyskyfire.thermo import FluidState, PHASE_UNKNOWN

nan = float("nan")
plenum = FluidState(T=95.0, p=120e5, h=nan, quality=nan, rho=nan,
                    phase=PHASE_UNKNOWN, T_sat=nan)

result = hole.solve(inlet_state=plenum, mass_flow=0.05, fluid="Oxygen")
result.pressure_loss      # Pa
result.velocity           # ideal jet velocity at the vena contracta, m/s

result = hole.mass_flow(inlet_state=plenum, p_out=100e5, fluid="Oxygen")
result.mass_flow, result.choked
```

The flow expands isentropically from the inlet to the throat, so:

- a gas chokes where its mass flux peaks;
- a liquid taken below its vapour pressure flashes (homogeneous equilibrium)
  and stops responding to the discharge pressure, as a cavitating orifice does;
- `Incompressible(rho)` gives exactly `Cd·A·sqrt(2·rho·Δp)`.

Asking `solve` for more flow than the orifice can pass raises
`OrificeCapacityError`, which carries the largest flux found.

`solve_inlet_pressure` answers the remaining question: the upstream pressure
that drives a flow into a known discharge pressure, at a given inlet
temperature. This is what a feed plenum must hold.

```python
result = hole.solve_inlet_pressure(mass_flow=0.05, p_out=100e5, T_in=95.0, fluid="Oxygen")
result.inlet_state.p
```

```{note}
An effective area quoted from test data is tied to the equation it was reduced
with. If that was the incompressible equation at upstream density, using it
here for a gas overstates the drop: by about 16 % for hydrogen at Δp/p ≈ 0.12.
```

## Size a hole for a pressure drop

```python
hole = Orifice.sized(mass_flow=0.05, inlet_state=plenum, dp=20e5,
                     fluid="Oxygen", cd=0.75, length=6e-3)
```

`required_area` returns the area alone, for other shapes. Both refuse a drop
large enough to choke the flow, since such a hole would pass the flow at a
smaller drop.

## Orifice plates

```python
from pyskyfire.ducts import OrificePlate

line = FluidState(T=300.0, p=10e5, h=nan, quality=nan, rho=nan,
                  phase=PHASE_UNKNOWN, T_sat=nan)

plate = OrificePlate(pipe_diameter=0.1, orifice_diameter=0.05, taps="flange")
result = plate.solve(inlet_state=line, mass_flow=20.0, fluid="Water")
result.differential_pressure   # tap-to-tap
result.pressure_loss           # permanent, after recovery
result.discharge_coefficient
```

By default `cd` is the Reader-Harris/Gallagher correlation, which needs a
viscosity. It is looked up from CoolProp, or passed as `viscosity=`
otherwise. Outside the geometry and Reynolds-number range of ISO 5167-2 the
correlation is extrapolated, and `OrificeCorrelationWarning` says so.
Pass `cd=` to use a known value instead.

The plate uses ISO 5167's expansibility factor for gases rather than the
isentropic nozzle expansion `Orifice` uses. A thin plate's jet keeps
contracting past the bore, and the two differ by about 1.6 % in the drop at a
2.5 % tap drop in nitrogen.

## In an engine network

`OrificeBlock` replaces a fixed-ratio `SimpleDuctBlock` with the physical
restriction, such as a metering orifice in a film-cooling feed line:

```python
import pyskyfire as psf

blocks.append(psf.network.OrificeBlock(
    name="film_meter",
    st_in="film_line",
    st_out="film_ring",
    orifice=psf.ducts.Orifice.circular(2.0e-3, cd=0.62),
    fluid="Methane",
))
```

Injector holes have their own block, `InjectorBlock`; see
{doc}`Design and solve an injector <injectors>`.

It emits `dp_<name>`, so pump `overcome` lists pick the loss up unchanged. For
a plate the signal is the permanent loss. The drop now scales with the flow
and density rather than staying a fixed fraction of the inlet pressure.
