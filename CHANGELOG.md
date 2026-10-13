# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

---

## [0.4.0] - 12-10-2026

### Added
- Pipes (`pyskyfire.ducts`). A `Pipe` is a bore swept along a `Path` of
  `PathStation`s, built from `PathNode` control vertices with `Path.from_nodes`
  or `generate_path_stations`, and joined by `Straight`, `Continuous` and
  `RoundedBend` connections. `Pipe.solve` marches stagnation pressure through
  wall friction, bends and area changes, and `PipeBlock` puts a pipe in an
  engine network. A bend whose centreline radius does not exceed the bore
  radius is rejected, since its inner side would fold.
- `Path.from_curve(curve, t_span, max_turn=...)` samples any curve-like object
  into a path: a SciPy spline such as `CubicSpline` or `make_interp_spline`, or
  any function mapping a parameter to points, vectorized or one point at a
  time, in 3D or in the x-y plane. Intervals are halved until none turns by
  more than `max_turn` (5 degrees by default), checking the turn inside each
  interval and across each station, so stations crowd into tight bends and a
  true corner raises instead of being cut. Each interval becomes the circular
  arc through its stations and the curve's midpoint: a `Continuous` carrying
  that arc's radius, length, centre and normal, or a `Straight`.
- `Pipe` places its own stations where the bore varies. It splits `Straight`
  and `Continuous` connections until the diameter changes by at most
  `max_diameter_change` (default 2 %) across each piece, placing the extra
  stations where the bore changes fastest and measuring inside each
  connection, so a bulge between equal end diameters is resolved too. The
  march evaluates a connection at its mean diameter, which without this
  underestimated the loss of a strong taper: 38 % on one 25-to-47 mm
  connection. A constant bore is left as supplied, and so is every
  `RoundedBend`. `Pipe.path`, `stations` and results refer to the resolved
  stations, which include the supplied ones. A diameter given per station is
  read as linear in arc length between them. `max_diameter_change=math.inf`
  keeps the supplied stations.
- `Pipe` and `AnnularDuct` charge the loss of a changing flow area, on top of
  wall friction: a 2 mm step from 25 to 50 mm loses 4.8 kPa, near the 4.7 kPa
  of a sudden expansion. A stretch over which the area keeps growing is one
  diffuser, and one over which it keeps shrinking one contraction, each ended
  where the area changes the other way or at a constant-section stretch longer
  than four hydraulic diameters. A diffuser loses the expansion term of
  Rennels' conical-diffuser correlation, and a contraction the convergence
  term of his conical-contraction correlation (both via `fluids`), for the end
  sections read as circles of the same area, at the cone angle averaged with
  weight where the flow slows or speeds up. That is Rennels' value for a cone
  and his sudden-expansion or sharp-contraction loss for a step, does not
  change as the pipe is split more finely, and is not watered down by a nearly
  parallel stretch at either end. `area_change_loss_coefficient` holds it per
  connection or interval, and `cone_angle` gives the local full angle of the
  cone of the same area.
- `Pipe(core_diameter=...)` puts a concentric core inside the bore: a
  centrebody, a turbine's nose, or the inner line of coaxial tubing. It is
  given as `diameter` is, and may start from a point and end, tapering or
  bluntly, where a growing core contracts the flow and an ending one diffuses
  it. Friction acts on the hydraulic diameter with the exact annulus laminar
  solution, a `RoundedBend` round a core uses the circular-pipe correlation
  on `solve(bend_length_scale=...)` with a `PipeCorrelationWarning`, and
  stations are refined on the area and the hydraulic diameter, which stay
  finite as the core shrinks to nothing. `Pipe.area`, `hydraulic_diameter`,
  `r_mean`, `core_radius` and `has_core` describe the section.
- `Pipe.solve(inlet_swirl=...)` carries swirl about the centreline, with or
  without a core, as `AnnularDuct` does: as angular momentum, worn down by
  friction along the helix the flow follows. `PipeBlock(inlet_swirl=...)`
  passes a fixed inlet swirl into the pipe.
- Annular ducts (`ducts.AnnularDuct`): an annulus about the shaft that may
  turn between axial and radial, such as a turbine inlet or exhaust, or a
  multistage pump's crossover. It is given as a mean line and a height, or
  re-paired from two walls with `AnnularDuct.from_walls`. Either way its span
  lines are normal to the mean line, so the area they sweep is the flow area.
  It shares `Pipe`'s march of stagnation pressure and carries swirl as angular
  momentum. Friction acts along the helix a swirling flow follows, and laminar
  friction uses the exact annulus solution
  (`common.friction.annulus_laminar_factor`, and a new `laminar_factor` on
  `darcy_friction_factor`). Meridional turns are charged as curved-duct
  friction. `viz.AnnularDuctComponent` draws it.
- `Pipe.solve` and `AnnularDuct.solve` return a `DuctResult`. Its `velocity` is
  the magnitude, through-flow and swirl together, alongside `area`,
  `meridional_velocity`, `swirl_velocity` and `swirl_angle`; `stations` and
  `connections` are None for an `AnnularDuct`. `p_stagnation` is the inlet
  pressure less the irreversible losses, and `p_static` subtracts the local
  `density*velocity**2/2`. `T_stagnation` is the inlet temperature, held
  because the duct is adiabatic, and `T_static` subtracts
  `velocity**2/(2*cp)`, so a duct solve needs a transport that provides `cp`.
  These are the names a volute uses too.
- `viz.PipeComponent` renders a solved `Pipe` directly, and
  `Pipe.display_path()` reopens lumped `RoundedBend` connections and
  `Continuous` arcs for drawing so a bend does not appear as a chord cut
  across its corner; `segments` counts drawn sides per quarter turn. Passing a
  `DuctResult` attaches its fields through `viz.duct_fields()`: station
  quantities as point fields, connection quantities as cell fields, with names
  matching the thrust chamber's so one `color_by` spans both. Refinement is
  display-only and does not touch any solver quantity. A pipe with a core is
  drawn as its core inside a translucent bore (`opacity`), the core only where
  it runs.
- `viz.PlotPipePressure` plots static and stagnation pressure, and
  `viz.PlotPipePressureLoss` the cumulative stagnation pressure loss, of one or
  more solved pipes against arc length. Markers sit on the solved stations.
- `examples/pipe/pipe_paths.py` holds all the code of the pipes how-to guide:
  rounded, continuous and flaring bends, a pipe along a SciPy cubic spline
  between two flanges, a quarter bend placed station by station with its
  curvature estimated and given, their solution and plots, and 3D scenes.
- Orifices (`ducts.Orifice`): round, annular, slot, or any area, with a
  supplied discharge coefficient, discharging into a large volume. The flow
  expands isentropically through the fluid's equation of state, so a gas
  chokes and a liquid below its vapour pressure cavitates (homogeneous
  equilibrium); a constant density gives `Cd·A·sqrt(2·rho·Δp)`. `solve` gives
  the drop at a flow, `mass_flow` the flow at a discharge pressure,
  `solve_inlet_pressure` the inlet pressure that drives a flow into a known
  discharge pressure, and `Orifice.sized` / `required_area` the hole for a
  drop. A flow the inlet cannot drive raises `OrificeCapacityError`.
- Orifice plates (`ducts.OrificePlate`) by ISO 5167: Reader-Harris/Gallagher
  discharge coefficient by default (warning outside its range), ISO
  expansibility for gases, and the permanent loss after recovery alongside
  the tap-to-tap drop.
- `network.OrificeBlock` puts either in an engine network, emitting its
  irreversible loss as `dp_<name>`.
- Volutes and collectors in `pyskyfire.ducts`, for the shaped manifold on a pump
  discharge and the plain constant-area toroidal manifold sometimes used on a
  cooling-circuit inlet. `EllipticalSlottedSection` is a bore fed through a
  fixed-width slot, with each wall-to-bore fillet tangent to both, so the bore
  centre follows from the bore size alone and the inlet height stays fixed as
  the passage grows. `ConstantAngularMomentum` sizes that passage from Aungier's
  Eq. 11-5 and `ConstantArea` gives the unshaped collector. The area law is
  required rather than defaulted, since a shaped volute and a uniform collector
  are different pieces of hardware, sized by different rules and costing
  different amounts of pressure. `Volute.solve` runs Aungier's one-dimensional
  loss analysis (Eqs. 5-105 to 5-109): meridional and tangential velocity
  dump and skin friction, each reported separately, with the exit cone
  charged as below. The sizing parameter defaults to `SP = 1.0`, which zeroes
  the tangential term. `VoluteBlock` adapts the solver to an engine network.
  A shaped volute wants zero area at the tongue, which no section family can
  provide, so stations where the area law could not be satisfied are flagged
  through `geometry.clipped` rather than failing.
- The volute bore is elliptical, and its aspect ratio doubles as the way past
  the area floor at the tongue. `EllipticalSlottedSection.aspect_ratio` defaults
  to one, a circle; above or below one it gives a bore oblong across the slot or
  outward from it. `size` is the semi-axis across the slot whatever the aspect,
  so the slot still floors it at `inlet_width / 2`.
  A shaped volute wants a vanishing area at the tongue and no bore spanning the
  slot can shrink that far, but flattening one shrinks its area roughly as one
  over the aspect ratio without narrowing it below the slot. `Volute` works
  down three regimes as the tongue is approached -- solve the size, then, once
  that has bottomed out, solve the aspect ratio up to `max_aspect_ratio`
  (default four), and only then clip -- which roughly halves the clipped run and
  cuts the area excess over the law near the tongue by an order of magnitude.
  `VoluteGeometry` carries `aspect` and `flattened`; `clipped` means the
  flattest allowed section was still too big. Each fillet is placed by solving
  its tangency to the ellipse, which bisects rather than converging faster on
  purpose: the root slides onto the end of its bracket and meets it tangentially
  as the bore closes onto the slot, where Newton stalls and a fixed-point
  iteration on the same residual diverges. `SectionEllipseArc` joins `SectionArc`
  in `common.section_geometry`.
- `Volute(tilt=...)` inclines the inlet plane away from radial, for a pump
  whose discharge points downward or a manifold following a nozzle's slope. It
  is constant around the circumference and rotates the section plane rigidly,
  so it decides where the passage sits and leaves every loss scalar unchanged.
  `Volute.sweep()` places the sections in machine coordinates.
- `Volute` takes an explicit `rotation`, one of `"cw"`, `"clockwise"`,
  `"ccw"` or `"counterclockwise"` in any capitalisation, naming which way the
  passage winds as it collects: counterclockwise seen from `+z` with the axis
  pointing at the viewer. It has to match the sense of the machine feeding the
  volute once both are placed, because a volute wound against the incoming
  swirl collects the flow backwards. Like `tilt`, it mirrors where the passage
  sits and leaves every loss scalar unchanged. `Impeller` takes the same
  argument about its own `+x` shaft, where it mirrors the blades and leaves the
  velocity triangles alone, since those are magnitudes in the rotating frame.
  Tangential at a volute outlet reverses with the sense, so a clockwise
  volute's exit cone leaves the other way round; radial and the elevation do
  not, being fixed by the machine axis rather than by which way the flow runs.
- `Volute(exit_cone=...)` attaches an `ExitCone`, the duct carrying collected
  flow from the full-collection plane to a pipe flange of a stated diameter, at
  any angle between tangential and radial and any elevation out of the
  volute's own plane. `first_angle` turns within the volute's plane from
  tangential, fixed at the volute outlet since the duct has moved round the
  volute by the time it reaches its flange; `"radial"` solves for the turn that
  points the flange straight out from the machine axis wherever it ends up,
  which is closer than a right angle away from tangential by exactly how far
  the duct has travelled round, and a numeric angle beyond that bound is
  rejected as pointing back inward. `second_angle` lifts the outlet out of that
  plane on top of the first angle. `arc_length_parameter` sets the duct's
  length in outlet diameters, with a floor of one; too tight a bend for that
  length is rejected rather than built self-intersecting. The flange diameter
  sets station 7, so `Volute.solve()` and `VoluteBlock` default
  `exit_diameter` from an attached cone. `Volute.sweep_exit_cone()` and
  `VoluteComponent` continue the swept duct seamlessly from the volute's own
  last ring -- bit-identical, not merely tracing the same curve at different
  points along it -- through to a circle at the flange; a tilted volute's
  cutwater edge carries its own tilted section frame into the duct rather than
  snapping to the untilted default. The example gives each of its three
  volutes a different outlet: tangential, radial, and an intermediate angle
  with elevation.
- The exit cone is charged as the conical transition it is, through the
  Rennels diffuser or contraction correlation according to which way it
  tapers. Aungier's Eq. 5-110 is the Borda-Carnot loss of a *sudden*
  expansion and carries no angle dependence; on a pump discharge diffusing 3:1
  through a near-optimal 8-degree cone it overcharges by roughly seven times,
  45 bar against the 6 to 8 bar that Rennels, Crane and Idelchik independently
  agree on, and for a flange narrower than the passage it is the wrong shape of
  formula altogether. It is kept only where there is no cone to resolve, the
  case it suits, and the two agree as the included angle opens toward 180
  degrees. `VoluteResult` reports `omega_exit_cone` and `omega_exit_friction`,
  the two halves of the one Rennels coefficient evaluated at the real friction
  factor and at zero, and `omega_exit_bend`, the loss of turning the duct from
  the resolved geometry. The three are disjoint and summed into `omega_total`,
  so `arc_length_parameter` is a real design variable, with an optimum trading
  diffusion and bend loss against wetted wall.
- Volute pressures: `stagnation_pressure()`, `static_pressure()` and
  `exit_cone_pressures()` give pressure along the passage and out the duct, and
  `volute_fields()` exposes them as `p_stagnation` and `p_static`. The lumped
  analysis says nothing about where round the volute each loss is spent, so the
  passage spreads them linearly with azimuth: an assumption made so it can be
  plotted at all, with the endpoints being what the model fixes. The duct's
  own ramp is between two stations the model resolves.
- `VoluteComponent` draws a resolved volute passage in a `Scene`, swept from the
  tongue through one revolution, so the step where the passage meets itself is
  the cutwater. `volute_fields()` colours it by any per-station
  quantity. Capping an end needs the section to be star-shaped about a known point,
  which `VoluteSection.fan_origin` supplies: for a slotted bore the area
  centroid will not do, since a long throat drags it inboard of the fillets.
  It returns `passage`, `inlet` and `outlet` child nodes, so the openings can
  be told from the walls. The slot mouth and the flange face are painted
  `opening_color`, red by default; they carry no fields, which is what keeps
  them that colour whichever field the scene is coloured by.
- `examples/volute/volute_report.py` designs two radial pump-discharge volutes
  and one manifold inclined 45 degrees, gives each its own scene, and reports
  their geometry and loss breakdown. Separate scenes rather than one shared one
  so each passage keeps its own colour scale: they differ in through-flow
  velocity by about a factor of three, and a shared scale flattened the two
  slower ones to a single tone.
- `Manifold`, the distributing or collecting ring, alongside `ManifoldLaw`
  with `UniformArea` and `ConstantVelocity`. It is not a volute with a
  different law: nothing feeds it with swirl, so it has no angular momentum to
  conserve and its circumferential velocity is a cost rather than a design
  variable. What replaces a volute's direction is a topology. `branches=1` runs
  one way round the whole ring and so needs a `rotation`; `branches=2` splits at
  a dead point, runs both ways and meets at the port, which is mirror-symmetric
  and refuses a rotation because there is no sense to name. Only one branch is
  resolved, and `sweep` repeats it to close the ring. Multiple ports arrive as
  more branches rather than as a second class.
- `ConstantVelocity` holds the through-flow velocity, and with it the dynamic
  head, uniform along a branch, so every channel tapping the ring sees the same
  static pressure. `ConstantVelocity.for_duty` turns a duty into the area the
  law is stated in, so the law itself carries no fluid state. Holding *static*
  pressure uniform is a different rule again and needs the branch pressure drop
  to couple against, so it waits for that coupling.
- `HalfMoonSection` is the manifold cross-section: the chamber wall itself on
  the inboard side, a half-ellipse outboard, and a circular fillet of a stated
  radius blending the two at each end. It reads as a D on a cylindrical wall
  and as a kidney where the wall curves. The slit and the width are arc lengths
  along the wall rather than axial extents, because the slit is room reserved
  on the wall surface for channel breakout and a sloped wall runs further than
  its station span does. The slit is a reservation, not modelled geometry:
  pyskyfire draws the cavity the channels open into and leaves the breakout
  itself to CAD. An anchor fixes the slit and it stays put at every azimuth,
  with `slit_placement` setting where it sits about the anchor and
  `exposure_placement` dividing the spare wall either side of it in a ratio
  held all the way round. Sizing follows the volute's three regimes: the
  exposed wall closes first, then the ellipse flattens toward
  `max_aspect_ratio`, then the area law is clipped.
- `VoluteSection.maximum_size` reports the largest section a family can build,
  unbounded by default. A manifold section has only so much chamber wall to
  spread onto, so a law asking for more area than exists is told exactly
  that, instead of the size search being chased past the end of the wall.
- `ManifoldComponent` renders a ring: its walls, the slit reserved for channel
  breakout, and its stub. Sections say which of their boundary pieces is the
  opening through `opening_index`, so a slotted section and a half-moon one can
  both be swept apart into walls and opening. Only the faces flow actually
  crosses are painted as openings: the slit and the flange the stub ends on.
  The rest -- the dome, the stub's own tube, and the wall outside the slit --
  is wall, takes the default colour and carries the fields.
- `Manifold.sweep_slit` and `HalfMoonSection.slit_boundary` draw the reserved
  slit apart from the chamber-facing face it sits on, so a layout can be
  checked against where the cooling passes actually start and end, and
  `Manifold.sweep_margins` and `HalfMoonSection.margin_boundaries` return the
  chamber-facing wall either side of it, tiling that face with the slit
  between them.
- Manifolds can be drawn on a finer grid than their area law was resolved on.
  A swept ring facets in azimuth long before it facets around its section -- it
  is metres round and centimetres across -- so `Manifold.sweep`, `sweep_slit`,
  `sweep_margins` and `fan_origins` take an `n_stations`, and
  `ManifoldComponent` takes one too. It is drawing resolution and not a finer
  design: the sizes the law resolved are reproduced exactly at their own
  stations and read between them, everything following from a size is computed
  at the new angles, and field values stay on the grid they were solved on,
  interpolated onto the refined mesh the way a pipe's are along an opened bend.
  `Manifold.sweep_stations` says how many stations a sweep will return.
- A `Port` protocol for the duct carrying collected flow to its flange, with
  `ExitCone` and the new `RadialPort` as its two implementations and
  `ResolvedPort` as the placed duct they both produce. A volute's duct leaves
  tangentially and turns; a chamber manifold has no such direction, so its duct
  leaves along the way the section faces and runs straight. Splitting the
  request from the placement means either passage can carry either duct -- a
  manifold fed tangentially is a real arrangement, and it does not cost a
  second implementation of the same geometry. `PortSite` is what the passage
  knows about the attachment, handed over at resolve time; each port asks for
  the fields it needs and names the missing one when it cannot proceed.
- `RadialPort` is the chamber-manifold stub: a straight tube leaving the dome
  its section is built on. It is aimed from the ellipse centre, where that
  section's two axes cross, and steered by a single `angle` in the section
  plane. The angle is read either from the way the section faces
  (`angle_reference="normal"`, so zero points straight out of the dome however
  the wall leans there) or from the plane the manifold rings
  (`"manifold_plane"`, so the same number is the same hardware direction at
  every station); the two differ by the section's own lean. `length` is
  measured from that same centre rather than from the surface, so the flange
  does not move when the angle or the area law does, and a stub too short to
  put its whole flange clear of the dome is refused. The stub takes the duct's
  own frame, and its first ring is trimmed generator by generator back onto
  the curved dome, so the tube sits on the passage rather than touching it at
  one point.

  How wide the stub may be follows from where it is aimed. It is not the
  ellipse's own semi-axis: past the fillet tangency the ellipse runs outboard of
  the metal, so what bounds the stub is how far the part of it that is really
  passage boundary reaches across the direction chosen. That closes down as the
  stub turns toward the fillets, and a stub aimed past them is refused outright
  rather than allowed to leave through the side. `diameter_fraction` is taken
  against that room, defaulting to 0.95 of it. The stub stays centred on
  the ray it was aimed along, so a tilted one is sized by whichever side runs
  out of dome first.
- `ThrustChamber.outer_wall` reports the cold-side envelope of the channelled
  wall, the counterpart to `contour`: where the contour is the face the hot gas
  sees, this is the outermost surface the channels and their walls reach, which
  is what a manifold wrapped around the chamber sits against. It is read the
  same way -- `r(x)`, `dr_dx(x)`, `normal(x)` -- and carries two radii, because
  a channelled wall is not a cylinder. `envelope` is the crest of the bumps and
  is what anything needing to clear the chamber is sized against; `equivalent`
  is the plain cylinder enclosing the same area, and is what belongs in area,
  volume and mass bookkeeping. They coincide for a squared channel, whose roof
  really is a cylinder, and separate only where the exterior is scalloped.
  Circuits interleaved around one circumference, such as an RL10-style half
  pass and return, contribute to the envelope by whichever reaches furthest and
  to the equivalent by the area each owns over its own share of the angle.
- `OuterWall.profile(x)` takes the wall off a resolved chamber as an
  `InnerProfile`, measured by arc length from the anchor, ready to hand to a
  section.
- `CoolingCircuit(outer_walls=[...])` adds cold-side liners outboard of the
  channel section, the way the silver insert sits on the hot side. They raise
  the outer wall by their own thickness and never enter the channel-wall
  thickness the section is built from, so a squared channel's roof and a
  rounded tube's own material are left exactly as they were.
- General path-based heat exchangers (`pyskyfire.heat_exchanger`) with
  prescribed chamber gas, coupled counterflow/coflow fluids, equal-flow copied
  channels, explicit shared-wall segment links, ideal series headers, and
  wall-temperature reports. Hydrogen/oxygen and paired ring/spoke baffle
  examples. Heated channels correct Gnielinski for the wall temperature: a
  heated fluid's turbulent Nusselt number is multiplied by `(T_b/T_w)^0.45`,
  Gnielinski's factor for gases
  (`common.physics.gnielinski_temperature_correction`, and
  `bulk_temperature`/`wall_temperature` on `h_internal`). Bulk properties
  alone miss how the fluid thins against a hot wall: for cryogenic oxygen
  against a baffle wall near 1000 K the factor is about 0.4-0.6. The wall
  temperature is the wetted interface's, lagging one coupling iteration, so the
  converged coefficients use the converged walls. Cooled fluids and the laminar
  branch are left uncorrected, as is a circuit's own `heat_transfer` callback.
  `HeatExchangerCircuit(property_correction=False)` keeps bulk properties.
- `HeatExchangerComponent` brings transported channel geometry and wall/flow
  fields into the shared scene viewer, preserving chamber installation
  transforms. It draws each channel's wetted passage, coloured by
  `T_coolant_wall`; `result.walls` carries both interface temperatures. The
  baffle example includes a translucent chamber cooling layout and an offline
  combined scene embedded in its report.
- Shared `ducts.Path`, `PathStation`, `PathPlacement`, and transported section
  frames.
- Gnielinski convection in `common.physics` and shared acceleration pressure
  change in `common.momentum`.
- Duties: what a turbomachine is asked to do, separately from how it is built.
  `pump.PumpDuty` takes flow, head and speed as pumps are specified
  (`from_head`), a mass flow and pressure rise at one density
  (`from_pressure_rise`), a medium and two pressures (`from_states`), or two
  network stations (`from_stations`), and stores the same thing each way: the
  isentropic head, the inlet flow and the flow at the isentropic outlet.
  `turbine.TurbineDuty` does the same for an expansion and adds the spouting
  velocity, pressure ratio, `u/c0`, and Balje's specific speed and diameter.
  A turbine duty always carries its states, because the gas expands on the way
  through and blade heights are sized on that. `Impeller.from_duty` sizes an
  impeller for a pump duty. The shared fields live on `thermo.MachineDuty`.
- Equations of state: `thermo.Incompressible(rho)` and `thermo.CoolPropEOS`
  answer the questions a machine design asks of a fluid: the state at
  `(p, T)`, at `(p, h)`, and after an isentropic process to another pressure,
  plus `entropy`, `state_hs` and `speed_of_sound`. `CoolPropEOS` resolves
  through `FluidStateResolver`, so the saturation dome is handled as it is in
  the cooling march. A `(p, T)` pair on the saturation line reads as saturated
  liquid, as at a pump inlet. `CoolPropEOS.state_hs` searches along the
  isenthalp where CoolProp's own `(h, s)` flash refuses. That happens past a
  fluid's temperature range: methane's ends at 625 K, below a methane
  expander's turbine inlet.
- Centrifugal impeller design (`pyskyfire.pump`). `Impeller` sizes a radial
  impeller from flow, head and speed at a supplied hydraulic efficiency, which
  is the impeller's own, Gülich's eta_h,La; `H` is the head at the impeller
  outlet. The design is closed on its head: exactly one of the blade outlet
  angle and the outlet diameter is solved so the predicted swirl is the swirl
  the head needs. By default the angle is solved and the diameter comes from
  `psi`, `d2`, `u2`, a `meridional`, or Eq. (3.26) at `fT` (default 1.0);
  given `beta2_deg`, the diameter is solved. `b2` can be given, and
  `blade_count` must be at least 3, the range of the slip correlation.
  `Impeller.blade_power` is the work the blades do; mechanical and disk
  friction losses are charged where they happen, not lumped into an overall
  efficiency. `Impeller.diagnostics` (`pump.ImpellerDiagnostics`) reports what
  the head was closed on, the resolved and correlated head coefficients, the
  slip prediction, the deviation angle and the de Haller ratio `w2/w1`, against
  Gülich's guidance and never enforced. `build_geometry()` draws the pressure
  and suction surfaces with pyTurbo Aero (`pip install 'pyskyfire[impeller]'`),
  with optional splitters (`SplitterSettings`), and `viz.ImpellerComponent`
  places them in a scene. `examples/impeller/impeller_report.py` designs a
  low-specific-speed impeller with splitters and reports it.
- Generated impeller blades are stacked on the trailing edge, so it lies in one
  axial plane. Each span section wraps by the integral of its own metal-angle
  law over its own conformal length, and the hub section is much longer in
  `m'` than the tip section, so the two wrap differently by construction: 216
  and 121 degrees in the impeller example. pyTurbo builds every section from a
  common leading edge, which would push that 95 degree difference onto the
  trailing edge, more than the 60 degree blade pitch. The leading edge carries
  the difference as sweep instead, the way a blade development is stacked.
  Stacking is a rigid rotation of each section about the shaft, so metal
  angles, wrap, thickness and the passage are unchanged, and a splitter is
  rotated with the main blade it follows.
- `BladeSettings(start=(hub, tip))` slants an impeller's leading edge, starting
  the blade further downstream on the hub than at the tip. A hub streamline is
  far longer in `m'` than a tip streamline, since `dm' = dm/r` and the hub runs
  close to the axis, so a leading edge placed on a single station makes the hub
  section wrap much further than the tip section -- 178° against 101° in the
  impeller example -- and the difference shows up as leading-edge sweep. A
  slanted leading edge is how a blade development is normally drawn, and it
  evens the wrap out: `start=(0.34, 0.05)` takes that spread from 77° to 14°.
  A section can only begin on a span line, so the slant is expressed as which
  hub and shroud stations pair, and `MeridionalContours.with_leading_edge`
  re-pairs the passage to carry it. The inlet velocity triangle and the inlet
  metal angles come from the leading edge rather than the inlet plane, with
  `cm` from unblocked continuity through the component of the through-flow that
  crosses it, taking the flow direction from the two wall tangents so that the
  area a leading edge passes does not depend on how the stations were paired to
  place it; eye sizing still uses the inlet plane and `phi1`. The example uses
  a slanted leading edge.
- Impeller slip, from Gülich's Table 3.2 (`pump.outlet_slip`): a slip factor
  from the blade outlet angle, blade count, blade blockage and inlet diameter.
  The outlet triangles carry the swirl the blades are predicted to leave,
  and `beta_flow_deg` the flow angle it implies.
- Disk friction, from Gülich's Table 3.6 (`pump.DiskFriction`,
  `pump.disk_friction_power`): the shear on a wheel's outer faces, with the
  face roughness, and no leakage through the side gaps yet.
- Vaned radial diffuser (`pump.VanedDiffuser`), designed from the impeller it
  follows by Gülich's chapter 7.9. The gap to the impeller is his
  pressure-pulsation minimum (Table 10.2), the vane angle the approach flow,
  the throat set by a deceleration or by `f_a3`, and each channel opened to
  the optimum area ratio for its length, so its recovery is his closed-form
  optimum, corrected for curvature. The vanes are thin, of his thickness, a
  leading log spiral that draws the design throat and then turning so the
  channel widens steadily to the design exit; both are checked on the drawn
  vanes. Losses follow Table 3.8 (2) plus the channel jets
  mixing out over the outlet ring. His quality criteria are reported on
  `checks`, never enforced. `outlet` is what the volute behind it collects.
  With `network.DiffuserBlock` and `viz.DiffuserComponent`.
- `common.MeridionalPassage`: the wetted flow path of a turbomachine or annular
  duct, as paired hub and shroud `(x, r)` curves. It reports area and
  hydraulic diameter per station, samples by span, and revolves its walls into
  surfaces. It also records which wall stretches rotate with the shaft
  (`WallMotion`) and where a blade edge stops short of a wall and by how much
  (`Clearance`). Efficiency evaluators need both and the curves alone don't
  say them. `pump.MeridionalContours` is its radial-impeller form, keeping
  only the checks that a radial wheel satisfies.
  `MeridionalPassage.end('inlet' | 'outlet')` returns the passage end as an
  `AnnulusEnd` (radii, area, wall directions), for whatever joins it there.
  `BladeSurface` and `cylindrical_to_cartesian` live in `pyskyfire.common`,
  shared by impellers and turbines.
- Axial turbine design (`pyskyfire.turbine`). `AxialTurbine` is a meanline
  design at constant mean radius from a `TurbineDuty` and a supplied
  efficiency. Each stage is given a flow coefficient and either a reaction or
  an exit swirl. Work is shared between stages by `work_split`. The entropy
  rise the efficiency implies is shared among the rows by `loss_split`,
  which sets each station's density and so its blade height. It reports
  the velocity triangles, static states and Mach numbers at every row edge,
  and draws the result as a `TurbineDefinition`. A `TurbineDefinition` is
  hardware only: `BladeRow`s (count, chords, annulus radii, metal angles,
  stagger, thicknesses, clearances, admission and a `twist` law, free vortex
  from the meanline) in a `MeridionalPassage`. The passage records the hub
  turning under rotors, the casing over shrouded ones, and every blade gap.
  `build_geometry()` draws each row's sections with pyTurbo
  (`pip install 'pyskyfire[turbine]'`) and wraps them
  onto the annulus. It also measures what the drawing leaves to the built
  blade, at every span: each row's throat (from the trailing edge) and gauging
  angle, the inlet width, the narrowest gap in the passage and where it is
  (with `pinch` relating it to the inlet width and throat), and the actual
  maximum thickness. `RowSettings` takes a `SectionProfile` for stator rows
  and one for rotor rows: thicknesses, trailing-edge wedges, a
  `ThicknessDistribution` along the camber line, and where the chord lies
  between the metal angles. `SectionProfile.impulse()` starts an impulse
  rotor. `viz.TurbineComponent` shows the passage and every blade.
  `viz.plot_velocity_triangles`, `plot_turbine_meridional`,
  `plot_turbine_cascade` and `plot_turbine_spanwise` plot the design.
  `examples/turbine/turbine_report.py` designs a single-stage impulse oxygen
  turbine and reports it with these plots and the 3D scene.
- A `turbine` install extra, pinning the same pyTurbo Aero as `impeller`.
- Injectors (`pyskyfire.injector`). An `Injector` is a list of `ElementSet`s,
  each one element type at supplied `(r, θ)` face positions and orientations.
  `ImpingingElement` is the first element type: `Jet`s, each with its own
  feed, `ducts.Orifice`, inclination and azimuth, aimed at one impingement
  point, which covers like and unlike doublets, triplets and pentads. A feed
  is whatever supplies a hole (feed hardware such as a `ducts.Manifold`, or a
  name), so a film ring fed separately from the main plenum is its own feed.
  Each feed is solved as one equivalent orifice, since all its holes share a
  plenum and a face; `solve` finds each feed's flow from its plenum pressure
  or its plenum pressure from its flow, with per-set flows, the mixture
  ratio, and each element's resultant momentum angle, element mixture ratio,
  Rupe ratio and momentum-flux ratio, reported and not enforced. `checks`
  reports webs between holes, overlaps, the clearance from the holes to the
  chamber wall and, separately, from the impingement points to it (zero for a
  jet aimed at the wall on purpose), and hole L/d. A single jet aimed
  shallowly at the wall is an impinging element too, for film cooling. With
  `network.InjectorBlock` (one per feed), `viz.InjectorComponent` (the wetted
  bores, coloured by feed and ending at the face plate, with per-hole flow
  fields and optional arrows to the impingement points) and
  `viz.plot_injector_face`. Example in `examples/injector/`.
- Modular 3D scenes (`viz.Scene`) with thrust-chamber, impeller, and pipe
  components; SI rigid placement and visual inlet/outlet anchors; optional
  cooling-result scalar fields in desktop and offline HTML output. Cooling
  channel instances share geometry built once per resolved seed. Standalone
  pages carry each mesh and each field's values once, however many parts
  share them, so a 280-channel chamber's page is about 0.7 MB, and drawing it
  through every station (`stride=1`) about 1.7 MB. Parts are drawn in one
  default colour unless a component has a reason to set its own:
  `circuit_colors` on the chamber and heat exchanger and `colors` on the
  turbine paint them apart, inlet and outlet openings are red, and injector
  bores are coloured by feed. See the 3D scene guide and
  `examples/visualisation/scene.py`.
- Clicking a part in the standalone 3D viewer selects the top-level component
  that owns it and reports, in a box beside the click, its name, kind, anchors
  in scene coordinates, and the `position=`/`rotation=` arguments that reproduce
  its placement. The angles are derived in Python through `Transform`'s own
  extrinsic xyz convention, so they can be pasted back into `scene.add()`
  without reinterpretation.
- Groups and multiple selection in the standalone 3D viewer. Shift+click adds
  or removes components from the selection, and the info box then lists them.
  Selected parts are tinted a calm blue, as in CAD; the tint is mixed in, so a
  field's colours still show through. Right-click lists every part under the
  cursor, front to back: each component hit, its circuits and other named
  parts, and the channel or blade actually under the cursor. Hovering an entry
  tints the part lightly, so a circuit buried inside a chamber can be found
  and selected. The info box is dragged by its heading and stays where it is
  put. A selection is saved as a named group from the info box. The groups
  panel offers those, one group per fluid and one per tag, and clicking a
  group selects its members. The legend's range selector gives any group the
  colour scale to itself and fades the rest of the scene, so a pressure drop
  through one circuit is not lost in the span of the whole engine. Groups
  saved in the viewer are kept in the browser's local storage per scene
  title, where the browser allows it.
- A field dropdown in the standalone viewer's legend recolors the scene by any
  field its components carry, with a "None" entry restoring component colors.
  The HTML already carries every field, so this costs no payload and no rerun:
  `color_by` chooses the initial selection rather than the only one, and
  `clim` pins that field's range while the rest range over their own values.
  Components lacking the selected field fall back to their material color.
- `Scene(show_axes=True)` draws a corner orientation triad labelled x, y, z in
  both the desktop and standalone HTML viewers. It marks the scene frame's
  directions only: it sits in a fixed corner viewport, so it carries no ticks
  and shows neither the origin's screen position nor the scene's scale. The
  RL10 report's engine scene turns it on, so the duct's parked orientation can
  be read against the chamber axis.
- `Scene.add(..., tags=...)` and `SceneNode.tags` label parts with free-form
  tags. Every tag is a group in the viewer, so groupings written in the script
  come back each time the results are regenerated.
- Solved components record their fluid. `DuctResult`, `VoluteResult` and the
  injector's `FeedResult` carry `fluid`, the CoolProp name their transport or
  feed condition named, and `thermo.fluid_name` reads one from a string, a
  `thermo.Fluid`, or anything carrying one on `.fluid`. Components given such
  a result, a `RegenResult` or a `HeatExchangerResult` record the fluid as
  `metadata["fluid"]` on the part carrying it, under one name however CoolProp
  was asked for it (`H2`, `HEOS::Hydrogen[1.000]` and `Hydrogen` are one
  group). A mixture is named by mass fraction, as propellants are mixed,
  rather than CoolProp's mole fractions: `Ethanol[0.610]&Water[0.390]` reads
  `Ethanol 80 wt% & Water 20 wt%`. Impellers and turbines take their results
  as plain fields and record none.
- `validation/RL10/plumbing/`, where explicit RL10 tubing accumulates one
  module per physical component. `bending_table.py` converts a fabrication
  bending table into `PathNode` control vertices, absorbing each corner's
  `R tan(angle / 2)` tangent distance and rolling bend planes about the
  intervening straight. `turbine_discharge_duct.py` is the Figure G1 duct: it
  reproduces the tabulated tangent distances and the 51.7 in actual duct
  length, and takes the 2.402 in bore from the drawing's OD and wall callouts.
  The duct is drawn in the RL10 report's 3D engine scene as a `PipeComponent`,
  at true size but parked beside the chamber, since its mounting point on the
  turbine housing is not modelled. Nothing is wired into `sizer_sim.py` yet,
  because a `SimpleDuctBlock` pressure ratio covers a whole network leg and can
  only be retired once every component along that leg is modelled.
- `solids.GRCop42_AM` and `solids.GRCop84_AM`: laser-powder-bed-fusion GRCop
  from Chen et al. (2023), Materials Today Communications 36, 106665. GRCop-42
  conductivity is their Eq. 3, the mean of eight builds; GRCop-84, for which
  they give no equation, is a cubic fitted to the two builds in their Fig. 2(b).
  Both are measured over 298-973 K, held at the 298 K value below it and
  continued linearly at the 973 K slope up to the solidus. The instantaneous
  CTE is the derivative of their thermal-strain fits (Eqs. 9 and 10).
- `solids.SilicaDeposit`: the silica layer a silicon fuel additive such as
  TEOS lays down on the hot wall, for use as a `role="liner"` wall. Its
  conductivity is dense vitreous silica without radiative transfer, from
  Carwile and Hoge (1966), U.S. Army Natick Laboratories TR 67-7-PR (DTIC
  AD0636276): their Table 1 over 80-660 K and their estimated non-radiative
  curve, digitised, up to 1960 K. A porous deposit conducts less, so its
  thickness is the dense-silica thickness with the same `t/k`.
- Explicit wall roles: `Wall(..., role="main")` is the default;
  `role="liner"` identifies combustion-side layers. Channel thickness sums
  all main layers, while liners contribute to hot-side offset and conduction.
  The RL10 silver throat insert is explicitly marked as a liner.
- `regen.ChannelKeepOut` reserves equally spaced sectors of the channel layer
  over an axial span, e.g. for hardware passing through the wall. It is given
  to `ThrustChamber(keep_outs=...)` and is not a circuit: it is not solved,
  drawn or exported. Surface channels in its span share only the circumference
  it leaves free, so all of them narrow and the thermal solve sees it. The
  numerical layout pins each sector as an obstacle. It first opens room with a
  continuous map and then relaxes, so neighbouring channels bend aside without
  overlapping or crossing. A `taper` runs the width out beyond each end so the
  bend is gradual. Circuit grids get a station on every corner of the width
  law. The analytical layout raises `AnalyticalLayoutError` when keep-outs are
  present, since constant phase shifts cannot bend around them.
- `SurfacePlacement` takes `channel_angular_width`. Stating the pitch holds it
  where the channel count changes, so a short stretch where only one of two
  interleaved passes runs keeps its width instead of doubling to fill the room.
- `redistribute_channels` takes `fixed`, for entries that push but never move,
  and `shape_offsets`, a detour its bending and slope terms treat as straight.
- Shared lazy cooling-channel section planes plus analytical and numerical
  layout resolvers, including validation for exactly patternable interleaves.
- Design-input functions: geometry inputs declare the quantities they need
  by parameter name, drawn from a tiered vocabulary of contour, layout and
  channel-section quantities. See the new *Design-Input Functions* explanation
  page for the slot table and the full glossary.
- An advanced-example Onshape JSON exporter and matching FeatureScript that
  loft unique channel seeds and circularly patterns their instances. The
  FeatureScript builds the true constant-width-rib section, whose flat rib
  faces make the channel widen with radius while the rib does not; the JSON
  schema is at version 8 and carries a `sectionType` per channel and the two
  resolved construction arcs per station. The exporter rejects cross sections
  its FeatureScript cannot sketch instead of exporting arcs that would be
  mis-built.
- The Onshape export always writes the hot wall, with circuit hand-off
  transitions and Bézier blends between unlike neighbours, tuned by
  `--hot-wall-trim-fraction` and `--hot-wall-transition-fraction`; the chamber
  must carry a contour. The FeatureScript lofts the wall one exported curve at
  a time along each `run` of sections with constant profile topology, then
  knits the strips. Pairing whole profiles by proximity instead failed past
  roughly ninety curves per profile, and twisted the wall against the channels
  once a helix advanced more than half a channel pitch between stations. The
  two profiles either side of a hand-off share no curves, so the export cuts
  both onto one set of azimuths and they pair up as well. Plain hot-wall spans
  revolve one spline per smooth run of the chamber contour, so a revolved wall
  is a single smooth face rather than hundreds of short segments. Runs are cut
  where the meridian's tangent jumps, so a throat corner or an interpolated,
  faceted contour keeps its shape.

### Changed
- **Breaking:** engine-network classes move from `pyskyfire.common` to the new
  `pyskyfire.network`. `EngineNetwork`, `Station`, `FluidBlock`,
  `SignalBlock`, `PumpBlock`, `RegenBlock`, `TurbineBlock`, `SimpleDuctBlock`,
  `MassFlowSplitterBlock`, `MassFlowMergerBlock` and `TransmissionBlock` are
  now `psf.network.<name>`; replace `psf.common.EngineNetwork` with
  `psf.network.EngineNetwork` and so on. No backwards aliases are provided.
- **Breaking:** `pyskyfire.skycea` is now `pyskyfire.thermo`, which also owns
  `Fluid` (formerly in `common`) and fluid-state resolution (formerly in
  `regen.coolant_state`). `CoolantTransport` becomes `FluidTransport`,
  `CoolantState` becomes `FluidStateResolver`, `BulkState` becomes `FluidState`,
  `CoolantCapability` becomes `FluidCapability`, and `GasState` becomes
  `CombustionGasState`. `TransportProperties`, `SaturationState`,
  `Aerothermodynamics`, and `NozzleLoss` keep their names under `thermo`.
  `probe_coolant` becomes `thermo.state.probe_fluid`; film-property policy moves
  to `regen.coolant_properties.film_properties(resolver, T_film, p, getters)`.
  No backwards aliases are provided. Old stored pickle files may require
  regenerating due to module and class-name changes; compatibility is not guaranteed.
- Heat-transfer physics and friction now live in `common.physics` and
  `common.friction`; regen and ducts import the shared implementations, and
  `pyskyfire.regen` still re-exports the physics functions. The
  `pyskyfire.regen.physics` module is gone.
- **Breaking:** `pump.Impeller` is the new design described above. It no
  longer accepts the 0.3.0 `Impeller(Q, H, n, material)` call: it needs the
  shaft and inlet hub diameters, blade count and thickness and a hydraulic
  efficiency.
- **Breaking:** `network.PumpBlock` takes the impeller's `hydraulic_efficiency`
  in place of `eta`, plus an optional `mechanical_efficiency` and
  `disk_friction`. Its outlet is the impeller outlet; shaft power is blade
  power plus disk friction over the mechanical efficiency, and `post_process`
  reports that split. Losses downstream, such as a volute's, are charged by
  their own blocks rather than folded into an overall efficiency.
- **Breaking:** `enable_fin` now belongs to `CoolingCircuit` (default `True`),
  not `ThrustChamber`. Each circuit can also set `rib_material`; omitting it
  retains the coolant-facing main-wall material for the fin calculation.
  The RL10 case sets fin enhancement on its two circuits individually.
- **Breaking:** `CrossSectionSquared` now takes `rib_width` instead of
  `blockage_ratio`, and models the rib as having the same thickness from the
  channel floor to the channel roof rather than as a wedge. Pass
  `rib_width=lambda pitch_mid: f * pitch_mid` to reproduce a blockage ratio of
  `f` exactly. The 3D viewer draws the area-equivalent wedge, since it cannot
  yet draw the true section.
- **Breaking:** design-input functions are now called with keyword arguments
  and must name what they need. `lambda x: ...` keeps working, a function that
  ignores position becomes `lambda: ...`, and `lambda *args: ...` is rejected
  with an explanatory error. This covers `Wall(thickness=...)`,
  `CoolingCircuit(channel_height=..., roughness=...,
  hot_gas_surface_area_multiplier=...)`, and the placement classes'
  `helix_angle` and `channel_angular_width`.
- `SurfacePlacement.helix_angle` and `InternalPlacement.channel_angular_width`
  are now methods that evaluate the design input, rather than stored callables.
  They, and `dtheta_dx`, take the circuit as an optional trailing argument so
  that a placement shared by several circuits evaluates against each in turn.
- Renamed angular cooling-channel width APIs from `channel_width` to
  `channel_angular_width`; previously saved result files are migrated when loaded.
- **Breaking:** a cooling circuit now carries three explicitly named paths in
  place of `centerlines`/`coolant_centerlines`. `reference_path` is the ideal
  contour carrying the nominal winding, `cross_section_base_point_path` is
  where the cross section is anchored (the coolant-facing floor for a squared
  section, the inner lobe centre for a rounded one), and
  `cooling_channel_centroid_path` is the coolant-area centroid that sets
  passage length and curvature. Circuits with the same helix angle now share
  one winding and cannot drift apart. Saved geometry using the old definitions
  is rejected on load rather than renamed; regenerate it.
- `solids.GRCop42` is now the AM model above. Its conductivity was a
  sum-of-Gaussians fit that gave 421 W/m·K at room temperature, above pure
  copper, and a pure-copper-like peak near 10 000 W/m·K at 20 K. It now gives
  329 W/m·K at room temperature and 297 W/m·K at 973 K.
- The standalone 3D scene viewer rotates with `TrackballControls` instead of
  `OrbitControls`, so the view can be spun continuously to any orientation
  rather than stopping at the poles of a fixed up-axis, as in CAD. Pressing
  `R` restores the initial view, and Ctrl+drag (Cmd+drag on a Mac) pans, as
  right-drag does, for trackpads. `OrbitControls.js` is no longer bundled.
- Refactored the 3D channel visualizer to consume the shared layout and plane
  helpers; numerical redistribution remains its default.
- The advanced example's nozzle circuits now use an area-preserving helix that
  starts at 20 deg, holding the coolant section constant instead of letting it
  swell with the nozzle radius.

### Deprecated
- `make_engine_3d` and `Engine3DViewer` are deprecated in favour of `Scene`
  and `ThrustChamberComponent`; examples and RL10 post-processing use the scene API.

### Removed
- `solids.TEOS`; use `solids.SilicaDeposit`. It was a constant 0.5 W/m·K with
  no source, and was named after the additive rather than the layer it forms.
- `ChannelPlacement.compute_centerline_radius` and
  `compute_flow_centerline_radius`. A placement owns the circumferential
  packing and the winding; where the section sits radially is resolved from
  the cross section and the wall stack.
- The superseded drawing helpers in `pyskyfire.viz.engine_viz`
  (`get_section_scalars`, `rounded_lobe_arcs_yz`, `rounded_slice_arcs_uv`,
  `rounded_base_height` and their local helpers). The viewer and the exporter
  both build sections from `build_channel_construction`.

### Fixed
- Resolved channel geometry is built once per geometry instead of once per
  solver residual. A single-station section has to be interpolated out of the
  whole-circuit bundle, since it carries resolved fields such as the lobe
  sector radius and the winding rate that nothing can reconstruct from a
  radius alone. Assembling that bundle walks every station of every wall and
  evaluates each design-input function a station at a time, and the wall
  balance reached it on every residual it formed. The bundle is now cached
  against the path it was built from, single stations against the station they
  describe, and each section's boundary integral against the profile it
  measures; the profile setters and the construction steps invalidate all
  three. Results are bit-identical. The RL10 validation case goes from 30 s to
  4 s, and the advanced example's full cycle from 46 s to 36 s. The worst of it
  was an engine with a callable wall thickness, such as the RL10 silver
  insert: its thickness function was called 677,000 times per run, and
  two-thirds of those were for tier B quantities that a rounded section never
  reads.
- Squared channel ribs now point at the engine axis. Their faces used to be
  ruled by the wall normal swept along the channel, which showed up in a
  transverse CAD section as ribs leaning by up to 20 degrees on a sloping
  helical wall. They are ruled by the radial direction instead.
- A helical rib is no longer thinner than asked for. Both of its faces are
  helicoids, and one angular rate is a steeper helix further out, so the
  material between neighbouring channels closes in above the channel floor --
  about 4 percent over the channel height in the advanced example. The
  azimuthal offset is now sized so that the tightest point, at the channel
  roof, is exactly ``rib_width``. ``rib_width`` therefore means the minimum
  material through the solid, which is what a CAD measurement reports.
- Rounded channel lobes are now sized at the radius they actually end up at.
  Each tube is tangent to the sector that separates it from its neighbours and
  keeps that sector's share of the combustion-gas area, instead of being sized
  on the contour and then translated outward. Neighbouring tubes therefore meet
  with the intended amount of material between them in CAD. The RL10 numerical
  regression baseline moved by about one percent as a result; the experimental
  validation cases were unaffected.
- Corrected the rounded full-height conversion, which used
  `h_base*(sin(theta/2) + cos(theta/2))` where the geometry gives
  `h_base*(1 + sin(theta/2))`.
- Corrected squared-channel visualization to account for the ribs and to
  avoid counting hot-side wall layers twice when locating the wetted passage.

---

## [0.3.0] - 30-08-2026

### Added
- A new coupled cooling solver that balances hot-gas convection, conduction
  through multi-material walls, and coolant heat pickup. It supports separate
  wall, heat-flux, and coolant grids, reports convergence at every wall node,
  and returns a richer `RegenResult` for post-processing.
- Enthalpy-based coolant marching for pure fluids, including saturation-state,
  phase, vapour-quality, and boiling-region diagnostics. Mixtures and fluids
  without saturation data retain a temperature-based fallback.
- Standalone and regeneratively coupled subcritical film-cooling analysis based
  on the Grisson model, including liquid-film evaporation and dryout, gaseous
  film entrainment and radiation, and result plots and examples.
- A redesigned cooling-channel geometry system supporting straight, slanted,
  helical, and interleaved circuits; multiple channels per placement leaf;
  variable channel dimensions, and uneven number of channels in each interleaved circuit. 
- Interactive 3D views of the chamber and its actual cooling
  passages.
- Interactive engine-cycle network diagrams, mass-flow Sankey diagrams, and a
  much broader plotting suite for regenerative cooling, film cooling, wall
  layers, coolant phase, hot-gas fields, and engineering reference charts.
- Markdown blocks, responsive layouts, improved tables, and embedded
  interactive content in standalone HTML engineering reports.
- A full Sphinx documentation site with API reference, tutorials, design
  examples, theory explanations, and project contribution guidance.
- Automated tests and CI coverage for package imports, examples, solver
  physics, geometry, serialization, plotting, reports, and validation data.
- End-to-end minimal, advanced engine, mixture-ratio optimization, boiling,
  curvature-correction, and film-cooling examples.
- A comprehensively rebuilt RL10A-3-3A validation case covering engine sizing,
  cycle layout, regenerative cooling, channel geometry, published reference data, and a single consolidated report.
- Citation metadata in `CITATION.cff`, so releases are archived on Zenodo under
  a citable DOI.

### Changed
- Regenerative cooling was overhauled from the legacy heat-exchanger ODE into
  a station-marching coupled solution. That is, film cooling and regenerative cooling is available in the same engine. Coolant energy is advanced on enthalpy
  where possible, wall temperatures come from local nonlinear heat balances,
  and pressure loss follows the full friction-and-acceleration momentum balance. Performance was improved by utilising caching during solve. 
- Hot-side, coolant-side, curvature, channel-rib, and multi-layer wall models
  were revised and made consistent with the physical channel count used by the
  solver, geometry, volume calculations, and visualizations.
- `Aerothermodynamics` was rebuilt around the native `cea` package. Gas states
  are now evaluated on demand along the attached contour, cached by station,
  serializable without native solver handles, and available for equilibrium or
  frozen transport calculations.
- Thrust-chamber contours, circuit placement, cross-sections, and 3D geometry
  generation were substantially reworked for stable interpolation and a common
  representation of real channel centerlines and dimensions.
- Hot-gas plot resolution is now selected when plotting rather than when
  attaching a contour. Plots can use a solved result grid, a uniform node
  count, or explicit coordinates and reuse cached CEA states across properties.
- Visualization dependencies moved to the optional `pyskyfire[viz]` extra, so
  the core installation no longer requires the plotting and meshing stack.
- The supported runtime is now Python 3.12 or newer, and development,
  documentation, and test dependencies are organized into dedicated groups.

### Fixed
- Corrected coolant pressure-drop calculations to include acceleration from
  heating and area changes with the correct sign and magnitude.
- Corrected physical channel counting and coolant-volume calculations when a
  placement leaf represents more than one channel.
- Fixed geometry errors affecting ribs, rounded cross-sections,
  helical passages, contour splines and curvature, wall layers, and their 3D
  representations.
- Improved hot-gas property caching and low-temperature behavior, result
  serialization, report sizing and table layout, and solver diagnostics for
  unconverged wall nodes.

### Removed
- The legacy `regen.solver` implementation and its `steady_heating_analysis`
  entry point; use `coupled_steady_heating_analysis` and the new result model.
- Deprecated `CoolingCircuitGroup` and `WallGroup` wrappers; cooling circuits
  are passed directly to `ThrustChamber`, and wall layers are defined per
  circuit.
- The Cantera combustion path and obsolete combustion/nozzle-solver modules.
  Chemical-equilibrium and hot-gas calculations now use `Aerothermodynamics`
  backed by `cea`.
- The `npts` sizing-constructor argument, stored aerothermodynamic plotting
  grid, and deprecated `Nt` argument. Resolution now belongs to the requested
  solve or visualization.

---

## [0.2.1] - 08-09-2025
## Added
- Various bug fixes

---

## [0.2.0] - 08-09-2025
## Added
- New combustion properties driver CEA_Wrap, replaces RocketCEA and Cantera. The old interface to RocketCEA and Cantera, CombustionTransport, is replaced with Aerothermodynamics. Leverages the large amount of work that has gone into the NASA CEA database. This means that all of NASA CEAs propellants are now supported! Legacy CombustionTransport depricated. 
- New method of hot gas heat transfer estimation in the regen solver based on the new Aerothermodynamics class. 
- Revamped cross sections to use improved class structure.
- Fixed an issue in the CrossSectionSquared class where the rib width was impossible to control. 
- Added a new fluids interface that allows for fuel and oxidizer blends. An example would be an ethanol-water mix.
- Rewrote entire visualisation pipeline. HTML Reports are now the predominant way of viewing results. 
- Added minimal example
- Rewrote RL10 validation case

---

## [0.1.2] - 23-08-2025
### Added
- Packaging and distribution support: `pyskyfire` is now available via pip and PyPI.
- An updated `README.md` explaning the package.
- New logo, which is currently just "pyskyfire" written in the orbitron open-source font. 
- This `CHANGELOG.md` file to track notable changes.

---

## [0.1.0] - 29-06-2025
### Added
- First public release of **pyskyfire**.
- Rocket engine **multi-pass regenerative cooling solver**.
- **Full cycle solver** for expander-cycle engines and related architectures.
- **Visualisation tools** for engine performance and cooling analysis.
