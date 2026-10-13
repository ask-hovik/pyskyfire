"""Design a LOX/methane impinging injector with a film ring, drill it, and report it.

The injector has unlike doublets in four rings, fed from a main fuel plenum and
an oxidizer plenum, and an outer ring of single fuel jets aimed shallowly at the
chamber wall, fed from a separate film ring at a higher pressure. Holes are sized for a design point, then
rounded to the nearest available drill. The drilled injector is rated two ways:
the flows it passes at the design plenum pressures, and the plenum pressures
it needs to pass the design flows.

Run from the repository with:
    python examples/injector/injector_report.py --output-dir examples/injector/output
Requires pyskyfire[viz]. No browser or display is needed for HTML output.
"""
import argparse
from dataclasses import dataclass
import math
from pathlib import Path

import numpy as np
import plotly.graph_objects as go

from pyskyfire.ducts import Orifice
from pyskyfire.injector import ElementSet, FeedCondition, ImpingingElement, Injector, Jet
from pyskyfire.thermo import FluidState, PHASE_UNKNOWN
from pyskyfire.viz import InjectorComponent, Report, Scene, plot_injector_face

# Feeds. A feed can be any hashable object, such as a ducts.Manifold; names
# stand in here because the plenums themselves are not modelled.
MAIN, FILM, OX = "Main fuel plenum", "Film ring", "Oxidizer plenum"
FUEL, OXIDIZER = "Methane", "Oxygen"

# The operating point the injector is designed for.
P_C = 60e5                 # face pressure [Pa]
MDOT_OX = 6.0              # [kg/s]
MDOT_FU = 1.9              # all fuel, film included [kg/s]
FILM_SHARE = 0.10          # of the fuel, through the film ring
T_FU, T_OX = 120.0, 95.0   # plenum temperatures [K]

# Stiffness: plenum pressure over face pressure, per feed. The film ring is
# fed from a higher-pressure source than the main plenum.
PLENUM_RATIO = {MAIN: 1.15, FILM: 1.30, OX: 1.20}

CHAMBER_RADIUS = 0.065     # [m]

# Film cooling: single fuel jets from holes on a ring inboard of the wall,
# leaning outward so they strike the wall at a shallow angle and spread along
# it as a film. Each jet's impingement point is where it meets the wall.
N_FILM = 30
FILM_HOLE_RADIUS = 0.058   # [m]
FILM_LEAN = np.radians(15)  # from the chamber axis
FILM_STRIKE = (CHAMBER_RADIUS - FILM_HOLE_RADIUS)/np.tan(FILM_LEAN)  # downstream of the face [m]
DRILL_STEP = 0.05e-3       # drills available in 0.05 mm steps [m]

# Line colour for the one-series sweep; slot 1 of the validated palette.
SERIES = "#2a78d6"


def state(p, T):
    nan = float("nan")
    return FluidState(T=T, p=p, h=nan, quality=nan, rho=nan, phase=PHASE_UNKNOWN, T_sat=nan)


def rings(counts, radii, stagger=True):
    """Element centres on concentric rings, every other ring rotated half a pitch."""
    r = np.concatenate([np.full(n, radius) for n, radius in zip(counts, radii)])
    theta = np.concatenate([2*np.pi*(np.arange(n) + 0.5*(k % 2 if stagger else 0))/n
                            for k, n in enumerate(counts)])
    return r, theta


@dataclass
class Design:
    """The designed and drilled injectors, with the numbers that made them."""

    designed: Injector
    drilled: Injector
    diameters: dict      # feed -> (designed, drilled) diameter [m]
    inclinations: dict   # feed -> jet inclination [rad]


def _element_pair(fuel_hole, ox_hole, fuel_inclination, ox_inclination):
    return ImpingingElement([Jet(MAIN, fuel_hole, fuel_inclination, azimuth=0.0),
                             Jet(OX, ox_hole, ox_inclination, azimuth=np.pi)],
                            impingement_height=6e-3)


def _film_jet(hole):
    # One jet, placed half a turn round from the element's +y, so with the
    # element centre on the wall the hole sits inboard and leans outward.
    return ImpingingElement([Jet(FILM, hole, FILM_LEAN, azimuth=np.pi)],
                            impingement_height=FILM_STRIKE)


def _injector(pair, film, main_layout, film_layout):
    return Injector([
        # Jets in the radial plane, so each spray fan opens tangentially.
        ElementSet(pair, *main_layout, orientation=0.0, name="Main doublets"),
        # Film elements are centred on the wall, where each jet strikes it.
        ElementSet(film, *film_layout, orientation=0.0, name="Film jets"),
    ])


def build_design():
    main_layout = rings(counts=[6, 12, 18, 24], radii=[0.012, 0.024, 0.036, 0.048])
    film_layout = rings(counts=[N_FILM], radii=[CHAMBER_RADIUS])
    n_main, n_film = len(main_layout[0]), len(film_layout[0])

    plenum = {feed: P_C*ratio for feed, ratio in PLENUM_RATIO.items()}
    fuel_main = Orifice.sized(mass_flow=MDOT_FU*(1 - FILM_SHARE)/n_main,
                              inlet_state=state(plenum[MAIN], T_FU), dp=plenum[MAIN] - P_C,
                              fluid=FUEL, cd=0.75, length=6e-3)
    ox = Orifice.sized(mass_flow=MDOT_OX/n_main, inlet_state=state(plenum[OX], T_OX),
                       dp=plenum[OX] - P_C, fluid=OXIDIZER, cd=0.75, length=7e-3)
    film = Orifice.sized(mass_flow=MDOT_FU*FILM_SHARE/n_film,
                         inlet_state=state(plenum[FILM], T_FU), dp=plenum[FILM] - P_C,
                         fluid=FUEL, cd=0.72, length=4e-3)

    # Lean the oxidizer jet less than the fuel jet, so the two jets' sideways
    # momenta cancel and the resultant runs down the element axis:
    # mdot_o v_o sin(a_o) = mdot_f v_f sin(a_f).
    v_fu = fuel_main.solve(inlet_state=state(plenum[MAIN], T_FU),
                           mass_flow=MDOT_FU*(1 - FILM_SHARE)/n_main, fluid=FUEL).velocity
    v_ox = ox.solve(inlet_state=state(plenum[OX], T_OX), mass_flow=MDOT_OX/n_main,
                    fluid=OXIDIZER).velocity
    a_fu = np.radians(30)
    a_ox = math.asin(MDOT_FU*(1 - FILM_SHARE)*v_fu*math.sin(a_fu)/(MDOT_OX*v_ox))

    designed = _injector(_element_pair(fuel_main, ox, a_fu, a_ox), _film_jet(film),
                         main_layout, film_layout)

    def drill(hole):
        d = DRILL_STEP*round(hole.hydraulic_diameter/DRILL_STEP)
        return Orifice.circular(d, cd=hole.cd, length=hole.length)

    drilled = _injector(_element_pair(drill(fuel_main), drill(ox), a_fu, a_ox),
                        _film_jet(drill(film)), main_layout, film_layout)
    diameters = {MAIN: (fuel_main.hydraulic_diameter, drill(fuel_main).hydraulic_diameter),
                 OX: (ox.hydraulic_diameter, drill(ox).hydraulic_diameter),
                 FILM: (film.hydraulic_diameter, drill(film).hydraulic_diameter)}
    return Design(designed, drilled, diameters, {MAIN: a_fu, OX: a_ox, FILM: FILM_LEAN})


FLUID = {MAIN: FUEL, FILM: FUEL, OX: OXIDIZER}
TEMPERATURE = {MAIN: T_FU, FILM: T_FU, OX: T_OX}
DESIGN_FLOW = {MAIN: MDOT_FU*(1 - FILM_SHARE), FILM: MDOT_FU*FILM_SHARE, OX: MDOT_OX}


def at_pressures(injector, ratios=PLENUM_RATIO):
    """The flows the injector passes with each plenum at ``ratio·p_c``."""
    return injector.solve({f: FeedCondition(FLUID[f], TEMPERATURE[f], pressure=P_C*ratios[f])
                           for f in injector.feeds}, p_face=P_C, oxidizer=[OX])


def at_flows(injector):
    """The plenum pressures the injector needs to pass the design flows."""
    return injector.solve({f: FeedCondition(FLUID[f], TEMPERATURE[f], mass_flow=DESIGN_FLOW[f])
                           for f in injector.feeds}, p_face=P_C, oxidizer=[OX])


def film_sweep(injector, ratios):
    """Film share of the fuel as the film ring pressure varies."""
    shares = []
    for ratio in ratios:
        result = at_pressures(injector, {**PLENUM_RATIO, FILM: ratio})
        shares.append(result.feeds[FILM].mass_flow/result.fuel_mass_flow)
    return np.asarray(shares)


def plot_film_sweep(ratios, shares, design_ratio, design_share):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=ratios, y=100*shares, mode="lines", name="Film share",
                             line=dict(color=SERIES, width=2),
                             hovertemplate="Film ring at %{x:.3f}·p_c<br>"
                                           "Film share %{y:.2f} %<extra></extra>"))
    fig.add_trace(go.Scatter(x=[design_ratio], y=[100*design_share], mode="markers+text",
                             marker=dict(size=10, color=SERIES, line=dict(color="white", width=2)),
                             text=[f"Design ring pressure: {100*design_share:.1f} %"], textposition="top left",
                             textfont=dict(color="#3d3d3a"), showlegend=False,
                             hovertemplate="Design point<extra></extra>"))
    fig.update_xaxes(title_text="Film ring pressure / face pressure", gridcolor="#e4e4e0")
    fig.update_yaxes(title_text="Film share of the fuel flow [%]", gridcolor="#e4e4e0",
                     rangemode="tozero")
    fig.update_layout(title="Film share of the drilled injector against film ring pressure",
                      template="plotly_white", showlegend=False, hovermode="closest")
    return fig


def _g(value, digits=4):
    return f"{value:.{digits}g}"


def _table(header, rows):
    lines = ["| " + " | ".join(header) + " |", "|" + "---|"*len(header)]
    return "\n".join(lines + ["| " + " | ".join(str(c) for c in row) + " |" for row in rows])


def create_report(output_dir, *, screenshot=False):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    design = build_design()
    drilled = design.drilled
    rated = at_pressures(drilled)
    needed = at_flows(drilled)
    checks = drilled.checks(chamber_radius=CHAMBER_RADIUS)
    main, film = rated.sets["Main doublets"], rated.sets["Film jets"]

    scene = Scene(title="Impinging injector", show_axes=True)
    scene.add(InjectorComponent(drilled, result=rated, face_radius=CHAMBER_RADIUS,
                                show_impingement=True), name="injector")
    scene.save_html(output_dir / "injector-scene.html")

    report = Report("LOX/methane impinging injector")
    tab = report.add_tab("Injector")
    tab.add_markdown(f"""
        # LOX/methane impinging injector
        {len(drilled.sets['Main doublets'].r)} unlike doublets in four rings
        take the main fuel and all the oxidizer. An outer ring of
        {len(drilled.sets['Film jets'].r)} single fuel jets, leaning
        {np.degrees(FILM_LEAN):.0f}° outward, strikes the wall
        {FILM_STRIKE*1e3:.0f} mm downstream of the face and spreads along it
        as a film. The film ring is fed from its own source at
        {PLENUM_RATIO[FILM]:.2f}·p_c, higher than the main fuel plenum's
        {PLENUM_RATIO[MAIN]:.2f}·p_c. The scene shows the wetted surface: each
        hole's bore, coloured by feed, ending at the face plate. The grey
        arrows are not hardware; each points from a hole to the impingement
        point it is aimed at. Colour the scene by *Hole mass flow* or *Jet
        velocity* to see the drilled injector at the design plenum pressures.
        The grey disc is the chamber, {2e3*CHAMBER_RADIUS:.0f} mm across.
    """)
    tab.add_iframe(scene.data_url, caption="The drilled injector, seen into the chamber (+x)",
                   height="clamp(360px, calc(100vh - 260px), 720px)")
    tab.add_figure(plot_injector_face(drilled, face_radius=CHAMBER_RADIUS),
                   caption="Hover a hole for its set, element and diameter.")

    sizing = report.add_tab("Design")
    sizing.add_markdown(f"""
        ## Design point and choices
        The injector is designed to deliver {MDOT_OX:g} kg/s of oxygen at
        {T_OX:g} K and {MDOT_FU:g} kg/s of methane at {T_FU:g} K into a
        {P_C/1e5:g} bar face, {100*FILM_SHARE:g} % of the fuel through the film
        ring: an overall mixture ratio of {MDOT_OX/MDOT_FU:.3g}. Every choice
        below is supplied, not solved for. Each hole is sized for its share
        of its feed's flow at its feed's pressure drop, with
        `Orifice.sized`, and then rounded to the nearest
        {DRILL_STEP*1e3:g} mm drill.
    """)
    sizing.add_markdown(_table(
        ("Feed", "Fluid", "Plenum / p_c", "Cd", "L [mm]", "Holes", "Designed d [mm]",
         "Drilled d [mm]", "Inclination [°]"),
        [(f, FLUID[f], f"{PLENUM_RATIO[f]:.2f}", _g(drilled.equivalent_orifice(f).cd, 3),
          _g(drilled.sets["Main doublets" if f != FILM else "Film jets"]
             .element.orifices(f)[0].length*1e3, 3),
          drilled.hole_count(f), _g(design.diameters[f][0]*1e3, 4),
          _g(design.diameters[f][1]*1e3, 4), f"{np.degrees(design.inclinations[f]):.1f}")
         for f in drilled.feeds]))
    sizing.add_markdown(f"""
        The fuel jet of each unlike doublet leans 30° towards the element
        axis. The oxidizer jet carries more momentum, so it leans only
        {np.degrees(design.inclinations[OX]):.1f}°: then the two jets'
        sideways momenta cancel and the spray leaves along the axis. Equal
        angles would have tipped it off the axis towards the fuel hole's
        side.

        ## Mixing, drilled injector at the design pressures
    """)
    sizing.add_table({
        "Element mixture ratio, main doublets": main.diagnostics.mixture_ratio,
        "Rupe ratio ρ_o v_o² d_o / ρ_f v_f² d_f": main.diagnostics.rupe,
        "Momentum-flux ratio ρ_o v_o² / ρ_f v_f²": main.diagnostics.momentum_flux_ratio,
        "Resultant angle from axis, main doublets [°]":
            np.degrees(main.diagnostics.resultant_angle),
        "Film jet angle from axis [°]": np.degrees(film.diagnostics.resultant_angle),
        "Fuel jet velocity [m/s]": rated.feeds[MAIN].velocity,
        "Oxidizer jet velocity [m/s]": rated.feeds[OX].velocity,
        "Film jet velocity [m/s]": rated.feeds[FILM].velocity,
    }, caption="Reported, not enforced")
    stronger, weaker = ("oxidizer", "fuel") if main.diagnostics.rupe > 1 else ("fuel", "oxidizer")
    sizing.add_markdown(f"""
        Rupe found mixing most uniform with his ratio near one. At
        {main.diagnostics.rupe:.2f} the {stronger} jets overpower the
        {weaker} jets. Since ρv²d grows as about Δp^0.75 at a fixed flow,
        raising the {weaker} pressure drop or lowering the {stronger} one
        closes the gap. The stiffness split here is kept for its margin on
        each side, and the ratio is left as a report.

        ## Layout checks
    """)
    sizing.add_markdown("""
        The film jets are aimed at the wall on purpose, so their impingement
        points sit on it and leave no clearance. The holes themselves stay
        clear of it.
    """)
    sizing.add_table({
        "Narrowest web between holes [mm]": checks.min_web*1e3,
        "Between sets": " and ".join(dict.fromkeys(checks.min_web_sets)),
        "Overlapping hole pairs": checks.overlaps,
        "Outermost hole edge [mm]": checks.outer_radius*1e3,
        "Outermost impingement point [mm]": checks.target_radius*1e3,
        "Clearance from the holes to the wall [mm]": checks.wall_clearance*1e3,
        # Rounded so floating-point noise on the wall-aimed jets reads as zero.
        "Clearance from the impingement points to the wall [mm]":
            round(checks.target_clearance*1e3, 6) + 0.0,
        "Shortest hole L/d": checks.min_length_to_diameter,
    })

    drilling = report.add_tab("As drilled")
    drilling.add_markdown("""
        ## What rounding the holes to drills did
        Two questions, one injector. *At design pressures*: each plenum is
        held at its design pressure and the flows are whatever the drilled
        holes pass. *At design flows*: each feed passes its design flow and
        the plenum pressure is whatever it takes. The design column is what
        the holes were sized for.
    """)
    drilling.add_markdown(_table(
        ("Feed", "Design flow [kg/s]", "Flow at design pressures [kg/s]", "Change [%]",
         "Design plenum / p_c", "Plenum for design flow / p_c"),
        [(f, _g(DESIGN_FLOW[f]), _g(rated.feeds[f].mass_flow),
          f"{100*(rated.feeds[f].mass_flow/DESIGN_FLOW[f] - 1):+.2f}",
          f"{PLENUM_RATIO[f]:.3f}", f"{needed.feeds[f].plenum_pressure/P_C:.4f}")
         for f in drilled.feeds]))
    drilling.add_table({
        "Mixture ratio, design": MDOT_OX/MDOT_FU,
        "Mixture ratio, drilled at design pressures": rated.mixture_ratio,
        "Film share, design [%]": 100*FILM_SHARE,
        "Film share, drilled at design pressures [%]":
            100*rated.feeds[FILM].mass_flow/rated.fuel_mass_flow,
    })

    ratios = np.linspace(1.10, 1.50, 41)
    shares = film_sweep(drilled, ratios)
    drilled_share = rated.feeds[FILM].mass_flow/rated.fuel_mass_flow
    ring = report.add_tab("Film ring")
    ring.add_markdown(f"""
        ## Setting the film share with the ring pressure
        The film ring and the main fuel plenum are separate feeds, so the
        film holes' flow depends on the film ring's own pressure and nothing
        else. Holding the main plenums at their design pressures and moving
        the film ring shows how the film share responds. Across this range
        it runs from {100*shares[0]:.1f} % to {100*shares[-1]:.1f} % of the
        fuel. How the ring pressure is set, and what it costs the pump, is
        the feed system's question, not the injector's. The marked point is
        the drilled injector at the design ring pressure: {100*drilled_share:.2f} %
        against the {100*FILM_SHARE:g} % the holes were sized for, since
        rounding to drills changed both the film and the main fuel holes.
    """)
    ring.add_figure(plot_film_sweep(ratios, shares, PLENUM_RATIO[FILM], drilled_share))
    ring.add_markdown(_table(("Film ring / p_c", "Film share [%]"),
                             [(f"{r:.2f}", f"{100*s:.2f}") for r, s in
                              zip(ratios[::5], shares[::5])]))

    report_path = output_dir / "injector-report.html"
    report.save_html(str(report_path))
    if screenshot:
        plotter = scene.to_plotter(off_screen=True)
        try:
            plotter.window_size = (1200, 1000)
            reach = CHAMBER_RADIUS
            plotter.camera_position = [(2.2*reach, -1.6*reach, 1.4*reach), (0.0, 0.0, 0.0),
                                       (1, 0, 0)]
            plotter.screenshot(str(output_dir / "injector-preview.png"))
        finally:
            plotter.close()
    print(f"Drilled injector at design pressures: MR {rated.mixture_ratio:.3f}, "
          f"film share {100*rated.feeds[FILM].mass_flow/rated.fuel_mass_flow:.2f} %")
    print(f"Report: {report_path.resolve()}")
    return design, report_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent / "output")
    parser.add_argument("--screenshot", action="store_true", help="Also render a PyVista PNG preview")
    args = parser.parse_args()
    create_report(args.output_dir, screenshot=args.screenshot)
