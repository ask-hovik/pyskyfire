"""Design a single-stage impulse oxygen turbine and report it, with a 3D scene.

Run from the repository with:
    python examples/turbine/turbine_report.py --output-dir examples/turbine/output
Requires pyskyfire[impeller,viz]. No browser or display is needed for HTML output.
"""
import argparse
from pathlib import Path

import CoolProp.CoolProp as CP
import numpy as np

from pyskyfire.turbine import AxialTurbine, RowSettings, SectionProfile, Stage, TurbineDuty
from pyskyfire.viz import (Report, Scene, TurbineComponent, plot_turbine_cascade,
                           plot_turbine_meridional, plot_turbine_spanwise, plot_velocity_triangles)


def build_turbine():
    # An oxygen-rich turbine: 26 kg/s of oxygen at 700 K, as it might leave an
    # oxygen-rich preburner or a heat exchanger. The pressures and speed are
    # the right size for a LOX turbopump of this flow, not any one engine's.
    duty = TurbineDuty.from_states(fluid="Oxygen", mdot=26., p0_in=220e5, T0_in=700.,
                                   p_out=140e5, n=30_000)
    return AxialTurbine(
        duty,
        # Supplied, as a turbine this size and duty might reach. The design
        # takes it as given; it sets the work and, through the loss it
        # implies, the density each blade height is sized on.
        efficiency=0.8,
        # One impulse stage: the nozzle takes the whole static pressure drop
        # and the rotor only turns the flow. That leaves little pressure to
        # drive leakage or axial thrust across the rotor, in oxygen at 200 bar.
        # u/c0 = 0.45 is near where an impulse stage works best, and here it
        # also sends the flow out almost without swirl.
        stages=[Stage(flow_coefficient=0.4, reaction=0.0)],
        blade_jet_ratio=0.45,
        rows=RowSettings(
            aspect_ratio=1.0, zweifel=0.8,
            # The nozzle turns the flow from axial to nearly 80 degrees. Its
            # chord has to lean most of the way to the outlet angle before the
            # built vanes' throats leave the flow where the triangles ask.
            stator=SectionProfile(stagger_fraction=0.8),
            # The impulse rotor takes its flow in as steeply as it lets it
            # out, so its inlet is as narrow as its throat. Thin edges and a
            # mid-chord thickness peak keep the passage from narrowing any
            # further between them.
            rotor=SectionProfile.impulse(stagger_fraction=0.55),
            # A shrouded rotor: the tip shroud turns with the blades and forms
            # the outer wall over them, so there is no tip gap to the casing.
            tip_clearance=0.0,
        ),
    )


def _area_ratio(mach, gamma):
    """Ideal-gas ``A/A*``: how much larger than sonic a passage is at this Mach number."""
    return (1/mach)*((2/(gamma + 1))*(1 + 0.5*(gamma - 1)*mach**2))**((gamma + 1)/(2*(gamma - 1)))


def _gamma(state):
    return (CP.PropsSI("Cpmass", "P", state.p, "T", state.T, "Oxygen")
            / CP.PropsSI("Cvmass", "P", state.p, "T", state.T, "Oxygen"))


def _g(value, digits=4):
    return f"{value:.{digits}g}"


def _table(header, rows):
    lines = ["| " + " | ".join(header) + " |", "|" + "---|"*len(header)]
    return "\n".join(lines + ["| " + " | ".join(str(c) for c in row) + " |" for row in rows])


def create_report(output_dir, *, screenshot=False):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    turbine = build_turbine()
    geometry = turbine.definition.build_geometry(span_points=9, chord_points=97)
    duty = turbine.duty

    scene = Scene(title="Axial turbine", show_axes=True)
    scene.add(TurbineComponent(geometry), name="turbine")
    scene.save_html(output_dir / "turbine-scene.html")

    report = Report("Single-stage impulse oxygen turbine")
    tab = report.add_tab("Turbine")
    tab.add_markdown(f"""
        # Single-stage impulse oxygen turbine
        {_g(duty.mdot)} kg/s of oxygen expanding from {_g(duty.inlet.p/1e5)} bar and
        {_g(duty.inlet.T)} K to {_g(duty.outlet.p/1e5)} bar, at {duty.n:.0f} rpm.
        The grey row is the nozzle vanes, the red row the rotor blades. The
        rotor is shrouded, so the casing over it turns with it. The translucent
        surfaces are the hub and casing. Rotate and zoom to look at the blades.
    """)
    tab.add_iframe(scene.data_url, caption="Blade rows generated with pyTurbo Aero",
                   height="clamp(360px, calc(100vh - 260px), 720px)")

    design = report.add_tab("Design")
    design.add_table({
        "Mass flow [kg/s]": duty.mdot,
        "Inlet total pressure [bar]": duty.inlet.p/1e5,
        "Inlet total temperature [K]": duty.inlet.T,
        "Outlet total pressure [bar]": duty.outlet.p/1e5,
        "Pressure ratio": duty.pressure_ratio,
        "Isentropic head [kJ/kg]": duty.isentropic_head/1e3,
        "Spouting velocity c0 [m/s]": duty.spouting_velocity,
        "Shaft speed [rpm]": duty.n,
        "Specific speed Ns": duty.specific_speed,
        "Mean diameter [mm]": 2e3*turbine.mean_radius,
        "Blade speed u [m/s]": turbine.blade_speed,
        "u/c0": turbine.blade_speed/duty.spouting_velocity,
        "Efficiency (supplied)": turbine.efficiency,
        "Shaft power [kW]": duty.mdot*turbine.work/1e3,
        "Exit swirl [m/s]": turbine.exit_swirl,
        "Exit flow angle [deg]": np.degrees(turbine.exit_flow_angle),
    }, caption="Duty and meanline")
    design.add_markdown("""
        A pressure ratio of about 1.57 across 700 K oxygen gives a spouting
        velocity of about 400 m/s. One impulse stage takes it at a
        blade-to-jet ratio of 0.45, near where an impulse stage works best,
        which at 30 000 rpm puts the mean diameter at about 116 mm.

        At zero reaction the nozzle takes the whole static enthalpy drop: the
        flow leaves it at Mach 0.74, almost tangentially. The rotor turns the
        flow back through about 135° at constant relative speed, and it
        leaves nearly axially. The rotor still loses some static pressure,
        because half the stage's losses are charged to it.

        The meanline runs at one mean radius. Blade heights follow from
        continuity at each station's density. That density comes from the
        station's enthalpy and the entropy the supplied efficiency implies,
        shared evenly between the two rows.
    """)
    design.add_markdown(_table(
        ("Station", "c_x [m/s]", "c_u [m/s]", "α [°]", "β [°]", "p [bar]", "p0 [bar]", "T [K]",
         "ρ [kg/m³]", "M", "M_rel", "Height [mm]"),
        [(s.name, _g(s.c_x), _g(s.c_u), f"{np.degrees(s.alpha):.1f}", f"{np.degrees(s.beta):.1f}",
          _g(s.static.p/1e5, 5), _g(s.p0/1e5, 5), _g(s.static.T, 5), _g(s.static.rho),
          _g(s.mach, 3), _g(s.mach_relative, 3), _g(turbine.station_height(s)*1e3))
         for s in turbine.stations]))
    design.add_markdown(_table(
        ("Stage", "Reaction", "Flow coefficient", "Loading Δh0/u²", "Work [kJ/kg]"),
        [(k + 1, _g(s.reaction), _g(s.flow_coefficient), _g(s.loading), _g(s.work/1e3))
         for k, s in enumerate(turbine.stages)]))
    design.add_figure(plot_velocity_triangles(turbine))

    rows = report.add_tab("Blade rows")
    rows.add_markdown(_table(
        ("Row", "Blades", "Axial chord [mm]", "Chord [mm]", "Height [mm]", "Metal angles [°]",
         "Stagger [°]", "Pitch [mm]", "Throat [mm]", "Gauging [°]", "Gauging − metal [°]", "t_max/c"),
        [(rg.row.name, rg.row.blade_count, _g(rg.row.axial_chord*1e3), _g(rg.row.chord*1e3),
          f"{_g(rg.row.height_in*1e3)} → {_g(rg.row.height_out*1e3)}",
          f"{np.degrees(rg.row.inlet_metal_angle):.1f} → {np.degrees(rg.row.outlet_metal_angle):.1f}",
          f"{np.degrees(rg.row.stagger):.1f}", _g(rg.at_mid_span(rg.pitch)*1e3),
          _g(rg.at_mid_span(rg.throat)*1e3), f"{np.degrees(rg.at_mid_span(rg.gauging_angle)):.1f}",
          f"{np.degrees(rg.at_mid_span(rg.gauging_angle - rg.outlet_metal_angle)):+.1f}",
          _g(rg.at_mid_span(rg.max_thickness)/rg.row.chord, 3))
         for rg in geometry.rows]))
    rows.add_markdown("""
        Metal angles, stagger and chord are the drawing's values at mean radius.
        Pitch, throat and gauging angle are measured on the built mid-span
        section. The gauging angle, `acos(throat/pitch)`, is roughly where the
        flow leaves each row, so its distance from the outlet metal angle shows
        how the drawn blade differs from the triangles it was drawn for.

        Stagger decides most of that distance. The nozzle turns the flow to
        nearly 80°, and its chord leans 80 % of the way to the outlet angle
        (`stagger_fraction=0.8`) to bring its gauging angle within about 2°.
        The rotor is a symmetric impulse blade (`SectionProfile.impulse`),
        with its chord just past the mean of its metal angles.
    """)
    choking = []
    st = turbine.stations
    for rg, inlet, m_in, m_out in ((geometry.rows[0], st[0].static, st[0].mach, st[1].mach),
                                   (geometry.rows[1], st[1].static, st[1].mach_relative,
                                    st[2].mach_relative)):
        inlet_area = rg.at_mid_span(rg.inlet_width)*rg.row.height_in
        sonic_area = inlet_area/_area_ratio(m_in, _gamma(inlet))
        smallest = rg.at_mid_span(rg.narrowest*rg.narrowest_height)
        choking.append((rg.row.name, _g(inlet_area*1e6), _g(smallest*1e6),
                        _g(rg.at_mid_span(rg.narrowest_at), 2), _g(rg.pinch.min(), 3),
                        _g(m_in, 3), _g(m_out, 3), _g(sonic_area*1e6), _g(smallest/sonic_area, 3)))
    rows.add_markdown(_table(
        ("Row", "Inlet area / passage [mm²]", "Smallest area [mm²]", "At x/c_ax",
         "Pinch (worst span)", "Inlet Mach", "Exit Mach", "Sonic area [mm²]",
         "Smallest / sonic"), choking))
    rows.add_markdown("""
        **Does the flow choke?** Mach numbers are absolute for the nozzle and
        relative for the rotor. Areas are per blade passage: gap times blade
        height, taken where the gap is. *Sonic area* is how small the passage
        could get before the flow entering it reached Mach 1: the inlet area
        over the ideal-gas `A/A*` at the inlet Mach number, with γ from
        CoolProp at the row inlet. The smallest area stays above it in both
        rows, so neither chokes. This is a one-dimensional, loss-free estimate
        at mid-span; losses shrink the margin a little.
        The nozzle comes closest, since it is meant to accelerate its flow to
        Mach 0.74. The rotor, at a relative Mach of about 0.4, has room to
        spare. A pinch of one means the passage narrows no further than the
        smaller of its inlet width and throat.
    """)
    rows.add_figure(plot_turbine_meridional(geometry))
    rows.add_figure(plot_turbine_cascade(geometry))
    rows.add_figure(plot_turbine_spanwise(geometry),
                    caption="Free-vortex twist: the rotor twists more than the nozzle, since blade speed grows with radius.")

    report_path = output_dir / "turbine-report.html"
    report.save_html(str(report_path))
    if screenshot:
        plotter = scene.to_plotter(off_screen=True)
        try:
            plotter.window_size = (1400, 1000)
            plotter.screenshot(str(output_dir / "turbine-preview.png"))
        finally:
            plotter.close()
    print(f"Shaft power {duty.mdot*turbine.work/1e3:.0f} kW, "
          f"blades {[r.blade_count for r in turbine.definition.rows]}")
    print(f"Report: {report_path.resolve()}")
    return turbine, geometry, report_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent / "output")
    parser.add_argument("--screenshot", action="store_true", help="Also render a PyVista PNG preview")
    args = parser.parse_args()
    create_report(args.output_dir, screenshot=args.screenshot)
