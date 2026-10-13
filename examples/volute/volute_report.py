"""Design three volutes, compare their losses, and report each in its own 3D scene.

Two radial pump-discharge volutes at different flow conditions, and one
cooling-jacket inlet manifold whose inlet plane is inclined 45 degrees. All
three are sized on conservation of angular momentum at SP = 1 and analysed with
Aungier's meanline loss model.

Each also carries a different exit-cone outlet, to show the range the shape
can take: one runs tangentially straight off the full-collection plane, one
curves all the way round to point radially outward, and one takes an
intermediate angle while also lifting out of the volute's own plane.

Run from the repository with:
    python examples/volute/volute_report.py --output-dir examples/volute/output
Requires pyskyfire[viz]. No browser or display is needed for HTML output.
"""
import argparse
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import plotly.graph_objects as go

import pyskyfire as psf
from pyskyfire.ducts import (ConstantAngularMomentum, EllipticalSlottedSection,
                             ExitCone, Volute, VoluteInlet)
from pyskyfire.thermo.state import PHASE_UNKNOWN, FluidState
from pyskyfire.viz import Report, Scene, VoluteComponent

# Categorical slots, validated for colour-vision deficiency against a light
# surface. Series identity is also carried by the legend and by the tables.
SERIES = ("#2a78d6", "#eb6834", "#1baf7a")
# The categorical theme's first six slots, in its own fixed order.
MECHANISM_COLORS = {"meridional": "#2a78d6", "tangential": "#eb6834",
                    "friction": "#1baf7a", "exit_cone": "#eda100",
                    "exit_friction": "#e87ba4", "exit_bend": "#008300"}
MECHANISM_LABELS = {"meridional": "Meridional dump", "tangential": "Tangential dump",
                    "friction": "Skin friction", "exit_cone": "Exit cone diffusion",
                    "exit_friction": "Duct friction", "exit_bend": "Duct bend"}


@dataclass
class Case:
    """One designed volute together with its operating point and result."""

    name: str
    note: str
    outlet_note: str
    volute: Volute
    result: object
    fluid: str
    T: float
    p: float
    mass_flow: float

    @property
    def exit_diameter(self):
        return self.volute.exit_cone.final_diameter

    @property
    def outer_radius(self):
        """Largest distance any rendered vertex reaches from the axis [m].

        Used only to frame each scene's screenshot, so a coarse sample of both
        sweeps is enough; the exit cone can reach well past the bore itself
        once it curves toward radial.
        """
        points = [self.volute.sweep(n_points=6).reshape(-1, 3)]
        if self.volute.exit_geometry is not None:
            points.append(
                self.volute.sweep_exit_cone(n_points=6, n_stations=9).reshape(-1, 3)
            )
        return float(np.linalg.norm(np.concatenate(points), axis=-1).max())


def design(name, note, *, fluid, T, p, mass_flow, inlet_radius, inlet_width,
           swirl_velocity, fillet_radius, throat_length, exit_cone,
           area_law=None, rotation="ccw", aspect_ratio=1.0,
           max_aspect_ratio=4.0, tilt=0.0, roughness=6.3e-6):
    """Size a volute for a duty and solve it at that duty.

    The meridional velocity is not a free input: the mass flow has to pass
    through the inlet slot, so ``C_m`` follows from the duty and the slot the
    feeding hardware presents. Only the swirl is a design choice here, and it
    is what the pump or the tangential feed delivers.
    """
    transport = psf.thermo.FluidTransport(psf.thermo.Fluid(name, [fluid], [1.0]))
    density = transport.get_rho(T, p)
    meridional = mass_flow / (density * 2.0 * math.pi * inlet_radius * inlet_width)
    volute = Volute(
        EllipticalSlottedSection(inlet_width, fillet_radius, throat_length,
                                 aspect_ratio=aspect_ratio,
                                 max_aspect_ratio=max_aspect_ratio),
        VoluteInlet(radius=inlet_radius, C_m=meridional, C_u=swirl_velocity),
        area_law if area_law is not None else ConstantAngularMomentum(),
        rotation=rotation,
        exit_cone=exit_cone,
        tilt=tilt,
        n_stations=121,
    )
    state = FluidState(T=T, p=p, h=math.nan, quality=math.nan, rho=density,
                       phase=PHASE_UNKNOWN, T_sat=math.nan)
    result = volute.solve(inlet_state=state, mass_flow=mass_flow,
                          transport=transport, roughness=roughness)
    outlet_note = _outlet_note(volute.exit_geometry)
    return Case(name, note, outlet_note, volute, result, fluid, T, p, mass_flow)


def _outlet_note(exit_geometry):
    if exit_geometry.is_straight:
        return "Tangential: runs straight off the full-collection plane."
    turn = math.degrees(exit_geometry.turn_angle)
    if exit_geometry.second_angle == 0.0:
        return f"Curves {turn:.0f}° in-plane to point radially outward."
    lift = "up" if exit_geometry.second_angle > 0.0 else "down"
    return (f"Curves {turn:.0f}° while lifting {lift} "
           f"{abs(math.degrees(exit_geometry.second_angle)):.0f}° out of the plane.")


def build_cases():
    return [
        design(
            "LOX pump discharge",
            "Radial, circular bore. A dense propellant and a narrow discharge "
            "slot, so the flow enters shallow and the spiral grows slowly.",
            fluid="O2", T=100.0, p=9.0e6, mass_flow=32.0,
            inlet_radius=0.055, inlet_width=0.010, swirl_velocity=48.0,
            fillet_radius=0.0025, throat_length=0.002,
            aspect_ratio=1.0, max_aspect_ratio=4.0,
            exit_cone=ExitCone(final_diameter=0.04, first_angle="radial",
                               arc_length_parameter=2.5),
        ),
        design(
            "Methane pump discharge",
            "Radial, and stretched across the slot to stay axially compact. A "
            "third of the density at a similar head, so the same duty enters "
            "at a steeper angle and needs a fatter passage.",
            fluid="CH4", T=125.0, p=1.2e7, mass_flow=11.0,
            inlet_radius=0.048, inlet_width=0.008, swirl_velocity=40.0,
            fillet_radius=0.002, throat_length=0.002,
            aspect_ratio=1.6, max_aspect_ratio=5.0,
            exit_cone=ExitCone(final_diameter=0.038, first_angle="tangential",
                               arc_length_parameter=2.0),
        ),
        design(
            "Coolant inlet manifold",
            "Inclined 45 degrees to follow the nozzle wall it feeds, and "
            "reaching outward rather than across. Modest swirl at a large "
            "radius gives the steepest inlet of the three.",
            fluid="CH4", T=130.0, p=1.2e7, mass_flow=11.0,
            inlet_radius=0.090, inlet_width=0.008, swirl_velocity=16.0,
            fillet_radius=0.002, throat_length=0.002,
            aspect_ratio=1.0, max_aspect_ratio=4.0,
            exit_cone=ExitCone(final_diameter=0.05, first_angle=math.radians(30.0),
                               second_angle=math.radians(35.0), arc_length_parameter=2.5),
            tilt=math.radians(45.0),
        ),
    ]


def slug(name):
    return name.replace(" ", "-").lower()


def build_scene(case, color=None):
    """One scene per volute, each on its own machine axis at the origin.

    Kept separate rather than lined up side by side so that each passage gets
    its own colour scale. Sharing one across all three would flatten the two
    slower passages against the fastest, and they differ by a factor of three.
    """
    scene = Scene(title=case.name, color_by="velocity", show_axes=True)
    node = scene.add(
        VoluteComponent(case.volute, case.result, n_points=16, exit_points=17),
        name=slug(case.name),
    )
    if color is not None:
        node.color = color
    return scene


def loss_figure(cases):
    """Stacked bars: what each case's total loss is made of.

    Stacked rather than grouped because the six terms sum to the total that
    matters, and because two of them are legitimately zero in some cases --
    a grouped zero reads as a missing series, a stacked one simply has no
    height while the total stays honest.
    """
    figure = go.Figure()
    for key, label in MECHANISM_LABELS.items():
        figure.add_bar(
            name=label, x=[case.name for case in cases],
            y=[case.result.loss_breakdown()[key] for case in cases],
            marker_color=MECHANISM_COLORS[key],
            # A thin surface-coloured gap keeps adjacent segments apart.
            marker_line=dict(color="white", width=2),
        )
    figure.update_layout(
        barmode="stack", bargap=0.45,
        title="Loss coefficients, referred to the inlet velocity head",
        yaxis_title="Loss coefficient ω [-]", template="plotly_white",
        legend_title_text="Mechanism",
        annotations=[dict(
            text="Tangential dump is zero throughout: sizing at SP = 1 removes it. "
                 "A straight duct has no bend to charge.",
            xref="paper", yref="paper", x=0.0, y=1.06, showarrow=False,
            xanchor="left", font=dict(size=12, color="#52514e"),
        )],
    )
    figure.update_yaxes(gridcolor="#e8e8e6", zerolinecolor="#d0d0cd")
    figure.update_xaxes(showgrid=False)
    return figure


def _azimuth_figure(cases, values, title, y_title):
    figure = go.Figure()
    for case, color in zip(cases, SERIES, strict=True):
        figure.add_scatter(
            x=np.degrees(case.volute.geometry.theta), y=values(case),
            name=case.name, mode="lines", line=dict(color=color, width=2),
        )
    figure.update_layout(title=title, xaxis_title="Azimuth from the tongue [deg]",
                         yaxis_title=y_title, template="plotly_white")
    figure.update_xaxes(dtick=60, gridcolor="#e8e8e6")
    figure.update_yaxes(gridcolor="#e8e8e6")
    return figure


def meridional_section(case, station=-1):
    """One section placed in the meridional plane, as ``(radius, axial)`` [mm].

    Drawing in the section's own coordinates would hide the tilt, since tilt is
    a rotation of exactly that frame. Mapping into the plane through the
    machine axis puts each passage where it actually sits.
    """
    volute = case.volute
    ring = volute.section.polygon(volute.geometry.size[station], n_points=240)[0]
    u, v = ring[:, 0], ring[:, 1]
    tilt = volute.tilt
    radius = volute.inlet.radius + u*math.cos(tilt) - v*math.sin(tilt)
    axial = u*math.sin(tilt) + v*math.cos(tilt)
    return radius*1e3, axial*1e3


def section_figure(cases):
    """Collection-plane sections where they sit in the meridional plane."""
    figure = go.Figure()
    for case, color in zip(cases, SERIES, strict=True):
        radius, axial = meridional_section(case)
        figure.add_scatter(x=radius, y=axial, name=case.name, mode="lines",
                           line=dict(color=color, width=2))
        # The slot mouth: where the flow enters, and what the tilt turns.
        half = case.volute.section.inlet_width / 2.0
        tilt = case.volute.tilt
        figure.add_scatter(
            x=[(case.volute.inlet.radius + half*math.sin(tilt))*1e3,
               (case.volute.inlet.radius - half*math.sin(tilt))*1e3],
            y=[-half*math.cos(tilt)*1e3, half*math.cos(tilt)*1e3],
            mode="lines", line=dict(color=color, width=4), showlegend=False,
            hovertemplate=f"{case.name} inlet slot<extra></extra>",
        )
    figure.update_layout(
        title="Collection-plane sections in the meridional plane, to scale",
        xaxis_title="Radius from the machine axis [mm]",
        yaxis_title="Axial position [mm]",
        yaxis_scaleanchor="x", template="plotly_white",
    )
    figure.update_xaxes(gridcolor="#e8e8e6")
    figure.update_yaxes(gridcolor="#e8e8e6", zerolinecolor="#d0d0cd")
    return figure


def tongue_section_figure(case):
    """Sections through the flattened run, against the circle they replace."""
    geometry = case.volute.geometry
    section = case.volute.section
    figure = go.Figure()
    # Sample evenly in aspect rather than in azimuth: the flattened run is
    # short and most of it sits pinned at the limit, so picking by angle would
    # draw the same clipped section several times over.
    shades = ("#c7dcf5", "#8fbcea", "#5a9ade", "#2a78d6")
    wanted = np.linspace(geometry.aspect.max(), geometry.aspect.min(), len(shades))
    for target, color in zip(wanted, shades, strict=True):
        index = int(np.argmin(np.abs(geometry.aspect - target)))
        ring = section.polygon(geometry.size[index], n_points=240,
                               aspect=geometry.aspect[index])[0]
        figure.add_scatter(
            x=ring[:, 0]*1e3, y=ring[:, 1]*1e3, mode="lines",
            name=f"{math.degrees(geometry.theta[index]):.0f}° · aspect {geometry.aspect[index]:.2f}",
            line=dict(color=color, width=2),
        )
    # What a fixed circular bore would have had to use over the same run.
    circle = section.polygon(section.minimum_size, n_points=240, aspect=1.0)[0]
    figure.add_scatter(x=circle[:, 0]*1e3, y=circle[:, 1]*1e3, mode="lines",
                       name="smallest circle", line=dict(color="#eb6834", width=2, dash="dash"))
    figure.update_layout(
        title=f"{case.name}: sections through the flattened run near the tongue",
        xaxis_title="Outward from the inlet plane [mm]",
        yaxis_title="Across the inlet slot [mm]",
        yaxis_scaleanchor="x", template="plotly_white",
    )
    figure.update_xaxes(gridcolor="#e8e8e6", zerolinecolor="#d0d0cd")
    figure.update_yaxes(gridcolor="#e8e8e6", zerolinecolor="#d0d0cd")
    return figure


def pressure_figure(cases):
    """Stagnation and static pressure round the passage and out the duct.

    Plotted as a change from each volute's own inlet, since the three sit at
    different supply pressures and it is the drop that is being compared.
    """
    figure = go.Figure()
    for case, color in zip(cases, SERIES, strict=True):
        result = case.result
        supply = result.states[0].p
        azimuth = np.degrees(result.geometry.theta)
        duct_stagnation, duct_static = result.exit_cone_pressures(13)
        # The duct has no azimuth of its own; lay it out past the revolution so
        # the whole path reads left to right.
        duct_azimuth = 360.0 + np.linspace(0.0, 60.0, len(duct_stagnation))
        figure.add_scatter(
            x=np.concatenate((azimuth, duct_azimuth)),
            y=(np.concatenate((result.stagnation_pressure(), duct_stagnation)) - supply)/1e5,
            name=f"{case.name} — stagnation", mode="lines",
            line=dict(color=color, width=2),
        )
        figure.add_scatter(
            x=np.concatenate((azimuth, duct_azimuth)),
            y=(np.concatenate((result.static_pressure(), duct_static)) - supply)/1e5,
            name=f"{case.name} — static", mode="lines",
            line=dict(color=color, width=2, dash="dot"),
        )
    figure.add_vline(x=360.0, line=dict(color="#d0d0cd", width=1, dash="dash"))
    figure.update_layout(
        title="Pressure along the passage and out through the exit cone",
        xaxis_title="Azimuth from the tongue [deg] — past 360°, along the duct",
        yaxis_title="Change from the inlet [bar]", template="plotly_white",
    )
    figure.update_xaxes(dtick=60, gridcolor="#e8e8e6")
    figure.update_yaxes(gridcolor="#e8e8e6", zerolinecolor="#d0d0cd")
    return figure


def design_rows(cases):
    """Design inputs and resolved geometry, one column per volute."""
    rows = [
        ("Propellant", lambda c: c.fluid),
        ("Inlet temperature [K]", lambda c: f"{c.T:.0f}"),
        ("Inlet pressure [bar]", lambda c: f"{c.p/1e5:.0f}"),
        ("Mass flow [kg/s]", lambda c: f"{c.mass_flow:.1f}"),
        ("Density [kg/m³]", lambda c: f"{c.result.density:.1f}"),
        ("Inlet radius r₅ [mm]", lambda c: f"{c.volute.inlet.radius*1e3:.1f}"),
        ("Inlet slot width b₅ [mm]", lambda c: f"{c.volute.section.inlet_width*1e3:.1f}"),
        ("Fillet radius [mm]", lambda c: f"{c.volute.section.fillet_radius*1e3:.1f}"),
        ("Design aspect ratio [-]", lambda c: f"{c.volute.section.aspect_ratio:.2f}"),
        ("Aspect ratio limit [-]", lambda c: f"{c.volute.section.max_aspect_ratio:.2f}"),
        ("Inlet plane tilt [deg]", lambda c: f"{math.degrees(c.volute.tilt):.0f}"),
        ("Meridional velocity C_m5 [m/s]", lambda c: f"{c.volute.inlet.C_m:.2f}"),
        ("Swirl velocity C_u5 [m/s]", lambda c: f"{c.volute.inlet.C_u:.1f}"),
        ("Inlet flow angle α₅ [deg]", lambda c: f"{math.degrees(c.volute.inlet.alpha):.1f}"),
        ("Bore radius at collection plane [mm]", lambda c: f"{c.volute.geometry.size[-1]*1e3:.1f}"),
        ("Collection-plane area A₆ [cm²]", lambda c: f"{c.result.A_6*1e4:.2f}"),
        ("Collection-plane centroid r₆ [mm]", lambda c: f"{c.result.r_6*1e3:.1f}"),
        ("Equivalent diameter d_H [mm]", lambda c: f"{c.result.d_H*1e3:.1f}"),
        ("Largest radius reached [mm]", lambda c: f"{c.outer_radius*1e3:.1f}"),
        ("Flattened out to [deg]", lambda c: (
            f"{math.degrees(c.volute.geometry.theta[c.volute.geometry.flattened].max()):.0f}"
            if c.volute.geometry.flattened.any() else "-")),
        ("Area law clipped to [deg]", lambda c: f"{math.degrees(c.volute.geometry.tongue_clip_angle):.0f}"),
    ]
    return rows


def outlet_rows(cases):
    """The exit-cone duct: its shape inputs and what they resolved to."""
    rows = [
        ("First angle input", lambda c: (
            c.volute.exit_cone.first_angle
            if isinstance(c.volute.exit_cone.first_angle, str)
            else f"{math.degrees(c.volute.exit_cone.first_angle):.0f}°"
        )),
        ("Resolved first angle [deg]", lambda c: f"{math.degrees(c.volute.exit_geometry.first_angle):.1f}"),
        ("Radial bound at this arc length [deg]", lambda c: f"{math.degrees(c.volute.exit_geometry.radial_angle):.1f}"),
        ("Second angle (elevation) [deg]", lambda c: f"{math.degrees(c.volute.exit_geometry.second_angle):+.0f}"),
        ("Arc length parameter [outlet diameters]", lambda c: f"{c.volute.exit_cone.arc_length_parameter:.1f}"),
        ("Arc length [mm]", lambda c: f"{c.volute.exit_geometry.arc_length*1e3:.1f}"),
        ("Duct turn angle [deg]", lambda c: f"{math.degrees(c.volute.exit_geometry.turn_angle):.1f}"),
        ("Centreline bend radius [mm]", lambda c: (
            "straight" if c.volute.exit_geometry.is_straight
            else f"{c.volute.exit_geometry.bend_radius*1e3:.1f}"
        )),
        ("Flange diameter [mm]", lambda c: f"{c.exit_diameter*1e3:.1f}"),
        ("Flange azimuth from the volute outlet [deg]", lambda c: f"{math.degrees(c.volute.exit_geometry.end_azimuth):.1f}"),
    ]
    return rows


def performance_rows(cases):
    rows = [
        ("Sizing parameter SP [-]", lambda c: f"{c.result.SP:.4f}"),
        ("Inlet velocity C₅ [m/s]", lambda c: f"{c.result.C_5:.1f}"),
        ("Collection-plane velocity C₆ [m/s]", lambda c: f"{c.result.C_6:.1f}"),
        ("Exit cone velocity C₇ [m/s]", lambda c: f"{c.result.C_7:.1f}"),
        ("Reynolds number [-]", lambda c: f"{c.result.reynolds:.3g}"),
        ("Fanning friction factor c_f [-]", lambda c: f"{c.result.friction_factor:.4f}"),
        ("Mean path length [mm]", lambda c: f"{c.result.path_length*1e3:.1f}"),
    ]
    rows += [(f"ω {label.lower()} [-]",
              (lambda key: lambda c: f"{c.result.loss_breakdown()[key]:.4f}")(key))
             for key, label in MECHANISM_LABELS.items()]
    rows += [
        ("ω total [-]", lambda c: f"{c.result.omega_total:.4f}"),
        ("Total pressure loss [bar]", lambda c: f"{c.result.pressure_loss/1e5:.3f}"),
    ]
    return rows


def markdown_table(cases, rows):
    header = "| | " + " | ".join(case.name for case in cases) + " |"
    rule = "|---|" + "---:|"*len(cases)
    body = [f"| {label} | " + " | ".join(value(case) for case in cases) + " |"
            for label, value in rows]
    return "\n".join([header, rule, *body])


def create_report(output_dir, *, screenshot=False):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    cases = build_cases()
    scenes = [build_scene(case, color)
              for case, color in zip(cases, SERIES, strict=True)]
    for case, scene in zip(cases, scenes, strict=True):
        scene.save_html(output_dir / f"volute-scene-{slug(case.name)}.html")

    report = Report("Volute design study")

    overview = report.add_tab("Volutes")
    overview.add_markdown("""
        # Three volutes

        Each passage gets its own scene below. Rotate and zoom to inspect them.
        Every one is swept from its tongue through a full revolution and out
        along its outlet duct to the flange, so the step visible where the
        passage meets itself is the cutwater. Colour shows the through-flow
        velocity that passage carries at each azimuth; use the field selector
        in the legend to colour by passage area, bore radius, centroid radius,
        or the fraction of mass flow collected so far.

        They are kept in separate scenes rather than lined up together so each
        gets its own colour scale: the manifold runs at about a third of the
        LOX volute's velocity, and a shared scale flattens it to one dark tone.

        The two pump-discharge volutes are radial. The third is a cooling-jacket
        inlet manifold with its inlet plane inclined 45 degrees, which is what
        lets a manifold sit against a sloping nozzle wall. A bore rotated about
        its own axis looks much the same from outside, so the tilt is hard to
        read here — the meridional section on the **Design** tab shows it
        directly. Tilt is placement only: it moves the passage in space and
        changes nothing in the loss analysis.

        The outlet ducts differ too, and those *are* easy to read here:
    """)
    overview.add_markdown(
        "| Case | Outlet |\n|---|---|\n"
        + "\n".join(f"| {c.name} | {c.outlet_note} |" for c in cases)
    )
    for case, scene in zip(cases, scenes, strict=True):
        overview.add_markdown(f"### {case.name}\n\n{case.note}")
        overview.add_iframe(
            scene.data_url,
            caption=f"{case.name} — {case.outlet_note} Coloured by through-flow velocity.",
            height="clamp(320px, 62vh, 640px)",
        )
    overview.add_markdown(f"""
        All three are sized on conservation of angular momentum at a sizing
        parameter of 1, and all three realise it exactly:
        {", ".join(f"**{c.name}** SP = {c.result.SP:.4f}" for c in cases)}.
    """)

    design_tab = report.add_tab("Design")
    design_tab.add_markdown("## Inputs and resolved geometry\n\n"
                            + markdown_table(cases, design_rows(cases)))
    design_tab.add_markdown("""
        The meridional velocity is not an independent choice. The mass flow has
        to cross the inlet slot, so `C_m5` follows from the duty, the inlet
        radius and the slot width that the feeding hardware presents. Only the
        swirl is chosen. The inlet flow angle that results is what sets how fast
        the spiral has to grow: the passage area over its centroid radius rises
        at exactly `b₅·tan α₅` per radian, so a steep inlet needs a fat volute.

        That is why the manifold is the largest of the three despite carrying a
        third of the LOX volute's mass flow. It sits at a large radius with
        modest swirl, which makes its inlet the steepest.

        Near the tongue a shaped volute wants a passage area approaching zero,
        which no bore spanning the slot can provide: it can never be narrower
        than the slot feeding it. Its *shape* is not stuck, though. Squashing
        the bore radially while it still spans the slot shrinks the area
        roughly as one over the aspect ratio, so the passage keeps tracking the
        area law well past where a fixed circular bore would have had to give
        up. The two rows above show where each volute starts flattening and
        where, having run out of aspect ratio, it finally clips.
    """)
    design_tab.add_figure(tongue_section_figure(cases[0]))
    design_tab.add_markdown("""
        Each section is drawn where it sits in the flattened run, against the
        smallest circle that would otherwise have had to cover the whole of it.
        Both span the slot exactly — they have to — but the flattened ones
        enclose far less area, which is what buys the extra reach. The bore
        also stops reaching as far outward as it flattens, so the passage
        tucks closer to the inlet plane near the tongue as a side effect.
    """)
    design_tab.add_figure(section_figure(cases))
    design_tab.add_markdown("""
        Drawn where they actually sit relative to the machine axis, with the
        thick segment marking each inlet slot. The manifold's slot is the one
        leaning over at 45 degrees: the passage it feeds is carried up and out
        rather than straight out. Tilting the inlet plane rotates the section
        frame rigidly, so it changes packaging and nothing else — every area,
        centroid radius and loss coefficient in this report is identical to what
        the same volute would give sitting flat.

        This section stops at the full-collection plane; the outlet duct beyond
        it is described next, and shown in the scene above. Only a purely
        radial duct stays inside this same meridional plane the whole way out,
        so a general outlet is not something a single 2-D section can show.
    """)

    design_tab.add_markdown("## Outlet ducts\n\n" + markdown_table(cases, outlet_rows(cases)))
    design_tab.add_markdown("""
        Two angles fix the outlet direction. The **first angle** turns within
        the volute's own plane, measured from *tangential* — the direction the
        flow is already travelling as it leaves the full-collection plane. That
        reference point is fixed at the volute outlet, because by the time the
        duct reaches its flange it has moved round the volute and there is no
        longer a tangential direction to be measured against there. *Radial* is
        different: it always means straight out from the machine axis,
        evaluated at wherever the flange actually ends up. Having travelled
        forward by some azimuth, that direction is only "radial bound" away
        from tangential rather than a full right angle, which is why the LOX
        duct's resolved first angle sits well under 90°: asking for `"radial"`
        solves for exactly that bound.

        The **second angle** lifts the outlet out of the volute's own plane —
        the manifold's alone is nonzero, which is what carries its duct up and
        away from the machine as well as round it.

        The methane pump's duct is straight because `"tangential"` asks for no
        turn at all: a duct with nothing to turn through has no bend radius to
        speak of. The other two curve, at a bend radius set by how far a given
        turn is spread over the requested arc length — tighter than the duct's
        own radius would fold it through itself, and the resolver rejects that
        rather than build unrealisable geometry.
    """)

    geometry_tab = report.add_tab("Geometry")
    geometry_tab.add_figure(_azimuth_figure(
        cases, lambda c: c.volute.geometry.area*1e4,
        "Passage area around the volute", "Passage area [cm²]"))
    geometry_tab.add_markdown("""
        Area grows from the tongue to the full-collection plane as more of the
        flow is gathered. The flat start on each curve is the clipped region
        where the section family has bottomed out.
    """)
    geometry_tab.add_figure(_azimuth_figure(
        cases, lambda c: c.volute.geometry.aspect,
        "Bore aspect ratio around the volute", "Aspect ratio [-]"))
    geometry_tab.add_markdown("""
        Each volute holds its design shape for most of the revolution and only
        flattens near the tongue, where the area law wants less area than a
        bore spanning the slot can enclose. Flattening is the second lever: the
        bore cannot narrow below the slot, but it can be squashed radially,
        which shrinks the area roughly as one over the aspect ratio. The flat
        top of each curve is where even the flattest allowed section is still
        too big and the law is finally clipped.
    """)
    geometry_tab.add_figure(_azimuth_figure(
        cases, lambda c: c.result.collected_velocity(),
        "Through-flow velocity around the volute", "Velocity [m/s]"))
    geometry_tab.add_markdown("""
        Only the fraction of mass flow collected so far passes any given
        station, so this is what the local area actually has to carry. It rises
        from zero through the clipped region, then **falls** steadily: under
        conservation of angular momentum the product of velocity and centroid
        radius is what stays constant, so the flow slows as the passage carries
        it outward. A constant-velocity sizing law would flatten these curves
        instead, at the cost of no longer conserving angular momentum.
    """)
    geometry_tab.add_figure(_azimuth_figure(
        cases, lambda c: (c.result.collected_velocity()
                          * c.volute.geometry.centroid_radius),
        "Angular momentum around the volute", "Velocity × centroid radius [m²/s]"))
    geometry_tab.add_markdown("""
        The same data as the curve above, multiplied by the centroid radius.
        Flat past the tongue is the sizing law being satisfied.
    """)

    loss_tab = report.add_tab("Losses")
    loss_tab.add_figure(loss_figure(cases))
    loss_tab.add_figure(pressure_figure(cases))
    loss_tab.add_markdown("""
        Stagnation pressure only ever falls; static pressure is the part of it
        left over after the local velocity head. It starts at the tongue equal
        to the stagnation pressure, because nothing has been collected there
        and nothing is moving, dips where the passage runs fastest just past
        the tongue, and climbs back as the flow slows on its way out.

        Past 360° the curves follow the exit-cone duct, where static pressure
        recovers again wherever the duct diffuses. Read the volute half with
        its model in mind: Aungier charges one lumped coefficient per mechanism
        across the whole passage and says nothing about where round the volute
        each is spent, so spreading them linearly with azimuth is an assumption
        made so the passage can be plotted and coloured at all. The endpoints
        are what the analysis actually fixes. The duct half is on firmer
        ground, since its losses genuinely are spent between two stations the
        model resolves.
    """)
    loss_tab.add_markdown("## Performance\n\n"
                          + markdown_table(cases, performance_rows(cases)))
    loss_tab.add_markdown("""
        Every coefficient is referred to the inlet velocity head, so the total
        pressure loss is `ω_total · ρ · C₅² / 2`.

        The tangential dump term is zero in all three. That is by construction
        rather than by luck: sizing on conservation of angular momentum at
        SP = 1 is exactly the condition that zeroes it, and the same sizing
        parameter appears in both the design equation and the loss model. It
        becomes nonzero as soon as a volute runs away from the duty it was
        sized for, or is deliberately oversized.

        What remains is largely outside the volute's control. The meridional
        dump charges the whole through-slot velocity head, because that
        component turns into swirl inside the passage and is never recovered.
        It is fixed by `C_m5`, which the pump and diffuser upstream decide. The
        only real lever the volute designer has over it is to ask for a wider
        inlet slot.

        The three cases land in different places for that reason. The LOX
        volute has a shallow inlet and a long fast path, so friction leads. The
        methane pump volute has a steeper inlet at lower density, so the
        meridional dump dominates. The manifold has the steepest inlet of all,
        and also pays an exit-cone term because its passage decelerates the flow
        well below the speed the exit pipe carries.

        The exit cone is charged three ways. **Exit cone diffusion** is the
        Rennels conical loss, set by the angle the cone tapers at — so the
        flange diameter and the arc length decide it together, not the area
        ratio alone. The manifold's flange is narrower than its passage, so it
        contracts rather than diffuses, and is charged the contraction
        correlation instead: much the cheaper of the two, since nothing
        separates on the way in. **Duct friction** and **duct bend** come from
        the same resolved geometry on the Design tab, and are what make a
        radial outlet cost more than a tangential one at the same flange
        diameter.

        That leaves the arc length parameter buying two trades at once.
        Lengthening the duct shallows the cone and recovers more of the
        velocity head, at the cost of more wetted wall; it also spreads a given
        turn over more duct, trading the bend term for the friction term. Both
        have an optimum in them, and neither is visible if the diffusion term
        is written as a function of area ratio alone.
    """)

    report_path = output_dir / "volute-report.html"
    report.save_html(str(report_path))
    if screenshot:
        for case, scene in zip(cases, scenes, strict=True):
            plotter = scene.to_plotter(off_screen=True)
            try:
                plotter.window_size = (1100, 900)
                # Looking down and in from outside the spiral: high enough to
                # read the passage growing, low enough to show the outlet duct
                # leaving and, for the manifold, lifting out of the plane.
                reach = case.outer_radius
                plotter.camera_position = [
                    (2.1*reach, -2.4*reach, 1.5*reach), (0.0, 0.0, 0.0), (0, 0, 1)]
                plotter.screenshot(
                    str(output_dir / f"volute-preview-{slug(case.name)}.png"))
            finally:
                plotter.close()
    for case in cases:
        print(f"{case.name:<24} SP={case.result.SP:.4f}  "
              f"A6={case.result.A_6*1e4:6.2f} cm2  "
              f"omega={case.result.omega_total:.4f}  "
              f"dp={case.result.pressure_loss/1e5:.3f} bar")
    print(f"Report: {report_path.resolve()}")
    return cases, scenes, report_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent / "output")
    parser.add_argument("--screenshot", action="store_true", help="Also render a PyVista PNG preview")
    args = parser.parse_args()
    create_report(args.output_dir, screenshot=args.screenshot)
