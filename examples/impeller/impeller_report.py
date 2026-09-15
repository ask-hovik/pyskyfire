"""Size an impeller, generate pyTurbo surfaces, and embed a scene in a report.

Run from the repository with:
    python examples/impeller/impeller_report.py --output-dir examples/impeller/output
Requires pyskyfire[impeller,viz]. No browser or display is needed for HTML output.
"""
import argparse
from pathlib import Path

import plotly.graph_objects as go

from pyskyfire.pump import BladeSettings, Impeller, SplitterSettings
from pyskyfire.viz import Scene, ImpellerComponent, Report


def thickness(*, chord_fraction, span_fraction):
    """A smooth physical body thickness, with a thicker hub and thinner exit."""
    return 0.0012*(1-0.30*span_fraction)*(1-0.20*chord_fraction)


def build_impeller():
    impeller = Impeller(
        # A low specific speed duty: high head from a small, narrow wheel.
        Q=0.008, H=900., n=25_000,
        shaft_diameter=0.012, inlet_hub_diameter=0.020,
        phi1=0.20, blade_count=4, beta2_deg=30.,
        blade_thickness=thickness, incidence_deg=2.,
        volumetric_efficiency=0.97, hydraulic_efficiency=0.85,
        mechanical_efficiency=0.98, rho=1000.,
        # A slanted leading edge: the blade starts well down the hub and up
        # in the eye at the tip, so both wrap by a comparable amount.
        blade_settings=BladeSettings(start=(0.2, 0.01), end=1.0),
        # Four full blades keep the eye open; the splitters double the blade
        # count over the outer half, where the passage is too wide to guide.
        splitter=SplitterSettings(start=0.45),
    )
    impeller.build_geometry(span_points=33, chord_points=129)
    return impeller


def create_report(output_dir, *, screenshot=False):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    impeller = build_impeller()
    scene = Scene(title="Centrifugal impeller", show_axes=True)
    scene.add(ImpellerComponent(impeller), name="impeller")
    scene.save_html(output_dir / "impeller-scene.html")
    impeller.geometry.export_npz(output_dir / "impeller-geometry.npz")

    report = Report("Centrifugal impeller design")
    tab = report.add_tab("Impeller")
    tab.add_markdown("""
        # Centrifugal impeller
        Rotate and zoom to inspect the blades. The translucent surfaces show
        the hub and shroud flow boundaries. See Design for dimensions and assumptions.
    """)
    tab.add_iframe(scene.data_url, caption="Impeller generated with pyTurbo Aero",
                   height="clamp(360px, calc(100vh - 260px), 720px)")

    d = impeller.dimensions
    design = report.add_tab("Design")
    design.add_table({
        "Delivered flow [m³/s]": impeller.Q, "Head [m]": impeller.H,
        "Shaft speed [rpm]": impeller.n, "Shaft diameter [m]": d.shaft_diameter,
        "Inlet hub diameter [m]": d.inlet_hub_diameter, "Eye diameter [m]": d.d1,
        "Outlet diameter [m]": d.d2, "Outlet width [m]": d.b2,
        "Specific speed nq": d.nq, "Full blade count": d.blade_count,
        "Splitter start [main chord fraction]": impeller.blade_definition.splitter.start,
        "Internal flow [m³/s]": d.internal_flow,
        "Volumetric efficiency (supplied)": impeller.volumetric_efficiency,
        "Hydraulic efficiency (supplied)": impeller.hydraulic_efficiency,
        "Mechanical efficiency (supplied)": impeller.mechanical_efficiency,
        "Shaft power from duty and supplied efficiencies [W]": impeller.shaft_power,
    })
    design.add_markdown("""
        This duty sits at low specific speed: a small, narrow wheel turning
        fast for a lot of head. Two things follow from that, and they pull the
        blade count in opposite directions. The outlet is narrow relative to
        its diameter, so a passage bounded by only a few blades is far wider
        circumferentially than it is tall and guides the flow poorly. But the
        eye is small and its relative flow angle is shallow, so every blade
        that reaches the leading edge blocks a disproportionate share of the
        inlet throat.

        Splitters resolve the conflict. Four full blades set the inlet
        blockage, while a splitter in each passage starts downstream of the
        eye and doubles the blade count over the outer half of the wheel,
        where the pitch has grown and the guidance is needed. Splitters follow
        the main-blade camber, so they add no independent angle inputs.

        Sizing uses the preliminary head-coefficient and width correlations.
        Inlet angles come from unblocked continuity through the leading edge
        with zero prewhirl, plus 2° incidence. The outlet metal angle is a
        design input. These calculations do not include a slip or loss
        prediction. The efficiencies are supplied assumptions. Wrap is
        estimated separately at each span; the slanted leading edge keeps the
        hub and tip values close together.
        The thickness function specifies body thickness; the blade edges close.
    """)
    figure = go.Figure()
    for name, curve in (("Hub", impeller.geometry.hub_curve), ("Shroud", impeller.geometry.shroud_curve)):
        figure.add_scatter(x=curve[:,0], y=curve[:,1], name=name, mode="lines")
    figure.update_layout(title="Meridional passage", xaxis_title="Axial x [m]",
                         yaxis_title="Radius [m]", yaxis_scaleanchor="x", template="plotly_white")
    design.add_figure(figure)
    design.add_markdown("\n".join([
        "| Span | Inlet metal angle [deg] | Outlet metal angle [deg] | Wrap [deg] |",
        "|---|---:|---:|---:|",
        *(f"| {s.span_fraction:g} | {s.beta1_deg:.2f} | {s.beta2_deg:.2f} | {s.wrap_deg:.2f} |"
          for s in impeller.blade_definition.sections),
    ]))
    report_path = output_dir / "impeller-report.html"
    report.save_html(str(report_path))
    if screenshot:
        plotter = scene.to_plotter(off_screen=True)
        try:
            plotter.window_size = (1400, 1000)
            plotter.screenshot(str(output_dir / "impeller-preview.png"))
        finally:
            plotter.close()
    print(impeller)
    print(f"Body thickness maximum relative error: {impeller.geometry.thickness_max_relative_error:.2%}")
    print(f"Report: {report_path.resolve()}")
    return impeller, scene, report_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent / "output")
    parser.add_argument("--screenshot", action="store_true", help="Also render a PyVista PNG preview")
    args = parser.parse_args()
    create_report(args.output_dir, screenshot=args.screenshot)
