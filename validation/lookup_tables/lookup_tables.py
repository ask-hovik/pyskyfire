import pyskyfire as psf
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
report = psf.viz.Report("Lookup Tables")

# -------------- Theta vs Epsilon ----------------
tab1 = report.add_tab("Theta Vs Epsilon")

text = "The construction of a Rao nozzle depends on data for θn and θe. Those angles have been graphed on page 15 of C. Hyde, J. and S. Gill, G.˙ Liquid rocket engine nozzles. NASA Space Vehicle Design Criteria--Chemical Propulsion SP-8120, National Aeronautics and Space Administration (NASA), Washington, D.C., July 1976. Prepared under the direction of NASA Lewis Research Center."
image = os.path.join(script_dir, "thrust_chamber_dimensions.png")
theta_v_epsilon = psf.viz.PlotThetaVsEpsilon()

tab1.add_text(text)
tab1.add_image(image, caption="Thrust chamber geometrical definitions")
tab1.add_figure(theta_v_epsilon, caption="Variation of nozzle contour angles θn and θe against area ratio, ε. Percentages correspond to various length fractions of the bell nozzle compared to a full length 15◦ conical nozzle")

# -------------- Moody Diagram ----------------
tab2 = report.add_tab("Moody Diagram")

text = "The moody diagram shows friction factor as a function of Reynolds number. Of particular importance is the transition region between laminar and turbulent regimes. In Pyskyfire this is done by simple interpolation when 2300 < Re < 3500. \n \n There were some issues with the log functionality in my graphing code, so the plot is looking a bit weird. To be remedied in the future. "
moody_diagram = psf.viz.PlotMoodyDiagram()

tab2.add_text(text)
tab2.add_figure(moody_diagram, caption = "Moody diagram resulting from the linear interpolation approach. Laminar to turbulent transition is fixed at 2300 to 3500. The graph is plotted for different values of r/D.")

out_path = os.path.join(script_dir, "lookup_tables.html")
report.save_html(out_path)
print(f"Report saved to {out_path}")