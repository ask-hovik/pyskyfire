"""Responsive sizing checks for generated HTML reports."""

import plotly.graph_objects as go

from pyskyfire.viz.report import Report


def test_report_presizes_all_plots_before_activating_a_tab(tmp_path) -> None:
    report = Report("Responsive report")
    report.add_tab("Plots").add_figure(go.Figure(go.Scatter(y=[1, 2])))
    output_path = tmp_path / "report.html"

    report.save_html(output_path)

    document = output_path.read_text()
    assert "async function psfResizePlots(panel)" in document
    assert "requestAnimationFrame" in document
    assert ".tab-panel.psf-preparing {" in document
    assert "visibility: hidden;" in document
    assert "position: absolute;" in document
    assert "function psfAfterLayoutAndPaint()" in document
    assert "await psfAfterLayoutAndPaint();" in document
    assert "await Plotly.Plots.resize(plot);" in document
    assert "const plots = Array.from" in document
    assert "if (!plots.length) return;" in document
    assert "async function psfResizeAllPlots()" in document
    assert "await Promise.all(panels.map(panel => psfResizePlots(panel)));" in document
    assert "panel.classList.add('psf-preparing')" in document
    assert "panel.classList.remove('psf-preparing')" in document
    assert "await psfResizeAllPlots();" in document
    assert document.index("await psfResizeAllPlots();") < document.index(
        "psfActivateTab(idx);"
    )

    activate_tab = document[
        document.index("function psfActivateTab(idx)") : document.index(
            "window.addEventListener('DOMContentLoaded'"
        )
    ]
    assert "psfResizePlots" not in activate_tab
    assert "const PSF_PLOT_ASPECT_RATIO = 16 / 9;" in document
    assert "width / PSF_PLOT_ASPECT_RATIO" in document
    assert "Plotly.Plots.resize(plot)" in document
    assert "async function psfPlaceLegend(plot, width)" in document
    assert "const PSF_SIDE_LEGEND_MIN_WIDTH = 900;" in document
    assert "'legend.x': 1.02" in document
    assert "'legend.y': -0.18" in document
    assert "base.right + Math.ceil(bounds.width) + 24" in document
    assert "base.bottom + Math.ceil(bounds.height) + 24" in document
    assert "window.addEventListener('resize'" in document
    assert "clearTimeout(psfWindowResizeTimer);" in document
    assert "psfResizeAllPlots().catch(" in document
    assert "grid-template-columns: var(--sidebar-w) minmax(0, 1fr);" in document
    assert ".content {" in document
    assert "min-width: 0;" in document
    assert "width: max-content;" in document
    assert "max-width: 100%;" in document
    assert "overflow-wrap: anywhere;" in document
    assert "width: 34%;" not in document
