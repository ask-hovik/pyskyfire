selector_to_html = {"a[href=\"#pyskyfire-viz-plot-film-cooling-plotfilmboundarylayer\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.plot_film_cooling.PlotFilmBoundaryLayer<a class=\"headerlink\" href=\"#pyskyfire-viz-plot-film-cooling-plotfilmboundarylayer\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.viz.plot_film_cooling.PlotFilmBoundaryLayer\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.plot_film_cooling.PlotFilmBoundaryLayer\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.plot_film_cooling.</span></span><span class=\"sig-name descname\"><span class=\"pre\">PlotFilmBoundaryLayer</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">gaseous</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">template</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/stdtypes.html#str\" title=\"(in Python v3.15)\"><span class=\"pre\">str</span></a></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">'plotly_white'</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Bases: <a class=\"reference internal\" href=\"../core/PlotBase.html#pyskyfire.viz.core.PlotBase\" title=\"pyskyfire.viz.core.PlotBase\"><code class=\"xref py py-obj docutils literal notranslate\"><span class=\"pre\">pyskyfire.viz.core.PlotBase</span></code></a></p><p>Boundary-layer mass flux in the gaseous film phase.</p></dd>", "a[href=\"../core/PlotBase.html#pyskyfire.viz.core.PlotBase\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.core.PlotBase\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.core.</span></span><span class=\"sig-name descname\"><span class=\"pre\">PlotBase</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">fig</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><span class=\"pre\">plotly.graph_objects.Figure</span><span class=\"w\"> </span><span class=\"p\"><span class=\"pre\">|</span></span><span class=\"w\"> </span><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/constants.html#None\" title=\"(in Python v3.15)\"><span class=\"pre\">None</span></a></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">None</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Common helpers + Plotly passthrough.</p></dd>"}
skip_classes = ["headerlink", "sd-stretched-link"]

window.onload = function () {
    for (const [select, tip_html] of Object.entries(selector_to_html)) {
        const links = document.querySelectorAll(` ${select}`);
        for (const link of links) {
            if (skip_classes.some(c => link.classList.contains(c))) {
                continue;
            }

            tippy(link, {
                content: tip_html,
                allowHTML: true,
                arrow: true,
                placement: 'auto-start', maxWidth: 600, interactive: false,
                onShow(instance) {MathJax.typesetPromise([instance.popper]).then(() => {});},
            });
        };
    };
    console.log("tippy tips loaded!");
};
