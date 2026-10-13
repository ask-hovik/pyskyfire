selector_to_html = {"a[href=\"#rao-nozzle-contour-angles\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Rao Nozzle Contour Angles<a class=\"headerlink\" href=\"#rao-nozzle-contour-angles\" title=\"Link to this heading\">\u00b6</a></h1><p>A Rao-style bell nozzle replaces a longer conical nozzle with a curved contour.\nThe contour begins at the end of the throat blend with wall angle <span class=\"math notranslate nohighlight\">\\(\\theta_n\\)</span>\nand reaches the nozzle exit with wall angle <span class=\"math notranslate nohighlight\">\\(\\theta_e\\)</span>. Both angles depend on\nthe nozzle expansion ratio and the selected nozzle length.</p>", "a[href=\"#source\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Source<a class=\"headerlink\" href=\"#source\" title=\"Link to this heading\">\u00b6</a></h2><p>The angle curves were digitised from page 15 of <em>Liquid Rocket Engine Nozzles</em>,\nNASA Space Vehicle Design Criteria (Chemical Propulsion), NASA SP-8120, July\n1976. The original report is available from the <a class=\"reference external\" href=\"https://ntrs.nasa.gov/citations/19770009165\">NASA Technical Reports\nServer</a>.</p>", "a[href=\"#how-the-chart-is-generated\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">How the chart is generated<a class=\"headerlink\" href=\"#how-the-chart-is-generated\" title=\"Link to this heading\">\u00b6</a></h2><p><code class=\"docutils literal notranslate\"><span class=\"pre\">tools/generate_engineering_charts.py</span></code> creates the chart with\n<code class=\"docutils literal notranslate\"><span class=\"pre\">pyskyfire.viz.PlotThetaVsEpsilon</span></code>. The plotting class reads the digitised\n<span class=\"math notranslate nohighlight\">\\(\\theta_n\\)</span> and <span class=\"math notranslate nohighlight\">\\(\\theta_e\\)</span> curves in\n<code class=\"docutils literal notranslate\"><span class=\"pre\">src/pyskyfire/regen/data/theta_n.json</span></code> and\n<code class=\"docutils literal notranslate\"><span class=\"pre\">src/pyskyfire/regen/data/theta_e.json</span></code>, then plots the data for nozzle length\nfractions from 60% to 100%. These are the same data tables used by Pyskyfire\u2019s\nRao contour generator.</p>"}
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
