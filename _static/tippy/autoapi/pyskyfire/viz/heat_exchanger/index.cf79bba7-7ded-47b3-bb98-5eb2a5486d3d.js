selector_to_html = {"a[href=\"#module-pyskyfire.viz.heat_exchanger\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.heat_exchanger<a class=\"headerlink\" href=\"#module-pyskyfire.viz.heat_exchanger\" title=\"Link to this heading\">\u00b6</a></h1><p>Heat-exchanger swept geometry and wall/flow report panels.</p>", "a[href=\"#module-contents\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Module Contents<a class=\"headerlink\" href=\"#module-contents\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"#pyskyfire.viz.heat_exchanger.heat_exchanger_report\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.heat_exchanger.heat_exchanger_report\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.heat_exchanger.</span></span><span class=\"sig-name descname\"><span class=\"pre\">heat_exchanger_report</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">exchanger</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">result</span></span></em>, <em class=\"sig-param\"><span class=\"keyword-only-separator o\"><abbr title=\"Keyword-only parameters separator (PEP 3102)\"><span class=\"pre\">*</span></abbr></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">title</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">'Heat</span> <span class=\"pre\">exchanger'</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">installation</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">None</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">scene</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">None</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Create a Report with geometry, fluid states, wall temperatures, and balance.</p></dd>", "a[href=\"#pyskyfire.viz.heat_exchanger.plot_heat_exchanger_geometry\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.heat_exchanger.plot_heat_exchanger_geometry\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.heat_exchanger.</span></span><span class=\"sig-name descname\"><span class=\"pre\">plot_heat_exchanger_geometry</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">exchanger</span></span></em>, <em class=\"sig-param\"><span class=\"keyword-only-separator o\"><abbr title=\"Keyword-only parameters separator (PEP 3102)\"><span class=\"pre\">*</span></abbr></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">installation</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">None</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Draw the actual exposed tube envelopes, including transported sections.</p></dd>", "a[href=\"#functions\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Functions<a class=\"headerlink\" href=\"#functions\" title=\"Link to this heading\">\u00b6</a></h2>"}
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
