selector_to_html = {"a[href=\"#module-pyskyfire.viz.backends.pyvista\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.backends.pyvista<a class=\"headerlink\" href=\"#module-pyskyfire.viz.backends.pyvista\" title=\"Link to this heading\">\u00b6</a></h1><p>Desktop rendering of a scene, retaining named scalar arrays and node IDs.</p>", "a[href=\"#functions\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Functions<a class=\"headerlink\" href=\"#functions\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"#pyskyfire.viz.backends.pyvista.make_plotter\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.backends.pyvista.make_plotter\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.backends.pyvista.</span></span><span class=\"sig-name descname\"><span class=\"pre\">make_plotter</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">scene</span></span></em>, <em class=\"sig-param\"><span class=\"keyword-only-separator o\"><abbr title=\"Keyword-only parameters separator (PEP 3102)\"><span class=\"pre\">*</span></abbr></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">off_screen</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">False</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#module-contents\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Module Contents<a class=\"headerlink\" href=\"#module-contents\" title=\"Link to this heading\">\u00b6</a></h2>"}
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
