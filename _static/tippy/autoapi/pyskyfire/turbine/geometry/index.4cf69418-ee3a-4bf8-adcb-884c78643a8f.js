selector_to_html = {"a[href=\"TurbineGeometry.html#pyskyfire.turbine.geometry.TurbineGeometry\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.turbine.geometry.TurbineGeometry\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.turbine.geometry.</span></span><span class=\"sig-name descname\"><span class=\"pre\">TurbineGeometry</span></span></dt><dd><p>Every row\u2019s built blade, and the passage they stand in.</p></dd>", "a[href=\"#classes\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Classes<a class=\"headerlink\" href=\"#classes\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"#module-pyskyfire.turbine.geometry\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.turbine.geometry<a class=\"headerlink\" href=\"#module-pyskyfire.turbine.geometry\" title=\"Link to this heading\">\u00b6</a></h1><p>Built turbine blades, and what can only be measured on them.</p><p>A <a class=\"reference internal\" href=\"TurbineGeometry.html#pyskyfire.turbine.geometry.TurbineGeometry\" title=\"pyskyfire.turbine.geometry.TurbineGeometry\"><code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">TurbineGeometry</span></code></a> is what <code class=\"xref py py-meth docutils literal notranslate\"><span class=\"pre\">TurbineDefinition.build_geometry()</span></code>\nreturns: one <a class=\"reference internal\" href=\"RowGeometry.html#pyskyfire.turbine.geometry.RowGeometry\" title=\"pyskyfire.turbine.geometry.RowGeometry\"><code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">RowGeometry</span></code></a> per row, and the passage they sit in. Each\nrow carries one blade\u2019s surfaces, which the scene repeats round the wheel,\nand the numbers the drawing left to the built section:</p>", "a[href=\"RowGeometry.html#pyskyfire.turbine.geometry.RowGeometry\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.turbine.geometry.RowGeometry\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.turbine.geometry.</span></span><span class=\"sig-name descname\"><span class=\"pre\">RowGeometry</span></span></dt><dd><p>One row\u2019s blade and its measured section properties.</p></dd>"}
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
