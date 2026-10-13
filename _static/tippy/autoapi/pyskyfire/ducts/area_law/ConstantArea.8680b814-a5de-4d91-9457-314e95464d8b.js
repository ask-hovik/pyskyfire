selector_to_html = {"a[href=\"#pyskyfire.ducts.area_law.ConstantArea\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.area_law.ConstantArea\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.ducts.area_law.</span></span><span class=\"sig-name descname\"><span class=\"pre\">ConstantArea</span></span></dt><dd><p>Bases: <code class=\"xref py py-obj docutils literal notranslate\"><span class=\"pre\">AreaLaw</span></code></p><p>Hold the passage area uniform all the way around.</p><p>This is Aungier\u2019s \u201ccollector\u201d rather than a volute. It collects the flow\nwithout conserving angular momentum, so its sizing parameter is whatever\nthe chosen area happens to produce and its tangential-velocity loss is\ncorrespondingly larger. Aungier\u2019s Eq. 5-108 additionally averages that loss\nagainst a total loss of tangential head, because a uniform passage tries to\ndiffuse the entering flow to rest where collection begins.</p></dd>", "a[href=\"#pyskyfire.ducts.area_law.ConstantArea.__post_init__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.area_law.ConstantArea.__post_init__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__post_init__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire-ducts-area-law-constantarea\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.ducts.area_law.ConstantArea<a class=\"headerlink\" href=\"#pyskyfire-ducts-area-law-constantarea\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.ducts.area_law.ConstantArea.residual\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.area_law.ConstantArea.residual\">\n<span class=\"sig-name descname\"><span class=\"pre\">residual</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">theta</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">area</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">centroid_radius</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">inlet_width</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">inlet</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Return the signed miss at one azimuth.</p><p>Negative means the trial section is too small. Must increase with\nsection size so that the solver can bracket a root.</p></dd>"}
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
