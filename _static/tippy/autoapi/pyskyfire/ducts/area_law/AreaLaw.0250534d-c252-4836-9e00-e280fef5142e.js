selector_to_html = {"a[href=\"#pyskyfire-ducts-area-law-arealaw\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.ducts.area_law.AreaLaw<a class=\"headerlink\" href=\"#pyskyfire-ducts-area-law-arealaw\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.ducts.area_law.AreaLaw.residual\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.area_law.AreaLaw.residual\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">abstractmethod</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">residual</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">theta</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">area</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">centroid_radius</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">inlet_width</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">inlet</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Return the signed miss at one azimuth.</p><p>Negative means the trial section is too small. Must increase with\nsection size so that the solver can bracket a root.</p></dd>", "a[href=\"#pyskyfire.ducts.area_law.AreaLaw\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.area_law.AreaLaw\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.ducts.area_law.</span></span><span class=\"sig-name descname\"><span class=\"pre\">AreaLaw</span></span></dt><dd><p>Bases: <a class=\"reference external\" href=\"https://docs.python.org/3/library/abc.html#abc.ABC\" title=\"(in Python v3.15)\"><code class=\"xref py py-obj docutils literal notranslate\"><span class=\"pre\">abc.ABC</span></code></a></p><p>A rule relating azimuth to passage size.</p></dd>"}
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
