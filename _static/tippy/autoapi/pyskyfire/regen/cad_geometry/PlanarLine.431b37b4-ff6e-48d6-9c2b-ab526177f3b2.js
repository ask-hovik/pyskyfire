selector_to_html = {"a[href=\"#pyskyfire.regen.cad_geometry.PlanarLine\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.cad_geometry.PlanarLine\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.regen.cad_geometry.</span></span><span class=\"sig-name descname\"><span class=\"pre\">PlanarLine</span></span></dt><dd><p>Directed line, sharing the sampling/export contract of PlanarArc.</p></dd>", "a[href=\"#pyskyfire.regen.cad_geometry.PlanarLine.to_data\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.cad_geometry.PlanarLine.to_data\">\n<span class=\"sig-name descname\"><span class=\"pre\">to_data</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.regen.cad_geometry.PlanarLine.tangent\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.cad_geometry.PlanarLine.tangent\">\n<span class=\"sig-name descname\"><span class=\"pre\">tangent</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">fraction</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.regen.cad_geometry.PlanarLine.trimmed\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.cad_geometry.PlanarLine.trimmed\">\n<span class=\"sig-name descname\"><span class=\"pre\">trimmed</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">left</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">right</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.regen.cad_geometry.PlanarLine.point\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.cad_geometry.PlanarLine.point\">\n<span class=\"sig-name descname\"><span class=\"pre\">point</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">fraction</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.regen.cad_geometry.PlanarLine.transformed\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.cad_geometry.PlanarLine.transformed\">\n<span class=\"sig-name descname\"><span class=\"pre\">transformed</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">matrix</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">origin</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire-regen-cad-geometry-planarline\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.regen.cad_geometry.PlanarLine<a class=\"headerlink\" href=\"#pyskyfire-regen-cad-geometry-planarline\" title=\"Link to this heading\">\u00b6</a></h1>"}
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
