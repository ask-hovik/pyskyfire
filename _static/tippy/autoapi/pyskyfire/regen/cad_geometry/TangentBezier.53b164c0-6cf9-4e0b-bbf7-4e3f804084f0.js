selector_to_html = {"a[href=\"#pyskyfire.regen.cad_geometry.TangentBezier.to_data\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.cad_geometry.TangentBezier.to_data\">\n<span class=\"sig-name descname\"><span class=\"pre\">to_data</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.regen.cad_geometry.TangentBezier.point\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.cad_geometry.TangentBezier.point\">\n<span class=\"sig-name descname\"><span class=\"pre\">point</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">fraction</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.regen.cad_geometry.TangentBezier.tangent\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.cad_geometry.TangentBezier.tangent\">\n<span class=\"sig-name descname\"><span class=\"pre\">tangent</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">fraction</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.regen.cad_geometry.TangentBezier\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.cad_geometry.TangentBezier\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.regen.cad_geometry.</span></span><span class=\"sig-name descname\"><span class=\"pre\">TangentBezier</span></span></dt><dd></dd>", "a[href=\"#pyskyfire-regen-cad-geometry-tangentbezier\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.regen.cad_geometry.TangentBezier<a class=\"headerlink\" href=\"#pyskyfire-regen-cad-geometry-tangentbezier\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.regen.cad_geometry.TangentBezier.between\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.cad_geometry.TangentBezier.between\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">classmethod</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">between</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">first</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">second</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>"}
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
