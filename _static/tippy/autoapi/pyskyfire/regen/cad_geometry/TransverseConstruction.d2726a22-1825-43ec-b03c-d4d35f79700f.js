selector_to_html = {"a[href=\"#pyskyfire.regen.cad_geometry.TransverseConstruction.polygon\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.cad_geometry.TransverseConstruction.polygon\">\n<span class=\"sig-name descname\"><span class=\"pre\">polygon</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">index</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">n_arc_points</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">8</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire-regen-cad-geometry-transverseconstruction\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.regen.cad_geometry.TransverseConstruction<a class=\"headerlink\" href=\"#pyskyfire-regen-cad-geometry-transverseconstruction\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.regen.cad_geometry.TransverseConstruction\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.cad_geometry.TransverseConstruction\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.regen.cad_geometry.</span></span><span class=\"sig-name descname\"><span class=\"pre\">TransverseConstruction</span></span></dt><dd></dd>"}
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
