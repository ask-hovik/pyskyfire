selector_to_html = {"a[href=\"#pyskyfire.pump._pyturbo.build_geometry\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.pump._pyturbo.build_geometry\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.pump._pyturbo.</span></span><span class=\"sig-name descname\"><span class=\"pre\">build_geometry</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">definition</span></span></em>, <em class=\"sig-param\"><span class=\"keyword-only-separator o\"><abbr title=\"Keyword-only parameters separator (PEP 3102)\"><span class=\"pre\">*</span></abbr></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">span_points</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">chord_points</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#module-pyskyfire.pump._pyturbo\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.pump._pyturbo<a class=\"headerlink\" href=\"#module-pyskyfire.pump._pyturbo\" title=\"Link to this heading\">\u00b6</a></h1><p>Translate physical pump inputs to pyTurbo 1.3.8 construction controls.</p><p>The conformal (m\u2019,theta) plane uses dm\u2019=dm/r. Metal angles are from +m\u2019,\nand upstream wrap_angle means atan(dtheta/m\u2019), not total blade wrap.\nThickness controls are section-normal offsets in this plane: locally t/(2r).</p>", "a[href=\"#functions\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Functions<a class=\"headerlink\" href=\"#functions\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"#module-contents\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Module Contents<a class=\"headerlink\" href=\"#module-contents\" title=\"Link to this heading\">\u00b6</a></h2>"}
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
