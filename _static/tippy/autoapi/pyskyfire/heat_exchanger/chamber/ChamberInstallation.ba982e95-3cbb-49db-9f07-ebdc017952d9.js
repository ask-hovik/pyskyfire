selector_to_html = {"a[href=\"#pyskyfire.heat_exchanger.chamber.ChamberInstallation.solve\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.heat_exchanger.chamber.ChamberInstallation.solve\">\n<span class=\"sig-name descname\"><span class=\"pre\">solve</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">inlets</span></span></em>, <em class=\"sig-param\"><span class=\"o\"><span class=\"pre\">**</span></span><span class=\"n\"><span class=\"pre\">kwargs</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire-heat-exchanger-chamber-chamberinstallation\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.heat_exchanger.chamber.ChamberInstallation<a class=\"headerlink\" href=\"#pyskyfire-heat-exchanger-chamber-chamberinstallation\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.heat_exchanger.chamber.ChamberInstallation._volume_change\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.heat_exchanger.chamber.ChamberInstallation._volume_change\">\n<span class=\"sig-name descname\"><span class=\"pre\">_volume_change</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.heat_exchanger.chamber.ChamberInstallation\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.heat_exchanger.chamber.ChamberInstallation\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.heat_exchanger.chamber.</span></span><span class=\"sig-name descname\"><span class=\"pre\">ChamberInstallation</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">chamber</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">exchanger</span></span></em>, <em class=\"sig-param\"><span class=\"keyword-only-separator o\"><abbr title=\"Keyword-only parameters separator (PEP 3102)\"><span class=\"pre\">*</span></abbr></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">position</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">(0.0,</span> <span class=\"pre\">0.0,</span> <span class=\"pre\">0.0)</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">rotation</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">None</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">coefficient</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">None</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>An exchanger in chamber coordinates, without gas-solution feedback.</p><p>Volume is the sum of nonoverlapping swept outer envelopes, including\npassages. Users must avoid overlapping copies. Partially intersecting\nsections use deterministic triangle/axial quadrature, not a CAD Boolean.\n<code class=\"docutils literal notranslate\"><span class=\"pre\">volume_change.chamber_displaced</span></code> counts only x &lt;= contour.x_t;\n<code class=\"docutils literal notranslate\"><span class=\"pre\">displaced</span></code> includes nozzle volume too.</p></dd>"}
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
