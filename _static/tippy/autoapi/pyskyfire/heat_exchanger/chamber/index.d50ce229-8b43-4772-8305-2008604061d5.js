selector_to_html = {"a[href=\"#module-pyskyfire.heat_exchanger.chamber\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.heat_exchanger.chamber<a class=\"headerlink\" href=\"#module-pyskyfire.heat_exchanger.chamber\" title=\"Link to this heading\">\u00b6</a></h1><p>Placement and volume accounting inside a prescribed thrust-chamber flow.</p>", "a[href=\"ChamberInstallation.html#pyskyfire.heat_exchanger.chamber.ChamberInstallation\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.heat_exchanger.chamber.ChamberInstallation\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.heat_exchanger.chamber.</span></span><span class=\"sig-name descname\"><span class=\"pre\">ChamberInstallation</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">chamber</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">exchanger</span></span></em>, <em class=\"sig-param\"><span class=\"keyword-only-separator o\"><abbr title=\"Keyword-only parameters separator (PEP 3102)\"><span class=\"pre\">*</span></abbr></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">position</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">(0.0,</span> <span class=\"pre\">0.0,</span> <span class=\"pre\">0.0)</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">rotation</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">None</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">coefficient</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">None</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>An exchanger in chamber coordinates, without gas-solution feedback.</p><p>Volume is the sum of nonoverlapping swept outer envelopes, including\npassages. Users must avoid overlapping copies. Partially intersecting\nsections use deterministic triangle/axial quadrature, not a CAD Boolean.\n<code class=\"docutils literal notranslate\"><span class=\"pre\">volume_change.chamber_displaced</span></code> counts only x &lt;= contour.x_t;\n<code class=\"docutils literal notranslate\"><span class=\"pre\">displaced</span></code> includes nozzle volume too.</p></dd>", "a[href=\"VolumeChange.html#pyskyfire.heat_exchanger.chamber.VolumeChange\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.heat_exchanger.chamber.VolumeChange\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.heat_exchanger.chamber.</span></span><span class=\"sig-name descname\"><span class=\"pre\">VolumeChange</span></span></dt><dd></dd>", "a[href=\"#classes\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Classes<a class=\"headerlink\" href=\"#classes\" title=\"Link to this heading\">\u00b6</a></h2>"}
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
