selector_to_html = {"a[href=\"#pyskyfire.turbine.definition.FreeVortexTwist.flow_angles\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.turbine.definition.FreeVortexTwist.flow_angles\">\n<span class=\"sig-name descname\"><span class=\"pre\">flow_angles</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">r</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Inlet and outlet flow angles, in the row\u2019s own frame, at radius <code class=\"docutils literal notranslate\"><span class=\"pre\">r</span></code>.</p></dd>", "a[href=\"#pyskyfire.turbine.definition.FreeVortexTwist\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.turbine.definition.FreeVortexTwist\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.turbine.definition.</span></span><span class=\"sig-name descname\"><span class=\"pre\">FreeVortexTwist</span></span></dt><dd><p>Metal angles from hub to tip for a free-vortex row.</p><p>In a free vortex the axial velocity is the same at every radius and the\nswirl falls off as <code class=\"docutils literal notranslate\"><span class=\"pre\">1/r</span></code>, so <code class=\"docutils literal notranslate\"><span class=\"pre\">r*c_u</span></code> is constant across the span. A\nrotor\u2019s blade speed grows with <code class=\"docutils literal notranslate\"><span class=\"pre\">r</span></code>, so its relative angles twist more\nthan a stator\u2019s absolute ones.</p></dd>", "a[href=\"#pyskyfire-turbine-definition-freevortextwist\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.turbine.definition.FreeVortexTwist<a class=\"headerlink\" href=\"#pyskyfire-turbine-definition-freevortextwist\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.turbine.definition.FreeVortexTwist.__call__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.turbine.definition.FreeVortexTwist.__call__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__call__</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">span_fraction</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">r</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>"}
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
