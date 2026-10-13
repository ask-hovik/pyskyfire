selector_to_html = {"a[href=\"#pyskyfire-turbine-axial-rowsettings\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.turbine.axial.RowSettings<a class=\"headerlink\" href=\"#pyskyfire-turbine-axial-rowsettings\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"../definition/index.html#pyskyfire.turbine.definition.metal_angles\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.turbine.definition.metal_angles\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.turbine.definition.</span></span><span class=\"sig-name descname\"><span class=\"pre\">metal_angles</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">flow_in</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">flow_out</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">incidence</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">0.0</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">deviation</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">0.0</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Metal angles for a row turning the flow from <code class=\"docutils literal notranslate\"><span class=\"pre\">flow_in</span></code> to <code class=\"docutils literal notranslate\"><span class=\"pre\">flow_out</span></code>.</p><p>Incidence is added to the inlet flow angle. Deviation is added to the\noutlet in the direction the row turns the flow, since the leaving flow\nfalls short of the metal and the blade has to turn it further.</p></dd>", "a[href=\"#pyskyfire.turbine.axial.RowSettings.__post_init__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.turbine.axial.RowSettings.__post_init__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__post_init__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.turbine.axial.RowSettings\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.turbine.axial.RowSettings\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.turbine.axial.</span></span><span class=\"sig-name descname\"><span class=\"pre\">RowSettings</span></span></dt><dd><p>How the triangles become blades.</p></dd>"}
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
