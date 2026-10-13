selector_to_html = {"a[href=\"#pyskyfire.heat_exchanger.circuit.PrescribedGas\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.heat_exchanger.circuit.PrescribedGas\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.heat_exchanger.circuit.</span></span><span class=\"sig-name descname\"><span class=\"pre\">PrescribedGas</span></span></dt><dd><p>Fixed surrounding flow: callbacks receive assembly position and wall T.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">temperature</span></code> may be a constant or <code class=\"docutils literal notranslate\"><span class=\"pre\">temperature(position)</span></code>; coefficient\nmay be a constant or <code class=\"docutils literal notranslate\"><span class=\"pre\">coefficient(position,</span> <span class=\"pre\">wall_temperature)</span></code> [W/m\u00b2/K].</p></dd>", "a[href=\"#pyskyfire.heat_exchanger.circuit.PrescribedGas.sample\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.heat_exchanger.circuit.PrescribedGas.sample\">\n<span class=\"sig-name descname\"><span class=\"pre\">sample</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">position</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">wall_temperature</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire-heat-exchanger-circuit-prescribedgas\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.heat_exchanger.circuit.PrescribedGas<a class=\"headerlink\" href=\"#pyskyfire-heat-exchanger-circuit-prescribedgas\" title=\"Link to this heading\">\u00b6</a></h1>"}
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
