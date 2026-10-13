selector_to_html = {"a[href=\"#pyskyfire-heat-exchanger-circuit-inlet\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.heat_exchanger.circuit.Inlet<a class=\"headerlink\" href=\"#pyskyfire-heat-exchanger-circuit-inlet\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.heat_exchanger.circuit.Inlet\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.heat_exchanger.circuit.Inlet\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.heat_exchanger.circuit.</span></span><span class=\"sig-name descname\"><span class=\"pre\">Inlet</span></span></dt><dd><p>Bank inlet: temperature [K], static pressure [Pa], TOTAL mass flow [kg/s].</p></dd>", "a[href=\"#pyskyfire.heat_exchanger.circuit.Inlet.__post_init__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.heat_exchanger.circuit.Inlet.__post_init__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__post_init__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>"}
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
