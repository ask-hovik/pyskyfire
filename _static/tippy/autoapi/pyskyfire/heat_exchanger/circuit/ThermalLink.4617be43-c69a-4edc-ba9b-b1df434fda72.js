selector_to_html = {"a[href=\"#pyskyfire-heat-exchanger-circuit-thermallink\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.heat_exchanger.circuit.ThermalLink<a class=\"headerlink\" href=\"#pyskyfire-heat-exchanger-circuit-thermallink\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.heat_exchanger.circuit.ThermalLink\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.heat_exchanger.circuit.ThermalLink\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.heat_exchanger.circuit.</span></span><span class=\"sig-name descname\"><span class=\"pre\">ThermalLink</span></span></dt><dd><p>Pair finite volumes on two paths through ONE shared wall.</p><p>Each (i,j) pair couples the segment beginning at station i on a to the\nsegment beginning at station j on b. Reverse j for counterflow. All copied\nchannels are paired one-to-one. Contact areas [m\u00b2 per pair per copy] may be\nexplicit, otherwise use the smaller exposed surface area. Wall thickness\nis evaluated at the midpoint s on circuit a; it is never counted twice.</p></dd>"}
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
