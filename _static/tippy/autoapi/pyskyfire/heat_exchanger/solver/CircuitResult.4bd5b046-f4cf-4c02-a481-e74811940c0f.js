selector_to_html = {"a[href=\"#pyskyfire.heat_exchanger.solver.CircuitResult\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.heat_exchanger.solver.CircuitResult\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.heat_exchanger.solver.</span></span><span class=\"sig-name descname\"><span class=\"pre\">CircuitResult</span></span></dt><dd><p>One representative channel; heat_rate is TOTAL bank heat per segment.</p></dd>", "a[href=\"#pyskyfire.heat_exchanger.solver.CircuitResult.energy_change\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.heat_exchanger.solver.CircuitResult.energy_change\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">energy_change</span></span></dt><dd></dd>", "a[href=\"#pyskyfire.heat_exchanger.solver.CircuitResult.total_heat\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.heat_exchanger.solver.CircuitResult.total_heat\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">total_heat</span></span></dt><dd></dd>", "a[href=\"#pyskyfire.heat_exchanger.solver.CircuitResult.p\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.heat_exchanger.solver.CircuitResult.p\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">p</span></span></dt><dd></dd>", "a[href=\"#pyskyfire.heat_exchanger.solver.CircuitResult.s\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.heat_exchanger.solver.CircuitResult.s\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">s</span></span></dt><dd></dd>", "a[href=\"#pyskyfire.heat_exchanger.solver.CircuitResult.outlet\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.heat_exchanger.solver.CircuitResult.outlet\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">outlet</span></span></dt><dd></dd>", "a[href=\"#pyskyfire.heat_exchanger.solver.CircuitResult.T\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.heat_exchanger.solver.CircuitResult.T\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">T</span></span></dt><dd></dd>", "a[href=\"#pyskyfire-heat-exchanger-solver-circuitresult\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.heat_exchanger.solver.CircuitResult<a class=\"headerlink\" href=\"#pyskyfire-heat-exchanger-solver-circuitresult\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.heat_exchanger.solver.CircuitResult.h\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.heat_exchanger.solver.CircuitResult.h\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">h</span></span></dt><dd></dd>"}
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
