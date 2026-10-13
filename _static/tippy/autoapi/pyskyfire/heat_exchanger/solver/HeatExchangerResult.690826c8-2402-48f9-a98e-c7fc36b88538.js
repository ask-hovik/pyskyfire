selector_to_html = {"a[href=\"#pyskyfire.heat_exchanger.solver.HeatExchangerResult.energy_balance_error\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.heat_exchanger.solver.HeatExchangerResult.energy_balance_error\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">energy_balance_error</span></span></dt><dd><p>Sum of all bank energy changes minus prescribed-gas heat [W].</p></dd>", "a[href=\"#pyskyfire.heat_exchanger.solver.HeatExchangerResult\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.heat_exchanger.solver.HeatExchangerResult\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.heat_exchanger.solver.</span></span><span class=\"sig-name descname\"><span class=\"pre\">HeatExchangerResult</span></span></dt><dd></dd>", "a[href=\"#pyskyfire-heat-exchanger-solver-heatexchangerresult\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.heat_exchanger.solver.HeatExchangerResult<a class=\"headerlink\" href=\"#pyskyfire-heat-exchanger-solver-heatexchangerresult\" title=\"Link to this heading\">\u00b6</a></h1>"}
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
