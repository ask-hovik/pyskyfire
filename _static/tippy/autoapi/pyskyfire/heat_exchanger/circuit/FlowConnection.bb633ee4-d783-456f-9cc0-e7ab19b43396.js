selector_to_html = {"a[href=\"#pyskyfire-heat-exchanger-circuit-flowconnection\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.heat_exchanger.circuit.FlowConnection<a class=\"headerlink\" href=\"#pyskyfire-heat-exchanger-circuit-flowconnection\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.heat_exchanger.circuit.FlowConnection\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.heat_exchanger.circuit.FlowConnection\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.heat_exchanger.circuit.</span></span><span class=\"sig-name descname\"><span class=\"pre\">FlowConnection</span></span></dt><dd><p>An ideal series header: all outlet flow enters the next bank.</p><p>The banks must have the same number of copies and use the same fluid.\nThe header has no heat or fitting loss; connecting geometry is not swept.</p></dd>"}
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
