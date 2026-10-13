selector_to_html = {"a[href=\"#pyskyfire.viz.network_viz._BlockSpec.is_signal_only\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.network_viz._BlockSpec.is_signal_only\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">is_signal_only</span></span><span class=\"property\"><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/functions.html#bool\" title=\"(in Python v3.15)\"><span class=\"pre\">bool</span></a></span></dt><dd></dd>", "a[href=\"#pyskyfire-viz-network-viz-blockspec\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.network_viz._BlockSpec<a class=\"headerlink\" href=\"#pyskyfire-viz-network-viz-blockspec\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.viz.network_viz._BlockSpec\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.network_viz._BlockSpec\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.network_viz.</span></span><span class=\"sig-name descname\"><span class=\"pre\">_BlockSpec</span></span></dt><dd></dd>"}
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
