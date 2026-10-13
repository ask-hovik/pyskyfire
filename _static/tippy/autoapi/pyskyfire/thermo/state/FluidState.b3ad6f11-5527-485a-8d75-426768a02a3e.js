selector_to_html = {"a[href=\"#pyskyfire.thermo.state.FluidState.is_two_phase\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.thermo.state.FluidState.is_two_phase\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">is_two_phase</span></span><span class=\"property\"><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/functions.html#bool\" title=\"(in Python v3.15)\"><span class=\"pre\">bool</span></a></span></dt><dd><p>True when the bulk sits inside the boiling dome.</p></dd>", "a[href=\"#pyskyfire.thermo.state.FluidState\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.thermo.state.FluidState\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.thermo.state.</span></span><span class=\"sig-name descname\"><span class=\"pre\">FluidState</span></span></dt><dd><p>Resolved bulk fluid state at one axial station.</p></dd>", "a[href=\"#pyskyfire-thermo-state-fluidstate\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.thermo.state.FluidState<a class=\"headerlink\" href=\"#pyskyfire-thermo-state-fluidstate\" title=\"Link to this heading\">\u00b6</a></h1>"}
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
