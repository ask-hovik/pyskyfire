selector_to_html = {"a[href=\"#pyskyfire.thermo.state.SaturationState\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.thermo.state.SaturationState\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.thermo.state.</span></span><span class=\"sig-name descname\"><span class=\"pre\">SaturationState</span></span></dt><dd><p>Saturation properties at one pressure.</p><p>Built once per axial station and cached, so that the inner wall solve \u2013\nwhich may evaluate the fluid side a couple of hundred times per node \u2013\nnever triggers a flash.</p></dd>", "a[href=\"#pyskyfire-thermo-state-saturationstate\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.thermo.state.SaturationState<a class=\"headerlink\" href=\"#pyskyfire-thermo-state-saturationstate\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.thermo.state.SaturationState.homogeneous_density\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.thermo.state.SaturationState.homogeneous_density\">\n<span class=\"sig-name descname\"><span class=\"pre\">homogeneous_density</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">quality</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/functions.html#float\" title=\"(in Python v3.15)\"><span class=\"pre\">float</span></a></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/functions.html#float\" title=\"(in Python v3.15)\"><span class=\"pre\">float</span></a></span></span></dt><dd><p>Homogeneous (no-slip) two-phase density [kg m\u207b\u00b3].</p></dd>", "a[href=\"#pyskyfire.thermo.state.SaturationState.h_fg\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.thermo.state.SaturationState.h_fg\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">h_fg</span></span><span class=\"property\"><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/functions.html#float\" title=\"(in Python v3.15)\"><span class=\"pre\">float</span></a></span></dt><dd><p>Latent heat of vaporisation [J kg\u207b\u00b9].</p></dd>"}
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
