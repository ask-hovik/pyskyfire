selector_to_html = {"a[href=\"#pyskyfire-injector-injector-setresult\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.injector.injector.SetResult<a class=\"headerlink\" href=\"#pyskyfire-injector-injector-setresult\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.injector.injector.SetResult\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.injector.injector.SetResult\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.injector.injector.</span></span><span class=\"sig-name descname\"><span class=\"pre\">SetResult</span></span></dt><dd><p>The flow through one element set.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">mass_flows</span></code> maps each of the set\u2019s feeds to the set\u2019s total from it\n[kg/s]. <code class=\"docutils literal notranslate\"><span class=\"pre\">diagnostics</span></code> is the element type\u2019s own mixing report for one\nelement; every element in a set runs identically.</p></dd>", "a[href=\"#pyskyfire.injector.injector.SetResult.element_mass_flow\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.injector.injector.SetResult.element_mass_flow\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">element_mass_flow</span></span><span class=\"property\"><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/functions.html#float\" title=\"(in Python v3.15)\"><span class=\"pre\">float</span></a></span></dt><dd><p>Mass flow through one element [kg/s].</p></dd>", "a[href=\"#pyskyfire.injector.injector.SetResult.mass_flow\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.injector.injector.SetResult.mass_flow\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">mass_flow</span></span><span class=\"property\"><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/functions.html#float\" title=\"(in Python v3.15)\"><span class=\"pre\">float</span></a></span></dt><dd><p>Total mass flow through the set [kg/s].</p></dd>"}
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
