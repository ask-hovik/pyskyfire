selector_to_html = {"a[href=\"#pyskyfire.turbine.axial.StageFlow\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.turbine.axial.StageFlow\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.turbine.axial.</span></span><span class=\"sig-name descname\"><span class=\"pre\">StageFlow</span></span></dt><dd><p>One stage\u2019s triangles: stator inlet (1), rotor inlet (2), rotor exit (3).</p></dd>", "a[href=\"#pyskyfire.turbine.axial.StageFlow.flow_coefficient\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.turbine.axial.StageFlow.flow_coefficient\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">flow_coefficient</span></span><span class=\"property\"><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/functions.html#float\" title=\"(in Python v3.15)\"><span class=\"pre\">float</span></a></span></dt><dd></dd>", "a[href=\"#pyskyfire.turbine.axial.StageFlow.work\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.turbine.axial.StageFlow.work\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">work</span></span><span class=\"property\"><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/functions.html#float\" title=\"(in Python v3.15)\"><span class=\"pre\">float</span></a></span></dt><dd><p>Euler work <code class=\"docutils literal notranslate\"><span class=\"pre\">u*(c_u2</span> <span class=\"pre\">-</span> <span class=\"pre\">c_u3)</span></code> [J kg\u207b\u00b9].</p></dd>", "a[href=\"#pyskyfire-turbine-axial-stageflow\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.turbine.axial.StageFlow<a class=\"headerlink\" href=\"#pyskyfire-turbine-axial-stageflow\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.turbine.axial.StageFlow.loading\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.turbine.axial.StageFlow.loading\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">loading</span></span><span class=\"property\"><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/functions.html#float\" title=\"(in Python v3.15)\"><span class=\"pre\">float</span></a></span></dt><dd><p><code class=\"docutils literal notranslate\"><span class=\"pre\">\u0394h0/u\u00b2</span></code>.</p></dd>", "a[href=\"#pyskyfire.turbine.axial.StageFlow.reaction\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.turbine.axial.StageFlow.reaction\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">reaction</span></span><span class=\"property\"><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/functions.html#float\" title=\"(in Python v3.15)\"><span class=\"pre\">float</span></a></span></dt><dd><p>Rotor share of the stage\u2019s static enthalpy drop.</p></dd>"}
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
