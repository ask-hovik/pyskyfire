selector_to_html = {"a[href=\"#pyskyfire.ducts.orifice.OrificeCapacityError\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.orifice.OrificeCapacityError\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">exception</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.ducts.orifice.</span></span><span class=\"sig-name descname\"><span class=\"pre\">OrificeCapacityError</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">message</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/stdtypes.html#str\" title=\"(in Python v3.15)\"><span class=\"pre\">str</span></a></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">max_mass_flow</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/functions.html#float\" title=\"(in Python v3.15)\"><span class=\"pre\">float</span></a></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Bases: <a class=\"reference external\" href=\"https://docs.python.org/3/builtins/exceptions.html#ValueError\" title=\"(in Python v3.15)\"><code class=\"xref py py-obj docutils literal notranslate\"><span class=\"pre\">ValueError</span></code></a></p><p>The inlet state cannot drive the requested mass flow through the orifice.</p><p>Either the flow chokes below the requested mass flow, or passing it would\nneed a discharge pressure below the searched floor. <code class=\"docutils literal notranslate\"><span class=\"pre\">max_mass_flow</span></code> is\nthe largest flow found [kg/s].</p></dd>", "a[href=\"#pyskyfire-ducts-orifice-orificecapacityerror\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.ducts.orifice.OrificeCapacityError<a class=\"headerlink\" href=\"#pyskyfire-ducts-orifice-orificecapacityerror\" title=\"Link to this heading\">\u00b6</a></h1>"}
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
