selector_to_html = {"a[href=\"#pyskyfire-ducts-orifice-orificeresult\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.ducts.orifice.OrificeResult<a class=\"headerlink\" href=\"#pyskyfire-ducts-orifice-orificeresult\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.ducts.orifice.OrificeResult.pressure_loss\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.orifice.OrificeResult.pressure_loss\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">pressure_loss</span></span><span class=\"property\"><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/functions.html#float\" title=\"(in Python v3.15)\"><span class=\"pre\">float</span></a></span></dt><dd><p>Irreversible pressure loss from inlet to outlet [Pa].</p></dd>", "a[href=\"#pyskyfire.ducts.orifice.OrificeResult\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.orifice.OrificeResult\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.ducts.orifice.</span></span><span class=\"sig-name descname\"><span class=\"pre\">OrificeResult</span></span></dt><dd><p>States and diagnostics of the flow through an orifice.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">throat_state</span></code> is the ideal isentropic state at the vena contracta and\n<code class=\"docutils literal notranslate\"><span class=\"pre\">velocity</span></code> the ideal jet velocity there [m/s]. The real jet is slower by\nthe velocity coefficient, which the discharge coefficient includes but\ndoes not separate out. <code class=\"docutils literal notranslate\"><span class=\"pre\">outlet_state</span></code> is at the discharge pressure with\nthe inlet enthalpy: the jet\u2019s kinetic energy has been dissipated.\n<code class=\"docutils literal notranslate\"><span class=\"pre\">choked</span></code> is True when the discharge pressure lies below the pressure at\nwhich the ideal mass flux peaks, so lowering it further adds no flow.</p></dd>"}
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
