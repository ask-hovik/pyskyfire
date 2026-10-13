selector_to_html = {"a[href=\"#pyskyfire.pump.impeller.OutletSlip\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.pump.impeller.OutletSlip\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.pump.impeller.</span></span><span class=\"sig-name descname\"><span class=\"pre\">OutletSlip</span></span></dt><dd><p>G\u00fclich\u2019s slip prediction at the impeller outlet, Table 3.2.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">blockage</span></code> is tau_2 (T3.2.3), <code class=\"docutils literal notranslate\"><span class=\"pre\">limiting_ratio</span></code> epsilon_lim (T3.2.4),\n<code class=\"docutils literal notranslate\"><span class=\"pre\">diameter_correction</span></code> k_w (T3.2.5) and <code class=\"docutils literal notranslate\"><span class=\"pre\">slip_factor</span></code> gamma (T3.2.6).\n<code class=\"docutils literal notranslate\"><span class=\"pre\">cm</span></code> is the unblocked meridional velocity and <code class=\"docutils literal notranslate\"><span class=\"pre\">cu</span></code> the absolute swirl\nthe blades leave (T3.2.7), both in m/s.</p></dd>", "a[href=\"#pyskyfire-pump-impeller-outletslip\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.pump.impeller.OutletSlip<a class=\"headerlink\" href=\"#pyskyfire-pump-impeller-outletslip\" title=\"Link to this heading\">\u00b6</a></h1>"}
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
