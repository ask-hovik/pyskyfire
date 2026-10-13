selector_to_html = {"a[href=\"#pyskyfire.injector.impinging.ImpingingDiagnostics\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.injector.impinging.ImpingingDiagnostics\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.injector.impinging.</span></span><span class=\"sig-name descname\"><span class=\"pre\">ImpingingDiagnostics</span></span></dt><dd><p>Mixing quantities of one impinging element.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">momentum</span></code> is the magnitude of the summed jet momentum [N] and\n<code class=\"docutils literal notranslate\"><span class=\"pre\">resultant_angle</span></code> its angle from the chamber axis [rad].\n<code class=\"docutils literal notranslate\"><span class=\"pre\">mixture_ratio</span></code> is the element\u2019s oxidizer over fuel mass flow, NaN\nwithout fuel or without any oxidizer feed named. <code class=\"docutils literal notranslate\"><span class=\"pre\">rupe</span></code> and\n<code class=\"docutils literal notranslate\"><span class=\"pre\">momentum_flux_ratio</span></code> are oxidizer over fuel, and NaN unless the element\nhas exactly one fuel and one oxidizer jet.</p></dd>", "a[href=\"#pyskyfire-injector-impinging-impingingdiagnostics\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.injector.impinging.ImpingingDiagnostics<a class=\"headerlink\" href=\"#pyskyfire-injector-impinging-impingingdiagnostics\" title=\"Link to this heading\">\u00b6</a></h1>"}
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
