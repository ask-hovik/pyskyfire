selector_to_html = {"a[href=\"#pyskyfire.pump.impeller.ImpellerDiagnostics\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.pump.impeller.ImpellerDiagnostics\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.pump.impeller.</span></span><span class=\"sig-name descname\"><span class=\"pre\">ImpellerDiagnostics</span></span></dt><dd><p>How the closed design sits against G\u00fclich\u2019s guidance. Nothing is enforced.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">solved</span></code> names what the head was closed on, <code class=\"docutils literal notranslate\"><span class=\"pre\">\"beta2\"</span></code> or <code class=\"docutils literal notranslate\"><span class=\"pre\">\"d2\"</span></code>.\n<code class=\"docutils literal notranslate\"><span class=\"pre\">psi</span></code> is the resolved head coefficient 2gH/u2\u00b2 and <code class=\"docutils literal notranslate\"><span class=\"pre\">psi_correlation</span></code>\nwhat Eq. (3.26) proposes at this specific speed and <code class=\"docutils literal notranslate\"><span class=\"pre\">fT</span></code>. <code class=\"docutils literal notranslate\"><span class=\"pre\">slip</span></code> is the\nmean-streamline prediction. <code class=\"docutils literal notranslate\"><span class=\"pre\">deviation_deg</span></code> is beta2B minus the flow\nangle, which G\u00fclich suggests keeping below 10 to 14 degrees, and\n<code class=\"docutils literal notranslate\"><span class=\"pre\">deceleration_ratio</span></code> is w2 over w1 at the shroud, whose de Haller limit\nhe puts at 0.7.</p></dd>", "a[href=\"#pyskyfire-pump-impeller-impellerdiagnostics\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.pump.impeller.ImpellerDiagnostics<a class=\"headerlink\" href=\"#pyskyfire-pump-impeller-impellerdiagnostics\" title=\"Link to this heading\">\u00b6</a></h1>"}
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
