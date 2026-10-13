selector_to_html = {"a[href=\"#pyskyfire-regen-coupled-solver-regenresult\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.regen.coupled_solver.RegenResult<a class=\"headerlink\" href=\"#pyskyfire-regen-coupled-solver-regenresult\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.regen.coupled_solver.RegenResult\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.coupled_solver.RegenResult\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.regen.coupled_solver.</span></span><span class=\"sig-name descname\"><span class=\"pre\">RegenResult</span></span></dt><dd><p>Regenerative-cooling result with optional film-cooling data.</p><p>The coolant-state fields (<code class=\"docutils literal notranslate\"><span class=\"pre\">coolant_enthalpy</span></code> onwards) describe the bulk\ncoolant along the circuit. <code class=\"docutils literal notranslate\"><span class=\"pre\">coolant_enthalpy</span></code> and <code class=\"docutils literal notranslate\"><span class=\"pre\">coolant_quality</span></code> are\nall-<code class=\"docutils literal notranslate\"><span class=\"pre\">nan</span></code> when the coolant fell back to the temperature march, and\n<code class=\"docutils literal notranslate\"><span class=\"pre\">coolant_quality</span></code> is <code class=\"docutils literal notranslate\"><span class=\"pre\">nan</span></code> at any station outside the boiling dome \u2013\nincluding supercritical stations, where no dome exists. Use\n<code class=\"docutils literal notranslate\"><span class=\"pre\">coolant_phase</span></code> to tell those cases apart.</p></dd>"}
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
