selector_to_html = {"a[href=\"#module-pyskyfire.common.momentum\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.common.momentum<a class=\"headerlink\" href=\"#module-pyskyfire.common.momentum\" title=\"Link to this heading\">\u00b6</a></h1><p>One-dimensional static-pressure changes shared by heated duct solvers.</p>", "a[href=\"#functions\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Functions<a class=\"headerlink\" href=\"#functions\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"#module-contents\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Module Contents<a class=\"headerlink\" href=\"#module-contents\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"#pyskyfire.common.momentum.acceleration_pressure_drop\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.momentum.acceleration_pressure_drop\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.common.momentum.</span></span><span class=\"sig-name descname\"><span class=\"pre\">acceleration_pressure_drop</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">rho_up</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">velocity_up</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">rho_down</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">velocity_down</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Trapezoidal integral of rho*u*du [Pa]; positive for acceleration.</p><p>Exact for constant density (Bernoulli) and constant area (mass-flux form).\nThis is reversible static-pressure change, separate from wall friction.</p></dd>"}
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
