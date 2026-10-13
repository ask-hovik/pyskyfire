selector_to_html = {"a[href=\"#submodules\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Submodules<a class=\"headerlink\" href=\"#submodules\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"solver/index.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.heat_exchanger.solver<a class=\"headerlink\" href=\"#module-pyskyfire.heat_exchanger.solver\" title=\"Link to this heading\">\u00b6</a></h1><p>Steady single-phase heat exchange with conservative segment heat rates.</p><p>A circuit solution represents identical parallel copies. Static pressure obeys\nDarcy friction plus rho*u*du; energy conserves h + u\u00b2/2. Exchanger heat rates\nare iterated against segment midpoint states, including counterflow links.</p>", "a[href=\"#module-pyskyfire.heat_exchanger\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.heat_exchanger<a class=\"headerlink\" href=\"#module-pyskyfire.heat_exchanger\" title=\"Link to this heading\">\u00b6</a></h1><p>Path-based heat exchangers with prescribed gas or coupled fluid circuits.</p>", "a[href=\"chamber/index.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.heat_exchanger.chamber<a class=\"headerlink\" href=\"#module-pyskyfire.heat_exchanger.chamber\" title=\"Link to this heading\">\u00b6</a></h1><p>Placement and volume accounting inside a prescribed thrust-chamber flow.</p>", "a[href=\"circuit/index.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.heat_exchanger.circuit<a class=\"headerlink\" href=\"#module-pyskyfire.heat_exchanger.circuit\" title=\"Link to this heading\">\u00b6</a></h1><p>Reusable channel banks, finite-volume thermal links, and prescribed gas.</p>"}
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
