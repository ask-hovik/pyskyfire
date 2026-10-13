selector_to_html = {"a[href=\"#functions\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Functions<a class=\"headerlink\" href=\"#functions\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"#pyskyfire.pump.pump_cli.main\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.pump.pump_cli.main\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.pump.pump_cli.</span></span><span class=\"sig-name descname\"><span class=\"pre\">main</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#module-pyskyfire.pump.pump_cli\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.pump.pump_cli<a class=\"headerlink\" href=\"#module-pyskyfire.pump.pump_cli\" title=\"Link to this heading\">\u00b6</a></h1><p>Command-line preliminary impeller sizing, with explicit design inputs.</p>", "a[href=\"#module-contents\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Module Contents<a class=\"headerlink\" href=\"#module-contents\" title=\"Link to this heading\">\u00b6</a></h2>"}
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
