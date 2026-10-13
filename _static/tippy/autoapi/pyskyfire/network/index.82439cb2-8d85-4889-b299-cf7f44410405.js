selector_to_html = {"a[href=\"#submodules\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Submodules<a class=\"headerlink\" href=\"#submodules\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"engine_network/index.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.network.engine_network<a class=\"headerlink\" href=\"#module-pyskyfire.network.engine_network\" title=\"Link to this heading\">\u00b6</a></h1><h2>Classes<a class=\"headerlink\" href=\"#classes\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"blocks/index.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.network.blocks<a class=\"headerlink\" href=\"#module-pyskyfire.network.blocks\" title=\"Link to this heading\">\u00b6</a></h1><h2>Classes<a class=\"headerlink\" href=\"#classes\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"#module-pyskyfire.network\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.network<a class=\"headerlink\" href=\"#module-pyskyfire.network\" title=\"Link to this heading\">\u00b6</a></h1><h2>Submodules<a class=\"headerlink\" href=\"#submodules\" title=\"Link to this heading\">\u00b6</a></h2>"}
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
