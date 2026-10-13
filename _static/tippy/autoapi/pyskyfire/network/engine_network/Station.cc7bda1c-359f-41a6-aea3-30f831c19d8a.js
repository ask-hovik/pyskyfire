selector_to_html = {"a[href=\"#pyskyfire.network.engine_network.Station\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.network.engine_network.Station\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.network.engine_network.</span></span><span class=\"sig-name descname\"><span class=\"pre\">Station</span></span></dt><dd><p>Thermodynamic state container used by the engine network.</p></dd>", "a[href=\"#pyskyfire-network-engine-network-station\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.network.engine_network.Station<a class=\"headerlink\" href=\"#pyskyfire-network-engine-network-station\" title=\"Link to this heading\">\u00b6</a></h1>"}
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
