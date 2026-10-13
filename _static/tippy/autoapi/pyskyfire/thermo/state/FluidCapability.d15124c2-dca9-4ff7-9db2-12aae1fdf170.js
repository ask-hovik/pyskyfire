selector_to_html = {"a[href=\"#pyskyfire.thermo.state.FluidCapability\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.thermo.state.FluidCapability\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.thermo.state.</span></span><span class=\"sig-name descname\"><span class=\"pre\">FluidCapability</span></span></dt><dd><p>What two-phase machinery is usable for a given fluid.</p></dd>", "a[href=\"#pyskyfire-thermo-state-fluidcapability\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.thermo.state.FluidCapability<a class=\"headerlink\" href=\"#pyskyfire-thermo-state-fluidcapability\" title=\"Link to this heading\">\u00b6</a></h1>"}
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
