selector_to_html = {"a[href=\"#pyskyfire-viz-components-impellercomponent\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.components.ImpellerComponent<a class=\"headerlink\" href=\"#pyskyfire-viz-components-impellercomponent\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.viz.components.ImpellerComponent.build\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.components.ImpellerComponent.build\">\n<span class=\"sig-name descname\"><span class=\"pre\">build</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.viz.components.ImpellerComponent\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.components.ImpellerComponent\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.components.</span></span><span class=\"sig-name descname\"><span class=\"pre\">ImpellerComponent</span></span></dt><dd><p>Display prepared impeller.geometry arrays, with the existing local x shaft axis.</p><p>Optional <code class=\"docutils literal notranslate\"><span class=\"pre\">results</span></code> maps part names (hub, shroud, blade_N) to dictionaries\nof mesh-aligned FieldData. <code class=\"docutils literal notranslate\"><span class=\"pre\">inlet</span></code> and <code class=\"docutils literal notranslate\"><span class=\"pre\">shaft</span></code> are axial centre datums;\nthere is no discrete outlet port on a radially discharging impeller.</p></dd>"}
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
