selector_to_html = {"a[href=\"#pyskyfire-viz-components-thrustchambercomponent\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.components.ThrustChamberComponent<a class=\"headerlink\" href=\"#pyskyfire-viz-components-thrustchambercomponent\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.viz.components.ThrustChamberComponent.build\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.components.ThrustChamberComponent.build\">\n<span class=\"sig-name descname\"><span class=\"pre\">build</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.viz.components.ThrustChamberComponent\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.components.ThrustChamberComponent\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.components.</span></span><span class=\"sig-name descname\"><span class=\"pre\">ThrustChamberComponent</span></span></dt><dd><p>Cooling-channel surfaces with optional RegenResult fields.</p><p>Each layout seed is constructed once. Instances share its immutable mesh\nand fields and carry separate rotations and identities. No hot wall is\nconstructed. <code class=\"docutils literal notranslate\"><span class=\"pre\">circuit_N.inlet/outlet</span></code> anchors refer to the first seed\u2019s\nfirst channel, at its end-section polygon centre, not a manifold port.\nEvery channel also exposes its own inlet/outlet anchors.</p><p>A circuit with a result records its coolant\u2019s CoolProp name as\n<code class=\"docutils literal notranslate\"><span class=\"pre\">metadata[\"fluid\"]</span></code>. <code class=\"docutils literal notranslate\"><span class=\"pre\">circuit_colors</span></code> paints each circuit apart;\nwithout it every channel takes the default colour.</p></dd>"}
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
