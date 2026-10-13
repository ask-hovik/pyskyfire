selector_to_html = {"a[href=\"#pyskyfire.viz.components.AnnularDuctComponent\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.components.AnnularDuctComponent\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.components.</span></span><span class=\"sig-name descname\"><span class=\"pre\">AnnularDuctComponent</span></span></dt><dd><p>Display an <code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">AnnularDuct</span></code>.</p><p>The duct is drawn as its hub and shroud revolved about the local +x\nshaft, with <code class=\"docutils literal notranslate\"><span class=\"pre\">inlet</span></code> and <code class=\"docutils literal notranslate\"><span class=\"pre\">outlet</span></code> anchors on the axis at its ends as a\nturbine\u2019s are. The shroud is translucent so the hub shows. A pipe with a\ncore is drawn by <code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">PipeComponent</span></code>.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">results</span></code> is a <code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">DuctResult</span></code>, or a mapping of\nfield names to FieldData with one value per station (point) or per\ninterval (cell). Both walls carry the same fields.</p></dd>", "a[href=\"#pyskyfire.viz.components.AnnularDuctComponent.build\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.components.AnnularDuctComponent.build\">\n<span class=\"sig-name descname\"><span class=\"pre\">build</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire-viz-components-annularductcomponent\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.components.AnnularDuctComponent<a class=\"headerlink\" href=\"#pyskyfire-viz-components-annularductcomponent\" title=\"Link to this heading\">\u00b6</a></h1>"}
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
