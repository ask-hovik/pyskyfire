selector_to_html = {"a[href=\"#pyskyfire-viz-scene-fielddata\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.scene.FieldData<a class=\"headerlink\" href=\"#pyskyfire-viz-scene-fielddata\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.viz.scene.FieldData\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.scene.FieldData\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.scene.</span></span><span class=\"sig-name descname\"><span class=\"pre\">FieldData</span></span></dt><dd><p>Scalar mesh values in base SI units. NaN represents unavailable data.</p><p>Units are descriptive metadata, never conversion inputs. Association is\n<code class=\"docutils literal notranslate\"><span class=\"pre\">point</span></code> or <code class=\"docutils literal notranslate\"><span class=\"pre\">cell</span></code> (one value per triangle).</p></dd>", "a[href=\"#pyskyfire.viz.scene.FieldData.__post_init__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.scene.FieldData.__post_init__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__post_init__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>"}
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
