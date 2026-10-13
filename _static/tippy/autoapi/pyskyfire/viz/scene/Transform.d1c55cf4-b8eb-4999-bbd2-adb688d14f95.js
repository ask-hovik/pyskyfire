selector_to_html = {"a[href=\"#pyskyfire.viz.scene.Transform\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.scene.Transform\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.scene.</span></span><span class=\"sig-name descname\"><span class=\"pre\">Transform</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">translation</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">None</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">rotation</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">None</span></span></em>, <em class=\"sig-param\"><span class=\"keyword-only-separator o\"><abbr title=\"Keyword-only parameters separator (PEP 3102)\"><span class=\"pre\">*</span></abbr></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">matrix</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">None</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Rigid local-to-parent transform, using column vectors.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">rotation</span></code> is a scipy Rotation or extrinsic xyz Euler angles in radians.\nA supplied <code class=\"docutils literal notranslate\"><span class=\"pre\">matrix</span></code> is mutually exclusive with translation and rotation.</p></dd>", "a[href=\"#pyskyfire-viz-scene-transform\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.scene.Transform<a class=\"headerlink\" href=\"#pyskyfire-viz-scene-transform\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.viz.scene.Transform.__matmul__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.scene.Transform.__matmul__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__matmul__</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">other</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.viz.scene.Transform.apply\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.scene.Transform.apply\">\n<span class=\"sig-name descname\"><span class=\"pre\">apply</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">points</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.viz.scene.Transform.inverse\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.scene.Transform.inverse\">\n<span class=\"sig-name descname\"><span class=\"pre\">inverse</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>"}
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
