selector_to_html = {"a[href=\"#pyskyfire.pump.geometry.ImpellerGeometry\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.pump.geometry.ImpellerGeometry\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.pump.geometry.</span></span><span class=\"sig-name descname\"><span class=\"pre\">ImpellerGeometry</span></span></dt><dd></dd>", "a[href=\"#pyskyfire.pump.geometry.ImpellerGeometry.__post_init__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.pump.geometry.ImpellerGeometry.__post_init__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__post_init__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire-pump-geometry-impellergeometry\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.pump.geometry.ImpellerGeometry<a class=\"headerlink\" href=\"#pyskyfire-pump-geometry-impellergeometry\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.pump.geometry.ImpellerGeometry.export_npz\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.pump.geometry.ImpellerGeometry.export_npz\">\n<span class=\"sig-name descname\"><span class=\"pre\">export_npz</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">path</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Save sampled geometry without serializing callables or pyTurbo objects.</p></dd>", "a[href=\"#pyskyfire.pump.geometry.ImpellerGeometry.rotation_sign\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.pump.geometry.ImpellerGeometry.rotation_sign\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">rotation_sign</span></span><span class=\"property\"><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/functions.html#float\" title=\"(in Python v3.15)\"><span class=\"pre\">float</span></a></span></dt><dd><p><code class=\"docutils literal notranslate\"><span class=\"pre\">+1</span></code> counterclockwise about <code class=\"docutils literal notranslate\"><span class=\"pre\">+x</span></code>, <code class=\"docutils literal notranslate\"><span class=\"pre\">-1</span></code> clockwise.</p><p>Blade surfaces are already built in this sense. It is the pitch offset\nbetween blades that still needs it, since stepping round the wheel has\nto step the same way the wheel is mirrored.</p></dd>"}
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
