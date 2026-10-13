selector_to_html = {"a[href=\"#pyskyfire-regen-channel-helpers-sectionplanes\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.regen.channel_helpers.SectionPlanes<a class=\"headerlink\" href=\"#pyskyfire-regen-channel-helpers-sectionplanes\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.regen.channel_helpers.SectionPlanes.__post_init__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.channel_helpers.SectionPlanes.__post_init__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__post_init__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/constants.html#None\" title=\"(in Python v3.15)\"><span class=\"pre\">None</span></a></span></span></dt><dd></dd>", "a[href=\"#pyskyfire.regen.channel_helpers.SectionPlanes\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.channel_helpers.SectionPlanes\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.regen.channel_helpers.</span></span><span class=\"sig-name descname\"><span class=\"pre\">SectionPlanes</span></span></dt><dd><p>Orthonormal planes located at the samples of a channel centerline.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">origins</span></code> and the three axes use the Cartesian coordinate convention of\nPySkyFire: X is the engine axis, while Y and Z span a radial plane.\n<code class=\"docutils literal notranslate\"><span class=\"pre\">normals</span></code> follow either the centerline tangent or the engine axis.\n<code class=\"docutils literal notranslate\"><span class=\"pre\">x_axes</span></code> point outward and form the zero-angle direction of a section.</p></dd>", "a[href=\"#pyskyfire.regen.channel_helpers.SectionPlanes.__len__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.channel_helpers.SectionPlanes.__len__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__len__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/functions.html#int\" title=\"(in Python v3.15)\"><span class=\"pre\">int</span></a></span></span></dt><dd></dd>", "a[href=\"#pyskyfire.regen.channel_helpers.SectionPlanes.y_axes\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.channel_helpers.SectionPlanes.y_axes\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">y_axes</span></span><span class=\"property\"><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><a class=\"reference external\" href=\"https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray\" title=\"(in NumPy v2.5)\"><span class=\"pre\">ndarray</span></a></span></dt><dd><p>Return each plane\u2019s second in-plane axis.</p></dd>"}
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
