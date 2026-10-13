selector_to_html = {"a[href=\"index.html#pyskyfire.common.solids.ArrayLike\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.solids.ArrayLike\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.common.solids.</span></span><span class=\"sig-name descname\"><span class=\"pre\">ArrayLike</span></span></dt><dd></dd>", "a[href=\"#pyskyfire.common.solids.TabulatedModel\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.solids.TabulatedModel\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.common.solids.</span></span><span class=\"sig-name descname\"><span class=\"pre\">TabulatedModel</span></span></dt><dd><p>Bases: <code class=\"xref py py-obj docutils literal notranslate\"><span class=\"pre\">PropertyModel</span></code></p><p>Property model based on tabulated data.</p><p>Performs 1D linear interpolation between discrete (T, Y) points.</p></dd>", "a[href=\"#pyskyfire.common.solids.TabulatedModel.__post_init__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.solids.TabulatedModel.__post_init__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__post_init__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire-common-solids-tabulatedmodel\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.common.solids.TabulatedModel<a class=\"headerlink\" href=\"#pyskyfire-common-solids-tabulatedmodel\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.common.solids.TabulatedModel.__call__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.solids.TabulatedModel.__call__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__call__</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">T</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference internal\" href=\"index.html#pyskyfire.common.solids.ArrayLike\" title=\"pyskyfire.common.solids.ArrayLike\"><span class=\"pre\">ArrayLike</span></a></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray\" title=\"(in NumPy v2.5)\"><span class=\"pre\">ndarray</span></a></span></span></dt><dd><p>Interpolate the property at temperature <code class=\"docutils literal notranslate\"><span class=\"pre\">T</span></code>.</p></dd>"}
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
