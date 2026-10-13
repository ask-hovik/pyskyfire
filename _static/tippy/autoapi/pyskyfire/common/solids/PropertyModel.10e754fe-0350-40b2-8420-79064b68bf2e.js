selector_to_html = {"a[href=\"index.html#pyskyfire.common.solids.ArrayLike\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.solids.ArrayLike\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.common.solids.</span></span><span class=\"sig-name descname\"><span class=\"pre\">ArrayLike</span></span></dt><dd></dd>", "a[href=\"#pyskyfire.common.solids.PropertyModel.__call__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.solids.PropertyModel.__call__\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">abstractmethod</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">__call__</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">T</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference internal\" href=\"index.html#pyskyfire.common.solids.ArrayLike\" title=\"pyskyfire.common.solids.ArrayLike\"><span class=\"pre\">ArrayLike</span></a></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray\" title=\"(in NumPy v2.5)\"><span class=\"pre\">ndarray</span></a></span></span></dt><dd></dd>", "a[href=\"#pyskyfire.common.solids.PropertyModel\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.solids.PropertyModel\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.common.solids.</span></span><span class=\"sig-name descname\"><span class=\"pre\">PropertyModel</span></span></dt><dd><p>Abstract base class for temperature-dependent property models.</p><p>All subclasses must implement <code class=\"docutils literal notranslate\"><span class=\"pre\">__call__(T)</span> <span class=\"pre\">-&gt;</span> <span class=\"pre\">np.ndarray</span></code>.</p></dd>", "a[href=\"#pyskyfire-common-solids-propertymodel\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.common.solids.PropertyModel<a class=\"headerlink\" href=\"#pyskyfire-common-solids-propertymodel\" title=\"Link to this heading\">\u00b6</a></h1>"}
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
