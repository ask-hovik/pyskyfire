selector_to_html = {"a[href=\"index.html#pyskyfire.common.solids.ArrayLike\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.solids.ArrayLike\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.common.solids.</span></span><span class=\"sig-name descname\"><span class=\"pre\">ArrayLike</span></span></dt><dd></dd>", "a[href=\"#pyskyfire.common.solids.Log10PolynomialModel\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.solids.Log10PolynomialModel\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.common.solids.</span></span><span class=\"sig-name descname\"><span class=\"pre\">Log10PolynomialModel</span></span></dt><dd><p>Bases: <code class=\"xref py py-obj docutils literal notranslate\"><span class=\"pre\">PropertyModel</span></code></p><p>Log-polynomial model (NIST cryogenic form).</p><p>Implements <span class=\"math notranslate nohighlight\">\\(y = 10^{P(\\log_{10} T)}\\)</span>,\nwhere <span class=\"math notranslate nohighlight\">\\(P(x) = \\sum_i a_i x^i\\)</span>.</p><p class=\"rubric\">Notes</p><p>Useful for fitting data that spans multiple orders of magnitude in\ntemperature or property values.</p></dd>", "a[href=\"#pyskyfire-common-solids-log10polynomialmodel\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.common.solids.Log10PolynomialModel<a class=\"headerlink\" href=\"#pyskyfire-common-solids-log10polynomialmodel\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.common.solids.Log10PolynomialModel.__call__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.solids.Log10PolynomialModel.__call__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__call__</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">T</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference internal\" href=\"index.html#pyskyfire.common.solids.ArrayLike\" title=\"pyskyfire.common.solids.ArrayLike\"><span class=\"pre\">ArrayLike</span></a></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray\" title=\"(in NumPy v2.5)\"><span class=\"pre\">ndarray</span></a></span></span></dt><dd><p>Evaluate the log-polynomial model.</p></dd>"}
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
