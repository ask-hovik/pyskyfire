selector_to_html = {"a[href=\"#pyskyfire.common.solids.SumOfGaussiansModel\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.solids.SumOfGaussiansModel\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.common.solids.</span></span><span class=\"sig-name descname\"><span class=\"pre\">SumOfGaussiansModel</span></span></dt><dd><p>Bases: <code class=\"xref py py-obj docutils literal notranslate\"><span class=\"pre\">PropertyModel</span></code></p><p>Empirical model as a sum of Gaussian peaks.</p><p>Implements <span class=\"math notranslate nohighlight\">\\(y(T) = \\sum_i b_i \\exp(-((T - c_i)/d_i)^2)\\)</span>.</p></dd>", "a[href=\"#pyskyfire.common.solids.SumOfGaussiansModel.__call__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.solids.SumOfGaussiansModel.__call__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__call__</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">T</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Evaluate the Gaussian sum model.</p></dd>", "a[href=\"#pyskyfire-common-solids-sumofgaussiansmodel\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.common.solids.SumOfGaussiansModel<a class=\"headerlink\" href=\"#pyskyfire-common-solids-sumofgaussiansmodel\" title=\"Link to this heading\">\u00b6</a></h1>"}
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
