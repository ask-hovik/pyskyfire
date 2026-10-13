selector_to_html = {"a[href=\"#pyskyfire.turbine.geometry.TurbineGeometry.__post_init__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.turbine.geometry.TurbineGeometry.__post_init__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__post_init__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.turbine.geometry.TurbineGeometry.row\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.turbine.geometry.TurbineGeometry.row\">\n<span class=\"sig-name descname\"><span class=\"pre\">row</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">name</span></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference internal\" href=\"RowGeometry.html#pyskyfire.turbine.geometry.RowGeometry\" title=\"pyskyfire.turbine.geometry.RowGeometry\"><span class=\"pre\">RowGeometry</span></a></span></span></dt><dd></dd>", "a[href=\"#pyskyfire-turbine-geometry-turbinegeometry\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.turbine.geometry.TurbineGeometry<a class=\"headerlink\" href=\"#pyskyfire-turbine-geometry-turbinegeometry\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"RowGeometry.html#pyskyfire.turbine.geometry.RowGeometry\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.turbine.geometry.RowGeometry\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.turbine.geometry.</span></span><span class=\"sig-name descname\"><span class=\"pre\">RowGeometry</span></span></dt><dd><p>One row\u2019s blade and its measured section properties.</p></dd>", "a[href=\"#pyskyfire.turbine.geometry.TurbineGeometry\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.turbine.geometry.TurbineGeometry\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.turbine.geometry.</span></span><span class=\"sig-name descname\"><span class=\"pre\">TurbineGeometry</span></span></dt><dd><p>Every row\u2019s built blade, and the passage they stand in.</p></dd>", "a[href=\"#pyskyfire.turbine.geometry.TurbineGeometry.rotation_sign\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.turbine.geometry.TurbineGeometry.rotation_sign\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">rotation_sign</span></span><span class=\"property\"><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/functions.html#float\" title=\"(in Python v3.15)\"><span class=\"pre\">float</span></a></span></dt><dd></dd>"}
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
