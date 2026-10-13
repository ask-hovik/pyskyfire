selector_to_html = {"a[href=\"#pyskyfire.ducts.manifold.ManifoldGeometry\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.manifold.ManifoldGeometry\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.ducts.manifold.</span></span><span class=\"sig-name descname\"><span class=\"pre\">ManifoldGeometry</span></span></dt><dd><p>Resolved passage geometry along one branch.</p></dd>", "a[href=\"#pyskyfire.ducts.manifold.ManifoldGeometry.dead_clip_angle\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.manifold.ManifoldGeometry.dead_clip_angle\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">dead_clip_angle</span></span><span class=\"property\"><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/functions.html#float\" title=\"(in Python v3.15)\"><span class=\"pre\">float</span></a></span></dt><dd><p>Angle up to which the law could not be satisfied [rad].</p><p>Zero when it is satisfied everywhere. A finite value means the passage\nnear the dead point is larger than the law wants, because the section\nfamily cannot shrink any further \u2013 the manifold\u2019s version of a\nvolute\u2019s tongue.</p></dd>", "a[href=\"#pyskyfire-ducts-manifold-manifoldgeometry\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.ducts.manifold.ManifoldGeometry<a class=\"headerlink\" href=\"#pyskyfire-ducts-manifold-manifoldgeometry\" title=\"Link to this heading\">\u00b6</a></h1>"}
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
