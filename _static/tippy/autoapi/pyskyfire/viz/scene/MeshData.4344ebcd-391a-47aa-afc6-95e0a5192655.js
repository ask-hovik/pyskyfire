selector_to_html = {"a[href=\"#pyskyfire.viz.scene.MeshData.__post_init__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.scene.MeshData.__post_init__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__post_init__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire-viz-scene-meshdata\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.scene.MeshData<a class=\"headerlink\" href=\"#pyskyfire-viz-scene-meshdata\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.viz.scene.MeshData\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.scene.MeshData\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.scene.</span></span><span class=\"sig-name descname\"><span class=\"pre\">MeshData</span></span></dt><dd><p>Surface geometry in local Cartesian coordinates; indices describe triangles.</p></dd>"}
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
