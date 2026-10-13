selector_to_html = {"a[href=\"#pyskyfire.viz.scene.Anchor.from_direction\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.scene.Anchor.from_direction\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">classmethod</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">from_direction</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">position</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">direction</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.viz.scene.Anchor\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.scene.Anchor\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.scene.</span></span><span class=\"sig-name descname\"><span class=\"pre\">Anchor</span></span></dt><dd><p>Named coordinate frame. The +x axis is the outward port direction.</p></dd>", "a[href=\"#pyskyfire.viz.scene.Anchor.position\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.scene.Anchor.position\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">position</span></span></dt><dd></dd>", "a[href=\"#pyskyfire-viz-scene-anchor\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.scene.Anchor<a class=\"headerlink\" href=\"#pyskyfire-viz-scene-anchor\" title=\"Link to this heading\">\u00b6</a></h1>"}
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
