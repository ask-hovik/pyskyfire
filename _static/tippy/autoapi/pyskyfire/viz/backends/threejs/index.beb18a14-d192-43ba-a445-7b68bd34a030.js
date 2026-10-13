selector_to_html = {"a[href=\"#pyskyfire.viz.backends.threejs.to_html\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.backends.threejs.to_html\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.backends.threejs.</span></span><span class=\"sig-name descname\"><span class=\"pre\">to_html</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">scene</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#module-contents\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Module Contents<a class=\"headerlink\" href=\"#module-contents\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"#pyskyfire.viz.backends.threejs.scene_payload\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.backends.threejs.scene_payload\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.backends.threejs.</span></span><span class=\"sig-name descname\"><span class=\"pre\">scene_payload</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">scene</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#module-pyskyfire.viz.backends.threejs\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.backends.threejs<a class=\"headerlink\" href=\"#module-pyskyfire.viz.backends.threejs\" title=\"Link to this heading\">\u00b6</a></h1><p>Offline HTML rendering with mesh references and original scalar values.</p><p>The scene is serialized directly instead of exporting a VTK rendering window\nthrough glTF, so fields and node identities survive the browser boundary.</p>", "a[href=\"#functions\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Functions<a class=\"headerlink\" href=\"#functions\" title=\"Link to this heading\">\u00b6</a></h2>"}
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
