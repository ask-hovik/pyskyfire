selector_to_html = {"a[href=\"#pyskyfire.viz.scene.SceneNode\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.scene.SceneNode\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.scene.</span></span><span class=\"sig-name descname\"><span class=\"pre\">SceneNode</span></span></dt><dd><p>A selectable part or group, independent of the renderer\u2019s actor layout.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">tags</span></code> are free-form labels. The HTML viewer offers every tag as a group\nof the nodes carrying it, as it does every <code class=\"docutils literal notranslate\"><span class=\"pre\">metadata[\"fluid\"]</span></code> a\ncomponent recorded from its results.</p></dd>", "a[href=\"#pyskyfire-viz-scene-scenenode\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.scene.SceneNode<a class=\"headerlink\" href=\"#pyskyfire-viz-scene-scenenode\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.viz.scene.SceneNode.id\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.scene.SceneNode.id\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">id</span></span></dt><dd></dd>", "a[href=\"#pyskyfire.viz.scene.SceneNode.world_transform\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.scene.SceneNode.world_transform\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">world_transform</span></span></dt><dd></dd>", "a[href=\"#pyskyfire.viz.scene.SceneNode.anchor\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.scene.SceneNode.anchor\">\n<span class=\"sig-name descname\"><span class=\"pre\">anchor</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">name</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Return a snapshot of a named anchor in scene coordinates.</p></dd>"}
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
