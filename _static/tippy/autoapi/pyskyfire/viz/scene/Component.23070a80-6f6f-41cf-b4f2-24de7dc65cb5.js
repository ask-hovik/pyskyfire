selector_to_html = {"a[href=\"#pyskyfire-viz-scene-component\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.scene.Component<a class=\"headerlink\" href=\"#pyskyfire-viz-scene-component\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.viz.scene.Component\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.scene.Component\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.scene.</span></span><span class=\"sig-name descname\"><span class=\"pre\">Component</span></span></dt><dd><p>Bases: <code class=\"xref py py-obj docutils literal notranslate\"><span class=\"pre\">Protocol</span></code></p><p>Base class for protocol classes.</p><p>Protocol classes are defined as:</p><p>Such classes are primarily used with static type checkers that recognize\nstructural subtyping (static duck-typing).</p><p>For example:</p><p>See PEP 544 for details. Protocol classes decorated with\n@typing.runtime_checkable act as simple-minded runtime protocols that check\nonly the presence of given attributes, ignoring their type signatures.\nProtocol classes can be generic, they are defined as:</p></dd>", "a[href=\"SceneNode.html#pyskyfire.viz.scene.SceneNode\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.scene.SceneNode\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.scene.</span></span><span class=\"sig-name descname\"><span class=\"pre\">SceneNode</span></span></dt><dd><p>A selectable part or group, independent of the renderer\u2019s actor layout.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">tags</span></code> are free-form labels. The HTML viewer offers every tag as a group\nof the nodes carrying it, as it does every <code class=\"docutils literal notranslate\"><span class=\"pre\">metadata[\"fluid\"]</span></code> a\ncomponent recorded from its results.</p></dd>", "a[href=\"#pyskyfire.viz.scene.Component.build\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.scene.Component.build\">\n<span class=\"sig-name descname\"><span class=\"pre\">build</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference internal\" href=\"SceneNode.html#pyskyfire.viz.scene.SceneNode\" title=\"pyskyfire.viz.scene.SceneNode\"><span class=\"pre\">SceneNode</span></a></span></span></dt><dd><p>Return a fresh hierarchy, optionally sharing immutable mesh data.</p></dd>"}
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
