selector_to_html = {"a[href=\"#pyskyfire.viz.engine_viz.Engine3DViewer.close\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.engine_viz.Engine3DViewer.close\">\n<span class=\"sig-name descname\"><span class=\"pre\">close</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/constants.html#None\" title=\"(in Python v3.15)\"><span class=\"pre\">None</span></a></span></span></dt><dd><p>Release PyVista rendering resources.</p></dd>", "a[href=\"#pyskyfire.viz.engine_viz.Engine3DViewer.save_html\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.engine_viz.Engine3DViewer.save_html\">\n<span class=\"sig-name descname\"><span class=\"pre\">save_html</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">path</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/stdtypes.html#str\" title=\"(in Python v3.15)\"><span class=\"pre\">str</span></a><span class=\"w\"> </span><span class=\"p\"><span class=\"pre\">|</span></span><span class=\"w\"> </span><a class=\"reference external\" href=\"https://docs.python.org/3/library/pathlib.html#pathlib.Path\" title=\"(in Python v3.15)\"><span class=\"pre\">Path</span></a></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/pathlib.html#pathlib.Path\" title=\"(in Python v3.15)\"><span class=\"pre\">Path</span></a></span></span></dt><dd><p>Write the interactive viewer to a standalone HTML file.</p></dd>", "a[href=\"#pyskyfire.viz.engine_viz.Engine3DViewer\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.engine_viz.Engine3DViewer\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.engine_viz.</span></span><span class=\"sig-name descname\"><span class=\"pre\">Engine3DViewer</span></span></dt><dd><p>Deprecated compatibility handle; new code should use <code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">Scene</span></code>.</p></dd>", "a[href=\"#pyskyfire.viz.engine_viz.Engine3DViewer.__post_init__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.engine_viz.Engine3DViewer.__post_init__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__post_init__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire-viz-engine-viz-engine3dviewer\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.engine_viz.Engine3DViewer<a class=\"headerlink\" href=\"#pyskyfire-viz-engine-viz-engine3dviewer\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.viz.engine_viz.Engine3DViewer.show\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.engine_viz.Engine3DViewer.show\">\n<span class=\"sig-name descname\"><span class=\"pre\">show</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd><p>Open the PyVista desktop viewer.</p></dd>"}
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
