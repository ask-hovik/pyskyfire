selector_to_html = {"a[href=\"#pyskyfire.viz.components.PipeComponent.build\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.components.PipeComponent.build\">\n<span class=\"sig-name descname\"><span class=\"pre\">build</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.viz.components.PipeComponent\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.components.PipeComponent\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.components.</span></span><span class=\"sig-name descname\"><span class=\"pre\">PipeComponent</span></span></dt><dd><p>Display an explicitly defined pipe; never route or modify a solver model.</p><p>The pipe supplies its own geometry. A <code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">Pipe</span></code> offers\n<code class=\"docutils literal notranslate\"><span class=\"pre\">display_path(segments=...)</span></code>, which opens its rounded and continuous bend\nconnections into arcs for drawing without touching the solved model;\n<code class=\"docutils literal notranslate\"><span class=\"pre\">segments</span></code> sets the drawn sides per quarter turn. Any other object is read through <code class=\"docutils literal notranslate\"><span class=\"pre\">centerline</span></code>, an (N,3)\nCartesian polyline in metres, and <code class=\"docutils literal notranslate\"><span class=\"pre\">radius</span></code>, a positive scalar or one\nradius per station, and is drawn exactly as supplied. The local frame is\nthe frame of those coordinates (no implicit recentering).</p><p>A pipe with a core is drawn as a <code class=\"docutils literal notranslate\"><span class=\"pre\">\"core\"</span></code> surface, one per stretch the\ncore runs along, inside a translucent <code class=\"docutils literal notranslate\"><span class=\"pre\">\"bore\"</span></code> of <code class=\"docutils literal notranslate\"><span class=\"pre\">opacity</span></code>, both\nchildren of the returned node, which carries the bore\u2019s end anchors.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">results</span></code> is a <code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">DuctResult</span></code>, or a mapping of\nfield names to FieldData: point-associated with one value per station, or\ncell-associated with one value per connection. Station values are\ninterpolated along refined bends and repeated around each ring; connection\nvalues fill the rings spanning that connection. Every surface carries the\nsame fields.</p></dd>", "a[href=\"#pyskyfire-viz-components-pipecomponent\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.components.PipeComponent<a class=\"headerlink\" href=\"#pyskyfire-viz-components-pipecomponent\" title=\"Link to this heading\">\u00b6</a></h1>"}
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
