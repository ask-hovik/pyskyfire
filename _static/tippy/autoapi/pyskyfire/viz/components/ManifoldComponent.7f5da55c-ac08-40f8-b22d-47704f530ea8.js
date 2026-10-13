selector_to_html = {"a[href=\"#pyskyfire.viz.components.ManifoldComponent._wall_fields\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.components.ManifoldComponent._wall_fields\">\n<span class=\"sig-name descname\"><span class=\"pre\">_wall_fields</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">stations</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">order</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Field values laid out to match the merged wall mesh\u2019s vertices.</p></dd>", "a[href=\"#pyskyfire.viz.components.ManifoldComponent\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.components.ManifoldComponent\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.components.</span></span><span class=\"sig-name descname\"><span class=\"pre\">ManifoldComponent</span></span></dt><dd><p>Display a resolved distributing or collecting manifold.</p><p>The manifold supplies its own geometry through <code class=\"docutils literal notranslate\"><span class=\"pre\">sweep(n_points=...)</span></code>,\nalready placed on the wall it was built against, and the local frame is\nthat frame: the machine axis is +z and the dead point sits at azimuth zero.\nA split manifold sweeps the whole ring, both branches together.</p><p>Three child nodes come back, so the opening can be told apart from the\nwalls at a glance:</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">opening_color</span></code> paints only those faces, red by default: the places flow\ncrosses the boundary, not the walls that carry it there. They hold no\nfields, so they keep that colour whichever field the scene is coloured by,\nrather than disappearing into the scale.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">results</span></code> is a mapping of field names to point-associated FieldData, one\nvalue per swept station of the resolved ring \u2013 which for a split manifold\nis both branches, not one.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">n_points</span></code> samples each boundary piece of the section. <code class=\"docutils literal notranslate\"><span class=\"pre\">n_stations</span></code>\nrefines the ring in azimuth for drawing, which is the direction a manifold\nfacets in first, and leaves the resolved design alone; values on the\nresolved grid are read between its stations onto the refined mesh. <code class=\"docutils literal notranslate\"><span class=\"pre\">None</span></code>\ndraws the ring exactly as it was resolved.</p></dd>", "a[href=\"#pyskyfire-viz-components-manifoldcomponent\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.components.ManifoldComponent<a class=\"headerlink\" href=\"#pyskyfire-viz-components-manifoldcomponent\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.viz.components.ManifoldComponent.build\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.components.ManifoldComponent.build\">\n<span class=\"sig-name descname\"><span class=\"pre\">build</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>"}
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
