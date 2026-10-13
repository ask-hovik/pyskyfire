selector_to_html = {"a[href=\"#pyskyfire-viz-components-volutecomponent\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.components.VoluteComponent<a class=\"headerlink\" href=\"#pyskyfire-viz-components-volutecomponent\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.viz.components.VoluteComponent.build\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.components.VoluteComponent.build\">\n<span class=\"sig-name descname\"><span class=\"pre\">build</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.viz.components.VoluteComponent._wall_fields\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.components.VoluteComponent._wall_fields\">\n<span class=\"sig-name descname\"><span class=\"pre\">_wall_fields</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">stations</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">mouth</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">exit_stations</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Field values laid out to match the merged wall mesh\u2019s vertices.</p></dd>", "a[href=\"#pyskyfire.viz.components.VoluteComponent\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.components.VoluteComponent\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.components.</span></span><span class=\"sig-name descname\"><span class=\"pre\">VoluteComponent</span></span></dt><dd><p>Display a resolved volute or collector passage.</p><p>The volute supplies its own geometry through <code class=\"docutils literal notranslate\"><span class=\"pre\">sweep(n_points=...)</span></code>, which\nplaces each station\u2019s cross-section in machine coordinates, already\ninclined by the volute\u2019s tilt. The local frame is that frame: the machine\naxis is +z, the tongue sits at azimuth zero, and no recentering is applied.</p><p>The passage is swept from the tongue through one revolution, and on through\nits <code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">ExitCone</span></code> duct to the flange when it has one.\nThe face at the tongue is a wall \u2013 the cutwater \u2013 so the step visible\nthere is the passage meeting itself.</p><p>Three child nodes come back, so the openings can be told apart from the\nwalls at a glance:</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">opening_color</span></code> paints both openings, red by default. They carry no\nfields, so they keep that colour whichever field the scene is coloured by,\nrather than disappearing into the scale.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">results</span></code> is a <code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">VoluteResult</span></code>, or a mapping of\nfield names to point-associated FieldData with one value per azimuthal\nstation. A supplied mapping holds its last value along the exit cone, which\na result does only for quantities that have no variation to carry there.</p></dd>"}
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
