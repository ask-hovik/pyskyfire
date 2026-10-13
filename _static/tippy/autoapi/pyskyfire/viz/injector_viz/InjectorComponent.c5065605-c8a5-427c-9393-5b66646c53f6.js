selector_to_html = {"a[href=\"#pyskyfire.viz.injector_viz.InjectorComponent.build\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.injector_viz.InjectorComponent.build\">\n<span class=\"sig-name descname\"><span class=\"pre\">build</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.viz.injector_viz.InjectorComponent\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.injector_viz.InjectorComponent\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.injector_viz.</span></span><span class=\"sig-name descname\"><span class=\"pre\">InjectorComponent</span></span></dt><dd><p>Display the wetted bores of an injector\u2019s holes, coloured by feed.</p><p class=\"rubric\">Notes</p><p>Only wetted surface is drawn as hardware. A hole\u2019s bore runs back from the\nface along its axis for the orifice\u2019s <code class=\"docutils literal notranslate\"><span class=\"pre\">length</span></code>, or two diameters where\nno length is given, and ends where it meets the face plane. The arrows\nare drawn in grey, without fields, as annotation. The <code class=\"docutils literal notranslate\"><span class=\"pre\">face</span></code> anchor sits\non the axis at the face, pointing into the chamber.</p></dd>", "a[href=\"#pyskyfire-viz-injector-viz-injectorcomponent\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.injector_viz.InjectorComponent<a class=\"headerlink\" href=\"#pyskyfire-viz-injector-viz-injectorcomponent\" title=\"Link to this heading\">\u00b6</a></h1>"}
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
