selector_to_html = {"a[href=\"#pyskyfire.viz.components.DiffuserComponent\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.components.DiffuserComponent\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.components.</span></span><span class=\"sig-name descname\"><span class=\"pre\">DiffuserComponent</span></span></dt><dd><p>Display a built vaned diffuser: its two walls and every vane.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">diffuser.build_geometry()</span></code> must have been called. The diffuser is drawn\nin the impeller\u2019s own frame, walls symmetric about the impeller outlet\nmid-plane at local x = 0, so it takes the impeller\u2019s rotation and a\ntranslation to that mid-plane. The <code class=\"docutils literal notranslate\"><span class=\"pre\">outlet</span></code> anchor is not a port: the\ndiffuser discharges radially all round into its volute.</p></dd>", "a[href=\"#pyskyfire-viz-components-diffusercomponent\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.components.DiffuserComponent<a class=\"headerlink\" href=\"#pyskyfire-viz-components-diffusercomponent\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.viz.components.DiffuserComponent.build\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.components.DiffuserComponent.build\">\n<span class=\"sig-name descname\"><span class=\"pre\">build</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>"}
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
