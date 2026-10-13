selector_to_html = {"a[href=\"#pyskyfire-pump-diffuser-diffusergeometry\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.pump.diffuser.DiffuserGeometry<a class=\"headerlink\" href=\"#pyskyfire-pump-diffuser-diffusergeometry\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.pump.diffuser.DiffuserGeometry\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.pump.diffuser.DiffuserGeometry\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.pump.diffuser.</span></span><span class=\"sig-name descname\"><span class=\"pre\">DiffuserGeometry</span></span></dt><dd><p>One vane and the walls it stands between, in the impeller\u2019s own frame.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">vane</span></code> is a <code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">BladeSurface</span></code> in (x, r, theta)\nabout the shaft (+x), in the wheel\u2019s rotation sense; <code class=\"docutils literal notranslate\"><span class=\"pre\">suction</span></code> is the\nface toward the impeller and <code class=\"docutils literal notranslate\"><span class=\"pre\">pressure</span></code> the one away from it. <code class=\"docutils literal notranslate\"><span class=\"pre\">hub</span></code>\nand <code class=\"docutils literal notranslate\"><span class=\"pre\">shroud</span></code> are the (x, r) wall curves from the impeller outlet to the\ndiffuser outlet, symmetric about the impeller\u2019s outlet mid-plane at x = 0,\nwith the hub at +x as on the impeller.</p></dd>"}
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
