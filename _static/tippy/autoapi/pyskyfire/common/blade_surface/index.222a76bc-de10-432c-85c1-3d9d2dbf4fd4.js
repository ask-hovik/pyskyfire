selector_to_html = {"a[href=\"#pyskyfire.common.blade_surface.cylindrical_to_cartesian\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.blade_surface.cylindrical_to_cartesian\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.common.blade_surface.</span></span><span class=\"sig-name descname\"><span class=\"pre\">cylindrical_to_cartesian</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">cyl</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"BladeSurface.html#pyskyfire.common.blade_surface.BladeSurface\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.blade_surface.BladeSurface\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.common.blade_surface.</span></span><span class=\"sig-name descname\"><span class=\"pre\">BladeSurface</span></span></dt><dd><p>One blade\u2019s closed-edge sides, (span, chord, 3), in (x,r,theta).</p></dd>", "a[href=\"#module-pyskyfire.common.blade_surface\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.common.blade_surface<a class=\"headerlink\" href=\"#module-pyskyfire.common.blade_surface\" title=\"Link to this heading\">\u00b6</a></h1><p>Blade surfaces in cylindrical coordinates, shared by impellers and turbines.</p><p>A blade is two sides, suction and pressure, each an array of shape\n<code class=\"docutils literal notranslate\"><span class=\"pre\">(span,</span> <span class=\"pre\">chord,</span> <span class=\"pre\">3)</span></code> in <code class=\"docutils literal notranslate\"><span class=\"pre\">(x,</span> <span class=\"pre\">r,</span> <span class=\"pre\">theta)</span></code>: metres, metres, radians, about the\nshaft (+x). The sides meet at the leading and trailing edges, so their first\nand last chord stations coincide.</p>", "a[href=\"#module-contents\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Module Contents<a class=\"headerlink\" href=\"#module-contents\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"#functions\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Functions<a class=\"headerlink\" href=\"#functions\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"#classes\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Classes<a class=\"headerlink\" href=\"#classes\" title=\"Link to this heading\">\u00b6</a></h2>"}
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
