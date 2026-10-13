selector_to_html = {"a[href=\"#pyskyfire.common.blade_surface.BladeSurface.__post_init__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.blade_surface.BladeSurface.__post_init__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__post_init__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire-common-blade-surface-bladesurface\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.common.blade_surface.BladeSurface<a class=\"headerlink\" href=\"#pyskyfire-common-blade-surface-bladesurface\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.common.blade_surface.BladeSurface\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.blade_surface.BladeSurface\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.common.blade_surface.</span></span><span class=\"sig-name descname\"><span class=\"pre\">BladeSurface</span></span></dt><dd><p>One blade\u2019s closed-edge sides, (span, chord, 3), in (x,r,theta).</p></dd>"}
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
