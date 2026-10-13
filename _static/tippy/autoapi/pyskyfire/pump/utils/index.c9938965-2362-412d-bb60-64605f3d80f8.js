selector_to_html = {"a[href=\"#module-pyskyfire.pump.utils\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.pump.utils<a class=\"headerlink\" href=\"#module-pyskyfire.pump.utils\" title=\"Link to this heading\">\u00b6</a></h1><h2>Functions<a class=\"headerlink\" href=\"#functions\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"#module-contents\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Module Contents<a class=\"headerlink\" href=\"#module-contents\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"#pyskyfire.pump.utils.tangential_velocity\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.pump.utils.tangential_velocity\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.pump.utils.</span></span><span class=\"sig-name descname\"><span class=\"pre\">tangential_velocity</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">n</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">d</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Convert rpm and diameter to tangential velocity</p></dd>", "a[href=\"#pyskyfire.pump.utils.interpolate_curve\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.pump.utils.interpolate_curve\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.pump.utils.</span></span><span class=\"sig-name descname\"><span class=\"pre\">interpolate_curve</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">points</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">num_points</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">200</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Given a list of points (tuples, e.g., [(x1, y1), (x2, y2), \u2026]) along a curve,\ninterpolate the curve to produce a new list of points that are equally spaced\nin arc-length.</p></dd>", "a[href=\"#functions\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Functions<a class=\"headerlink\" href=\"#functions\" title=\"Link to this heading\">\u00b6</a></h2>"}
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
