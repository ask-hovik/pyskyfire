selector_to_html = {"a[href=\"#pyskyfire-injector-injector-injectorchecks\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.injector.injector.InjectorChecks<a class=\"headerlink\" href=\"#pyskyfire-injector-injector-injectorchecks\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.injector.injector.InjectorChecks\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.injector.injector.InjectorChecks\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.injector.injector.</span></span><span class=\"sig-name descname\"><span class=\"pre\">InjectorChecks</span></span></dt><dd><p>What the injector can check about its own layout, reported not enforced.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">min_web</span></code> is the narrowest ligament between two hole outlines on the\nface [m], between holes in sets <code class=\"docutils literal notranslate\"><span class=\"pre\">min_web_sets</span></code>; <code class=\"docutils literal notranslate\"><span class=\"pre\">overlaps</span></code> counts\npairs whose outlines intersect. <code class=\"docutils literal notranslate\"><span class=\"pre\">outer_radius</span></code> is how far the\noutermost hole outline reaches from the axis, and <code class=\"docutils literal notranslate\"><span class=\"pre\">target_radius</span></code> the\nfurthest any jet\u2019s target lies from it [m]. <code class=\"docutils literal notranslate\"><span class=\"pre\">wall_clearance</span></code> is the\nchamber radius less <code class=\"docutils literal notranslate\"><span class=\"pre\">outer_radius</span></code>: room between the holes and the\nwall. <code class=\"docutils literal notranslate\"><span class=\"pre\">target_clearance</span></code> is the chamber radius less <code class=\"docutils literal notranslate\"><span class=\"pre\">target_radius</span></code>:\nzero where a jet is aimed at the wall on purpose, as a film-cooling jet\nis, and negative where a target lies beyond it. Both are NaN without a\nchamber radius. <code class=\"docutils literal notranslate\"><span class=\"pre\">min_length_to_diameter</span></code> is the shortest hole relative\nto its diameter.</p></dd>"}
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
