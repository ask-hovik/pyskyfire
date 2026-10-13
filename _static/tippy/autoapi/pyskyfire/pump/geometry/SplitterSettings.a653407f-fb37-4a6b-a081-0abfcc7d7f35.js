selector_to_html = {"a[href=\"#pyskyfire.pump.geometry.SplitterSettings.__post_init__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.pump.geometry.SplitterSettings.__post_init__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__post_init__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire-pump-geometry-splittersettings\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.pump.geometry.SplitterSettings<a class=\"headerlink\" href=\"#pyskyfire-pump-geometry-splittersettings\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.pump.geometry.SplitterSettings\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.pump.geometry.SplitterSettings\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.pump.geometry.</span></span><span class=\"sig-name descname\"><span class=\"pre\">SplitterSettings</span></span></dt><dd><p>One repeated splitter family following the main-blade camber.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">start</span></code> is the fraction of the main camber parameter where it begins;\n<code class=\"docutils literal notranslate\"><span class=\"pre\">pitch_fraction</span></code> locates it between adjacent main blades.</p></dd>"}
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
