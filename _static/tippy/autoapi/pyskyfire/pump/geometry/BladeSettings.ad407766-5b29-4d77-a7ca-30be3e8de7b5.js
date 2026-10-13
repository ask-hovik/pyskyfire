selector_to_html = {"a[href=\"#pyskyfire.pump.geometry.BladeSettings\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.pump.geometry.BladeSettings\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.pump.geometry.</span></span><span class=\"sig-name descname\"><span class=\"pre\">BladeSettings</span></span></dt><dd><p>Blade extent in normalized hub arc length and edge construction controls.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">start</span></code> is the leading edge. One station puts it on an existing span line;\na <code class=\"docutils literal notranslate\"><span class=\"pre\">(hub,</span> <span class=\"pre\">tip)</span></code> pair slants it, beginning the blade further downstream on\nthe hub as a blade development normally is drawn, and the passage is\nre-paired to carry that line. Both stations are always available as\n<code class=\"docutils literal notranslate\"><span class=\"pre\">start_hub</span></code> and <code class=\"docutils literal notranslate\"><span class=\"pre\">start_tip</span></code>.</p><p>Thickness is a body thickness; both edges close to a point. The leading\nedge is rounded by pyTurbo. <code class=\"docutils literal notranslate\"><span class=\"pre\">trailing_edge</span></code> selects its closure style.\n<code class=\"docutils literal notranslate\"><span class=\"pre\">angle_law</span></code> is used only to estimate total wrap.</p></dd>", "a[href=\"#pyskyfire.pump.geometry.BladeSettings.__post_init__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.pump.geometry.BladeSettings.__post_init__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__post_init__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.pump.geometry.BladeSettings.start_hub\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.pump.geometry.BladeSettings.start_hub\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">start_hub</span></span></dt><dd></dd>", "a[href=\"#pyskyfire-pump-geometry-bladesettings\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.pump.geometry.BladeSettings<a class=\"headerlink\" href=\"#pyskyfire-pump-geometry-bladesettings\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.pump.geometry.BladeSettings.start_tip\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.pump.geometry.BladeSettings.start_tip\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">start_tip</span></span></dt><dd></dd>"}
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
