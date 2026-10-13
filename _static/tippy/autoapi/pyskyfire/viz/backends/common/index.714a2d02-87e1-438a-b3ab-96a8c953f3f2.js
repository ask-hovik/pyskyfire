selector_to_html = {"a[href=\"#module-pyskyfire.viz.backends.common\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.backends.common<a class=\"headerlink\" href=\"#module-pyskyfire.viz.backends.common\" title=\"Link to this heading\">\u00b6</a></h1><p>Scalar display validation shared by desktop and HTML backends.</p>", "a[href=\"#module-contents\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Module Contents<a class=\"headerlink\" href=\"#module-contents\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"#pyskyfire.viz.backends.common.is_visible\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.backends.common.is_visible\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.backends.common.</span></span><span class=\"sig-name descname\"><span class=\"pre\">is_visible</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">node</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.viz.backends.common.field_styles\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.backends.common.field_styles\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.backends.common.</span></span><span class=\"sig-name descname\"><span class=\"pre\">field_styles</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">scene</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Color ranges for every field the HTML viewer\u2019s field selector can offer.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">clim</span></code> pins the range of the initially selected field only; the rest are\nranged over their own visible values. A field whose units disagree between\ncomponents is omitted rather than offered with a meaningless shared scale,\nexcept when it is <code class=\"docutils literal notranslate\"><span class=\"pre\">color_by</span></code>, where it still raises.</p></dd>", "a[href=\"#pyskyfire.viz.backends.common._style_for\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.backends.common._style_for\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.backends.common.</span></span><span class=\"sig-name descname\"><span class=\"pre\">_style_for</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">scene</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">name</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">clim</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.viz.backends.common.field_style\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.backends.common.field_style\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.backends.common.</span></span><span class=\"sig-name descname\"><span class=\"pre\">field_style</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">scene</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#functions\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Functions<a class=\"headerlink\" href=\"#functions\" title=\"Link to this heading\">\u00b6</a></h2>"}
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
