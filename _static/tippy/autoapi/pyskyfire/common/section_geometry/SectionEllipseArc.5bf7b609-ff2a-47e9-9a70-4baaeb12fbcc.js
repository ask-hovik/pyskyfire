selector_to_html = {"a[href=\"#pyskyfire-common-section-geometry-sectionellipsearc\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.common.section_geometry.SectionEllipseArc<a class=\"headerlink\" href=\"#pyskyfire-common-section-geometry-sectionellipsearc\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.common.section_geometry.SectionEllipseArc.evaluate\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.section_geometry.SectionEllipseArc.evaluate\">\n<span class=\"sig-name descname\"><span class=\"pre\">evaluate</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">fraction</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.common.section_geometry.SectionEllipseArc._angle\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.section_geometry.SectionEllipseArc._angle\">\n<span class=\"sig-name descname\"><span class=\"pre\">_angle</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">fraction</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.common.section_geometry.SectionEllipseArc.points\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.section_geometry.SectionEllipseArc.points\">\n<span class=\"sig-name descname\"><span class=\"pre\">points</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">index</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">count</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.common.section_geometry.SectionEllipseArc\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.section_geometry.SectionEllipseArc\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.common.section_geometry.</span></span><span class=\"sig-name descname\"><span class=\"pre\">SectionEllipseArc</span></span></dt><dd><p>Elliptical boundary, with the same sampling API as a circular arc.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">semi_x</span></code> and <code class=\"docutils literal notranslate\"><span class=\"pre\">semi_y</span></code> are the semi-axes along the section plane\u2019s own\ntwo axes, so equal semi-axes reproduce a <code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">SectionArc</span></code> exactly. The\nangles are the ellipse\u2019s parameter, not the polar angle to the point: a\npoint sits at <code class=\"docutils literal notranslate\"><span class=\"pre\">center</span> <span class=\"pre\">+</span> <span class=\"pre\">(semi_x*cos(t),</span> <span class=\"pre\">semi_y*sin(t))</span></code>.</p></dd>", "a[href=\"#pyskyfire.common.section_geometry.SectionEllipseArc.derivative\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.section_geometry.SectionEllipseArc.derivative\">\n<span class=\"sig-name descname\"><span class=\"pre\">derivative</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">fraction</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>"}
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
