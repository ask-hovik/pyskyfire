selector_to_html = {"a[href=\"#pyskyfire.common.section_geometry.SectionArc.evaluate\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.section_geometry.SectionArc.evaluate\">\n<span class=\"sig-name descname\"><span class=\"pre\">evaluate</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">fraction</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire-common-section-geometry-sectionarc\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.common.section_geometry.SectionArc<a class=\"headerlink\" href=\"#pyskyfire-common-section-geometry-sectionarc\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.common.section_geometry.SectionArc.points\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.section_geometry.SectionArc.points\">\n<span class=\"sig-name descname\"><span class=\"pre\">points</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">index</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">count</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.common.section_geometry.SectionArc.derivative\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.section_geometry.SectionArc.derivative\">\n<span class=\"sig-name descname\"><span class=\"pre\">derivative</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">fraction</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.common.section_geometry.SectionArc\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.section_geometry.SectionArc\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.common.section_geometry.</span></span><span class=\"sig-name descname\"><span class=\"pre\">SectionArc</span></span></dt><dd></dd>"}
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
