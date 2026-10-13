selector_to_html = {"a[href=\"#pyskyfire.turbine.definition.ThicknessDistribution.__post_init__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.turbine.definition.ThicknessDistribution.__post_init__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__post_init__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.turbine.definition.ThicknessDistribution.symmetric\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.turbine.definition.ThicknessDistribution.symmetric\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">classmethod</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">symmetric</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd><p>Thickest at mid-chord, thin at both edges: an impulse blade.</p></dd>", "a[href=\"#pyskyfire.turbine.definition.ThicknessDistribution\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.turbine.definition.ThicknessDistribution\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.turbine.definition.</span></span><span class=\"sig-name descname\"><span class=\"pre\">ThicknessDistribution</span></span></dt><dd><p>How a section\u2019s thickness is laid along its camber line.</p><p>The section is drawn as offsets from the camber line: control points on\nthe suction side and on the pressure side, each a fraction of the row\u2019s\n<code class=\"docutils literal notranslate\"><span class=\"pre\">max_thickness</span></code>, spaced from leading edge towards trailing edge. The\nbuilt section\u2019s actual thickness is measured, not assumed from these.</p></dd>", "a[href=\"#pyskyfire-turbine-definition-thicknessdistribution\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.turbine.definition.ThicknessDistribution<a class=\"headerlink\" href=\"#pyskyfire-turbine-definition-thicknessdistribution\" title=\"Link to this heading\">\u00b6</a></h1>"}
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
