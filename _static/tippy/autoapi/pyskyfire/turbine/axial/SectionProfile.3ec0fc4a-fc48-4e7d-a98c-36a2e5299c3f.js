selector_to_html = {"a[href=\"#pyskyfire.turbine.axial.SectionProfile\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.turbine.axial.SectionProfile\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.turbine.axial.</span></span><span class=\"sig-name descname\"><span class=\"pre\">SectionProfile</span></span></dt><dd><p>The section shape of one kind of row, in fractions of its chord.</p></dd>", "a[href=\"#pyskyfire-turbine-axial-sectionprofile\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.turbine.axial.SectionProfile<a class=\"headerlink\" href=\"#pyskyfire-turbine-axial-sectionprofile\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.turbine.axial.SectionProfile.__post_init__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.turbine.axial.SectionProfile.__post_init__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__post_init__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.turbine.axial.SectionProfile.stagger\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.turbine.axial.SectionProfile.stagger\">\n<span class=\"sig-name descname\"><span class=\"pre\">stagger</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">inlet_metal</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">outlet_metal</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Stagger at mean radius for a section between these metal angles [rad].</p></dd>", "a[href=\"#pyskyfire.turbine.axial.SectionProfile.impulse\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.turbine.axial.SectionProfile.impulse\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">classmethod</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">impulse</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"o\"><span class=\"pre\">**</span></span><span class=\"n\"><span class=\"pre\">changes</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>An impulse blade: symmetric thickness and equally thin edges.</p></dd>"}
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
