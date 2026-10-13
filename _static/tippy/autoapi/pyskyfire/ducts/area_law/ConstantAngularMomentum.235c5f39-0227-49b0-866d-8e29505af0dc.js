selector_to_html = {"a[href=\"#pyskyfire.ducts.area_law.ConstantAngularMomentum.residual\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.area_law.ConstantAngularMomentum.residual\">\n<span class=\"sig-name descname\"><span class=\"pre\">residual</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">theta</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">area</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">centroid_radius</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">inlet_width</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">inlet</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Return the signed miss at one azimuth.</p><p>Negative means the trial section is too small. Must increase with\nsection size so that the solver can bracket a root.</p></dd>", "a[href=\"#pyskyfire.ducts.area_law.ConstantAngularMomentum.__post_init__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.area_law.ConstantAngularMomentum.__post_init__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__post_init__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.ducts.area_law.ConstantAngularMomentum\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.area_law.ConstantAngularMomentum\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.ducts.area_law.</span></span><span class=\"sig-name descname\"><span class=\"pre\">ConstantAngularMomentum</span></span></dt><dd><p>Bases: <code class=\"xref py py-obj docutils literal notranslate\"><span class=\"pre\">AreaLaw</span></code></p><p>Size the passage to conserve the angular momentum it collects.</p><p>Imposing <code class=\"docutils literal notranslate\"><span class=\"pre\">(r*C_u)_c</span> <span class=\"pre\">=</span> <span class=\"pre\">(r*C_u)_5</span></code> together with mass conservation gives\nAungier\u2019s Eqs. 11-1 to 11-3, generalised by the sizing parameter in\nEq. 11-5:</p><p>so the area divided by the centroid radius grows linearly with azimuth.\nBecause both <code class=\"docutils literal notranslate\"><span class=\"pre\">A</span></code> and <code class=\"docutils literal notranslate\"><span class=\"pre\">r_c</span></code> follow from the section size, this is\nimplicit in the size and is solved rather than evaluated.</p></dd>", "a[href=\"#pyskyfire-ducts-area-law-constantangularmomentum\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.ducts.area_law.ConstantAngularMomentum<a class=\"headerlink\" href=\"#pyskyfire-ducts-area-law-constantangularmomentum\" title=\"Link to this heading\">\u00b6</a></h1>"}
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
