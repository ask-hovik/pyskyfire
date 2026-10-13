selector_to_html = {"a[href=\"#pyskyfire-regen-thrust-chamber-wall\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.regen.thrust_chamber.Wall<a class=\"headerlink\" href=\"#pyskyfire-regen-thrust-chamber-wall\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.regen.thrust_chamber.Wall\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.thrust_chamber.Wall\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.regen.thrust_chamber.</span></span><span class=\"sig-name descname\"><span class=\"pre\">Wall</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">material</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">thickness</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">name</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">None</span></span></em>, <em class=\"sig-param\"><span class=\"keyword-only-separator o\"><abbr title=\"Keyword-only parameters separator (PEP 3102)\"><span class=\"pre\">*</span></abbr></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">role</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">'main'</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Single structural wall or coating layer in a thrust chamber.</p><p>Encapsulates material and thickness information for heat-conduction\ncalculations.</p></dd>", "a[href=\"#pyskyfire.regen.thrust_chamber.Wall.__setstate__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.thrust_chamber.Wall.__setstate__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__setstate__</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">state</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Restore legacy walls with a main role and a thickness evaluator.</p></dd>", "a[href=\"#pyskyfire.regen.thrust_chamber.Wall.thickness\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.thrust_chamber.Wall.thickness\">\n<span class=\"sig-name descname\"><span class=\"pre\">thickness</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">x</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Return the layer thickness [m] at <code class=\"docutils literal notranslate\"><span class=\"pre\">x</span></code>, which may be an array.</p></dd>", "a[href=\"#pyskyfire.regen.thrust_chamber.Wall.bind_contour\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.thrust_chamber.Wall.bind_contour\">\n<span class=\"sig-name descname\"><span class=\"pre\">bind_contour</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">contour</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Attach the contour that tier A design inputs are evaluated against.</p><p>Called by <code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">CoolingCircuit</span></code>. A wall whose thickness only depends\non <code class=\"docutils literal notranslate\"><span class=\"pre\">x</span></code> never needs this.</p></dd>"}
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
