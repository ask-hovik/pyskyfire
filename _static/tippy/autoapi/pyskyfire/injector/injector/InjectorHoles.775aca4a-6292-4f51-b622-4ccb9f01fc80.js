selector_to_html = {"a[href=\"#pyskyfire.injector.injector.InjectorHoles\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.injector.injector.InjectorHoles\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.injector.injector.</span></span><span class=\"sig-name descname\"><span class=\"pre\">InjectorHoles</span></span></dt><dd><p>Every hole of an injector, placed in the injector frame.</p><p>Arrays have one row per hole. <code class=\"docutils literal notranslate\"><span class=\"pre\">set_name</span></code> and <code class=\"docutils literal notranslate\"><span class=\"pre\">element_index</span></code> say\nwhich element a hole belongs to; <code class=\"docutils literal notranslate\"><span class=\"pre\">feeds</span></code> and <code class=\"docutils literal notranslate\"><span class=\"pre\">orifices</span></code> are tuples\naligned with them. <code class=\"docutils literal notranslate\"><span class=\"pre\">footprint</span></code> is the radius of the hole\u2019s outline on\nthe face [m]: half its diameter, stretched by the inclination.</p></dd>", "a[href=\"#pyskyfire.injector.injector.InjectorHoles.__len__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.injector.injector.InjectorHoles.__len__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__len__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire-injector-injector-injectorholes\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.injector.injector.InjectorHoles<a class=\"headerlink\" href=\"#pyskyfire-injector-injector-injectorholes\" title=\"Link to this heading\">\u00b6</a></h1>"}
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
