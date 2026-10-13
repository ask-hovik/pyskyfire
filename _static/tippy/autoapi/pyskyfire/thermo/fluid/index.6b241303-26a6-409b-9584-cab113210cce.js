selector_to_html = {"a[href=\"Fluid.html#pyskyfire.thermo.fluid.Fluid\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.thermo.fluid.Fluid\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.thermo.fluid.</span></span><span class=\"sig-name descname\"><span class=\"pre\">Fluid</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">type</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">propellants</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">fractions</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">basis</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">'mass'</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">precision</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">3</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Mixture helper for fluids.</p><p>Stores a set of component names and composition.</p><p class=\"rubric\">Notes</p><p>Fractions are re-normalized internally when converting to mole basis to\nguard against small rounding errors.</p></dd>", "a[href=\"#classes\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Classes<a class=\"headerlink\" href=\"#classes\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"#module-pyskyfire.thermo.fluid\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.thermo.fluid<a class=\"headerlink\" href=\"#module-pyskyfire.thermo.fluid\" title=\"Link to this heading\">\u00b6</a></h1><h2>Classes<a class=\"headerlink\" href=\"#classes\" title=\"Link to this heading\">\u00b6</a></h2>"}
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
