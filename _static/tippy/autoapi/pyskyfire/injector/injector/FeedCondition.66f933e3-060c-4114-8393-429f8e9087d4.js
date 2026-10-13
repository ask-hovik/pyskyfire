selector_to_html = {"a[href=\"#pyskyfire.injector.injector.FeedCondition\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.injector.injector.FeedCondition\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.injector.injector.</span></span><span class=\"sig-name descname\"><span class=\"pre\">FeedCondition</span></span></dt><dd><p>The operating point of one feed.</p><p>Give the plenum <code class=\"docutils literal notranslate\"><span class=\"pre\">temperature</span></code> [K], the <code class=\"docutils literal notranslate\"><span class=\"pre\">fluid</span></code> (anything\n<a class=\"reference internal\" href=\"../../thermo/eos/index.html#pyskyfire.thermo.eos.as_equation_of_state\" title=\"pyskyfire.thermo.eos.as_equation_of_state\"><code class=\"xref py py-func docutils literal notranslate\"><span class=\"pre\">as_equation_of_state()</span></code></a> reads), and either the\nplenum <code class=\"docutils literal notranslate\"><span class=\"pre\">pressure</span></code> [Pa], to find the flow it drives, or the\n<code class=\"docutils literal notranslate\"><span class=\"pre\">mass_flow</span></code> [kg/s], to find the pressure it needs.</p></dd>", "a[href=\"#pyskyfire.injector.injector.FeedCondition.__post_init__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.injector.injector.FeedCondition.__post_init__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__post_init__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"../../thermo/eos/index.html#pyskyfire.thermo.eos.as_equation_of_state\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.thermo.eos.as_equation_of_state\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.thermo.eos.</span></span><span class=\"sig-name descname\"><span class=\"pre\">as_equation_of_state</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">fluid</span></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference internal\" href=\"../../thermo/eos/EquationOfState.html#pyskyfire.thermo.eos.EquationOfState\" title=\"pyskyfire.thermo.eos.EquationOfState\"><span class=\"pre\">EquationOfState</span></a></span></span></dt><dd><p>Read <code class=\"docutils literal notranslate\"><span class=\"pre\">fluid</span></code> as an equation of state.</p><p>An <a class=\"reference internal\" href=\"../../thermo/eos/EquationOfState.html#pyskyfire.thermo.eos.EquationOfState\" title=\"pyskyfire.thermo.eos.EquationOfState\"><code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">EquationOfState</span></code></a> is returned as it is. Anything else names a\nCoolProp fluid and is handed to <a class=\"reference internal\" href=\"../../thermo/eos/CoolPropEOS.html#pyskyfire.thermo.eos.CoolPropEOS\" title=\"pyskyfire.thermo.eos.CoolPropEOS\"><code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">CoolPropEOS</span></code></a>.</p></dd>", "a[href=\"#pyskyfire-injector-injector-feedcondition\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.injector.injector.FeedCondition<a class=\"headerlink\" href=\"#pyskyfire-injector-injector-feedcondition\" title=\"Link to this heading\">\u00b6</a></h1>"}
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
