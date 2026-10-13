selector_to_html = {"a[href=\"#module-pyskyfire.turbine.duty\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.turbine.duty<a class=\"headerlink\" href=\"#module-pyskyfire.turbine.duty\" title=\"Link to this heading\">\u00b6</a></h1><p>The duty a turbine is designed for: a gas expanding between two pressures.</p><p>Unlike a pump\u2019s, a turbine\u2019s duty cannot be given as head and flow alone.\nThe gas expands on its way through, so the volumetric flow at the exit can be\na good deal larger than at the inlet, and every blade height downstream is\nsized on that. <a class=\"reference internal\" href=\"TurbineDuty.html#pyskyfire.turbine.duty.TurbineDuty\" title=\"pyskyfire.turbine.duty.TurbineDuty\"><code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">TurbineDuty</span></code></a> therefore always carries the states: a\nmedium, a mass flow, the inlet state and the pressure the gas leaves at.</p>", "a[href=\"#classes\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Classes<a class=\"headerlink\" href=\"#classes\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"TurbineDuty.html#pyskyfire.turbine.duty.TurbineDuty\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.turbine.duty.TurbineDuty\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.turbine.duty.</span></span><span class=\"sig-name descname\"><span class=\"pre\">TurbineDuty</span></span></dt><dd><p>Bases: <a class=\"reference internal\" href=\"../../thermo/duty/MachineDuty.html#pyskyfire.thermo.duty.MachineDuty\" title=\"pyskyfire.thermo.duty.MachineDuty\"><code class=\"xref py py-obj docutils literal notranslate\"><span class=\"pre\">pyskyfire.thermo.duty.MachineDuty</span></code></a></p><p>What a turbine is asked to do. See <code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">MachineDuty</span></code>.</p><p>Build one with <a class=\"reference internal\" href=\"#pyskyfire.turbine.duty.TurbineDuty.from_states\" title=\"pyskyfire.turbine.duty.TurbineDuty.from_states\"><code class=\"xref py py-meth docutils literal notranslate\"><span class=\"pre\">from_states()</span></code></a> or <a class=\"reference internal\" href=\"#pyskyfire.turbine.duty.TurbineDuty.from_stations\" title=\"pyskyfire.turbine.duty.TurbineDuty.from_stations\"><code class=\"xref py py-meth docutils literal notranslate\"><span class=\"pre\">from_stations()</span></code></a>.</p></dd>"}
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
