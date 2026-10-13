selector_to_html = {"a[href=\"#classes\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Classes<a class=\"headerlink\" href=\"#classes\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"MachineDuty.html#pyskyfire.thermo.duty.MachineDuty\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.thermo.duty.MachineDuty\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.thermo.duty.</span></span><span class=\"sig-name descname\"><span class=\"pre\">MachineDuty</span></span></dt><dd><p>Head, flow and speed asked of one turbomachine.</p><p>Build one through a subclass\u2019s constructors rather than field by field;\nthey are what keep the fields consistent with one another.</p></dd>", "a[href=\"#module-pyskyfire.thermo.duty\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.thermo.duty<a class=\"headerlink\" href=\"#module-pyskyfire.thermo.duty\" title=\"Link to this heading\">\u00b6</a></h1><p>What a turbomachine is asked to do, before anything is decided about how.</p><p>A duty is the fluid\u2019s side of the contract: the mass flow, the state it\narrives in, the pressure it must leave at and the shaft speed. Everything a\ndesign sizes itself on follows from those four, and <a class=\"reference internal\" href=\"MachineDuty.html#pyskyfire.thermo.duty.MachineDuty\" title=\"pyskyfire.thermo.duty.MachineDuty\"><code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">MachineDuty</span></code></a> holds\nthe results in the form designs read:</p>"}
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
