selector_to_html = {"a[href=\"#module-pyskyfire.pump.duty\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.pump.duty<a class=\"headerlink\" href=\"#module-pyskyfire.pump.duty\" title=\"Link to this heading\">\u00b6</a></h1><p>The duty a pump is designed for: flow, head and speed.</p><p>A pump is usually specified as <code class=\"docutils literal notranslate\"><span class=\"pre\">Q</span></code>, <code class=\"docutils literal notranslate\"><span class=\"pre\">H</span></code> and <code class=\"docutils literal notranslate\"><span class=\"pre\">n</span></code>, but in an engine those\nfall out of a medium, a mass flow and two pressures. <a class=\"reference internal\" href=\"PumpDuty.html#pyskyfire.pump.duty.PumpDuty\" title=\"pyskyfire.pump.duty.PumpDuty\"><code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">PumpDuty</span></code></a> takes\neither and stores the same thing, so a design reads it the same way whichever\nthe user had to hand.</p>", "a[href=\"#classes\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Classes<a class=\"headerlink\" href=\"#classes\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"PumpDuty.html#pyskyfire.pump.duty.PumpDuty\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.pump.duty.PumpDuty\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.pump.duty.</span></span><span class=\"sig-name descname\"><span class=\"pre\">PumpDuty</span></span></dt><dd><p>Bases: <a class=\"reference internal\" href=\"../../thermo/duty/MachineDuty.html#pyskyfire.thermo.duty.MachineDuty\" title=\"pyskyfire.thermo.duty.MachineDuty\"><code class=\"xref py py-obj docutils literal notranslate\"><span class=\"pre\">pyskyfire.thermo.duty.MachineDuty</span></code></a></p><p>What a pump is asked to do. See <code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">MachineDuty</span></code>.</p><p>Build one with a constructor:</p><p>Pressures are total pressures, as the network carries them, so the head\nis total-to-total.</p></dd>"}
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
