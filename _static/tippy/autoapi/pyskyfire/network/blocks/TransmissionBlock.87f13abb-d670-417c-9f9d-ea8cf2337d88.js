selector_to_html = {"a[href=\"../engine_network/Station.html#pyskyfire.network.engine_network.Station\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.network.engine_network.Station\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.network.engine_network.</span></span><span class=\"sig-name descname\"><span class=\"pre\">Station</span></span></dt><dd><p>Thermodynamic state container used by the engine network.</p></dd>", "a[href=\"#pyskyfire-network-blocks-transmissionblock\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.network.blocks.TransmissionBlock<a class=\"headerlink\" href=\"#pyskyfire-network-blocks-transmissionblock\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.network.blocks.TransmissionBlock.compute\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.network.blocks.TransmissionBlock.compute\">\n<span class=\"sig-name descname\"><span class=\"pre\">compute</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">st</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">sg</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Sum input power signals.</p></dd>", "a[href=\"#pyskyfire.network.blocks.TransmissionBlock\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.network.blocks.TransmissionBlock\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.network.blocks.</span></span><span class=\"sig-name descname\"><span class=\"pre\">TransmissionBlock</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">name</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">sink_keys</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">source_keys</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Bases: <code class=\"xref py py-obj docutils literal notranslate\"><span class=\"pre\">SignalBlock</span></code></p><p>Sum multiple shaft-power signals into one required-power signal.</p><p>Useful for aggregating pump/compressor power demands prior to a turbine\nblock that must supply the total shaft power.</p></dd>"}
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
