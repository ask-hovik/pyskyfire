selector_to_html = {"a[href=\"#pyskyfire.ducts.orifice.OrificePlateResult\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.orifice.OrificePlateResult\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.ducts.orifice.</span></span><span class=\"sig-name descname\"><span class=\"pre\">OrificePlateResult</span></span></dt><dd><p>Bases: <code class=\"xref py py-obj docutils literal notranslate\"><span class=\"pre\">OrificeResult</span></code></p><p>An <code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">OrificeResult</span></code> with the plate\u2019s metering quantities.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">throat_state</span></code> is at the downstream tap, and\n<a class=\"reference internal\" href=\"#pyskyfire.ducts.orifice.OrificePlateResult.differential_pressure\" title=\"pyskyfire.ducts.orifice.OrificePlateResult.differential_pressure\"><code class=\"xref py py-attr docutils literal notranslate\"><span class=\"pre\">differential_pressure</span></code></a> is the tap-to-tap drop a flow meter reads.\n<code class=\"docutils literal notranslate\"><span class=\"pre\">outlet_state</span></code> is where the flow has recovered to the pipe bore, so\n<code class=\"xref py py-attr docutils literal notranslate\"><span class=\"pre\">pressure_loss</span></code> is the permanent loss, which is smaller.</p></dd>", "a[href=\"#pyskyfire-ducts-orifice-orificeplateresult\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.ducts.orifice.OrificePlateResult<a class=\"headerlink\" href=\"#pyskyfire-ducts-orifice-orificeplateresult\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.ducts.orifice.OrificePlateResult.differential_pressure\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.orifice.OrificePlateResult.differential_pressure\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">differential_pressure</span></span><span class=\"property\"><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/functions.html#float\" title=\"(in Python v3.15)\"><span class=\"pre\">float</span></a></span></dt><dd><p>Drop between the upstream and downstream taps [Pa].</p></dd>"}
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
