selector_to_html = {"a[href=\"#pyskyfire.ducts.port.Port.resolve\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.port.Port.resolve\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">abstractmethod</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">resolve</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">site</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference internal\" href=\"PortSite.html#pyskyfire.ducts.port.PortSite\" title=\"pyskyfire.ducts.port.PortSite\"><span class=\"pre\">PortSite</span></a></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference internal\" href=\"ResolvedPort.html#pyskyfire.ducts.port.ResolvedPort\" title=\"pyskyfire.ducts.port.ResolvedPort\"><span class=\"pre\">ResolvedPort</span></a></span></span></dt><dd><p>Place this port on a passage, given what the passage knows.</p></dd>", "a[href=\"#pyskyfire-ducts-port-port\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.ducts.port.Port<a class=\"headerlink\" href=\"#pyskyfire-ducts-port-port\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.ducts.port.Port\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.port.Port\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.ducts.port.</span></span><span class=\"sig-name descname\"><span class=\"pre\">Port</span></span></dt><dd><p>Bases: <a class=\"reference external\" href=\"https://docs.python.org/3/library/abc.html#abc.ABC\" title=\"(in Python v3.15)\"><code class=\"xref py py-obj docutils literal notranslate\"><span class=\"pre\">abc.ABC</span></code></a></p><p>A request for the duct carrying collected flow to its flange.</p></dd>", "a[href=\"PortSite.html#pyskyfire.ducts.port.PortSite\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.port.PortSite\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.ducts.port.</span></span><span class=\"sig-name descname\"><span class=\"pre\">PortSite</span></span></dt><dd><p>What a passage knows about where its port attaches.</p></dd>", "a[href=\"ResolvedPort.html#pyskyfire.ducts.port.ResolvedPort\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.port.ResolvedPort\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.ducts.port.</span></span><span class=\"sig-name descname\"><span class=\"pre\">ResolvedPort</span></span></dt><dd><p>A duct placed in machine coordinates.</p></dd>"}
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
