selector_to_html = {"a[href=\"#pyskyfire.ducts.port.PortSite.require\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.port.PortSite.require\">\n<span class=\"sig-name descname\"><span class=\"pre\">require</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">name</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Return a field a port needs, or say which port wanted it.</p></dd>", "a[href=\"#pyskyfire-ducts-port-portsite\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.ducts.port.PortSite<a class=\"headerlink\" href=\"#pyskyfire-ducts-port-portsite\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.ducts.port.PortSite\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.port.PortSite\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.ducts.port.</span></span><span class=\"sig-name descname\"><span class=\"pre\">PortSite</span></span></dt><dd><p>What a passage knows about where its port attaches.</p></dd>", "a[href=\"../../common/rotation/index.html#pyskyfire.common.rotation.rotation_sign\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.rotation.rotation_sign\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.common.rotation.</span></span><span class=\"sig-name descname\"><span class=\"pre\">rotation_sign</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">value</span></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/functions.html#float\" title=\"(in Python v3.15)\"><span class=\"pre\">float</span></a></span></span></dt><dd><p><code class=\"docutils literal notranslate\"><span class=\"pre\">+1.0</span></code> for a counterclockwise sense, <code class=\"docutils literal notranslate\"><span class=\"pre\">-1.0</span></code> for a clockwise one.</p></dd>"}
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
