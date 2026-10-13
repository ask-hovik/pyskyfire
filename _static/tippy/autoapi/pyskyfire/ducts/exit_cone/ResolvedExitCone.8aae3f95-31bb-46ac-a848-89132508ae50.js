selector_to_html = {"a[href=\"#pyskyfire-ducts-exit-cone-resolvedexitcone\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.ducts.exit_cone.ResolvedExitCone<a class=\"headerlink\" href=\"#pyskyfire-ducts-exit-cone-resolvedexitcone\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.ducts.exit_cone.ResolvedExitCone\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.exit_cone.ResolvedExitCone\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.ducts.exit_cone.</span></span><span class=\"sig-name descname\"><span class=\"pre\">ResolvedExitCone</span></span></dt><dd><p>Bases: <a class=\"reference internal\" href=\"../port/ResolvedPort.html#pyskyfire.ducts.port.ResolvedPort\" title=\"pyskyfire.ducts.port.ResolvedPort\"><code class=\"xref py py-obj docutils literal notranslate\"><span class=\"pre\">pyskyfire.ducts.port.ResolvedPort</span></code></a></p><p>An exit cone placed in machine coordinates.</p><p>Everything a duct has in common with any other is on\n<a class=\"reference internal\" href=\"../port/ResolvedPort.html#pyskyfire.ducts.port.ResolvedPort\" title=\"pyskyfire.ducts.port.ResolvedPort\"><code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">ResolvedPort</span></code></a>; what a cone adds is the pair\nof angles it was steered by, and where round the volute they carried it.\n<code class=\"docutils literal notranslate\"><span class=\"pre\">start_direction</span></code> is tangential at the volute outlet.</p></dd>", "a[href=\"../port/ResolvedPort.html#pyskyfire.ducts.port.ResolvedPort\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.port.ResolvedPort\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.ducts.port.</span></span><span class=\"sig-name descname\"><span class=\"pre\">ResolvedPort</span></span></dt><dd><p>A duct placed in machine coordinates.</p></dd>"}
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
