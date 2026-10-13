selector_to_html = {"a[href=\"#pyskyfire-ducts-pipe-pipecorrelationwarning\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.ducts.pipe.PipeCorrelationWarning<a class=\"headerlink\" href=\"#pyskyfire-ducts-pipe-pipecorrelationwarning\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.ducts.pipe.PipeCorrelationWarning\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.pipe.PipeCorrelationWarning\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">exception</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.ducts.pipe.</span></span><span class=\"sig-name descname\"><span class=\"pre\">PipeCorrelationWarning</span></span></dt><dd><p>Bases: <a class=\"reference external\" href=\"https://docs.python.org/3/builtins/exceptions.html#RuntimeWarning\" title=\"(in Python v3.15)\"><code class=\"xref py py-obj docutils literal notranslate\"><span class=\"pre\">RuntimeWarning</span></code></a></p><p>A duct loss correlation is being evaluated outside what it describes.</p></dd>"}
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
