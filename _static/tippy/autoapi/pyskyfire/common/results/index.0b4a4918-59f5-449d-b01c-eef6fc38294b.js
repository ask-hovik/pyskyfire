selector_to_html = {"a[href=\"#classes\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Classes<a class=\"headerlink\" href=\"#classes\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"#module-pyskyfire.common.results\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.common.results<a class=\"headerlink\" href=\"#module-pyskyfire.common.results\" title=\"Link to this heading\">\u00b6</a></h1><h2>Classes<a class=\"headerlink\" href=\"#classes\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"Results.html#pyskyfire.common.results.Results\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.results.Results\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.common.results.</span></span><span class=\"sig-name descname\"><span class=\"pre\">Results</span></span></dt><dd><p>Bases: <a class=\"reference external\" href=\"https://docs.python.org/3/library/collections.abc.html#collections.abc.MutableMapping\" title=\"(in Python v3.15)\"><code class=\"xref py py-obj docutils literal notranslate\"><span class=\"pre\">collections.abc.MutableMapping</span></code></a></p><p>Dictionary-like container with dot access and persistence helpers.</p><p>Behaves as a normal mutable mapping while allowing attribute-style\naccess (<code class=\"docutils literal notranslate\"><span class=\"pre\">res.foo</span></code> \u2261 <code class=\"docutils literal notranslate\"><span class=\"pre\">res[\"foo\"]</span></code>) and built-in save/load using\n<code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">cloudpickle</span></code>.</p><p class=\"rubric\">Notes</p><p>Dot assignment automatically writes into <code class=\"xref py py-attr docutils literal notranslate\"><span class=\"pre\">_store</span></code>. Internal\nattributes such as <code class=\"docutils literal notranslate\"><span class=\"pre\">_store</span></code> and <code class=\"docutils literal notranslate\"><span class=\"pre\">metadata</span></code> are protected from\naccidental overwrite.</p></dd>"}
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
