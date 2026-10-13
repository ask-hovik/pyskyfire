selector_to_html = {"a[href=\"#pyskyfire.ducts.path.PathNode.__post_init__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.path.PathNode.__post_init__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__post_init__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire-ducts-path-pathnode\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.ducts.path.PathNode<a class=\"headerlink\" href=\"#pyskyfire-ducts-path-pathnode\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.ducts.path.PathNode\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.path.PathNode\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.ducts.path.</span></span><span class=\"sig-name descname\"><span class=\"pre\">PathNode</span></span></dt><dd><p>Human-friendly route vertex consumed by <code class=\"xref py py-func docutils literal notranslate\"><span class=\"pre\">generate_path_stations()</span></code>.</p></dd>"}
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
