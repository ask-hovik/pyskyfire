selector_to_html = {"a[href=\"#pyskyfire-ducts-pipe-displaypath\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.ducts.pipe.DisplayPath<a class=\"headerlink\" href=\"#pyskyfire-ducts-pipe-displaypath\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.ducts.pipe.DisplayPath.connection_count\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.pipe.DisplayPath.connection_count\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">connection_count</span></span><span class=\"property\"><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/functions.html#int\" title=\"(in Python v3.15)\"><span class=\"pre\">int</span></a></span></dt><dd><p>Number of source connections the path was refined from.</p></dd>", "a[href=\"#pyskyfire.ducts.pipe.DisplayPath.__post_init__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.pipe.DisplayPath.__post_init__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__post_init__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.ducts.pipe.DisplayPath.station_count\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.pipe.DisplayPath.station_count\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">station_count</span></span><span class=\"property\"><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/functions.html#int\" title=\"(in Python v3.15)\"><span class=\"pre\">int</span></a></span></dt><dd><p>Number of source <code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">PathStation</span></code> objects the path came from.</p></dd>", "a[href=\"#pyskyfire.ducts.pipe.DisplayPath\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.pipe.DisplayPath\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.ducts.pipe.</span></span><span class=\"sig-name descname\"><span class=\"pre\">DisplayPath</span></span></dt><dd><p>A refined centerline for drawing, never for pressure-loss marching.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">points</span></code> [m] and <code class=\"docutils literal notranslate\"><span class=\"pre\">radius</span></code> [m] describe the swept circular bore, and\n<code class=\"docutils literal notranslate\"><span class=\"pre\">core_radius</span></code> [m] the core inside it, zero where there is none. The\npath carries its own provenance so results can follow it: <code class=\"docutils literal notranslate\"><span class=\"pre\">station</span></code> is\nthe fractional <code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">PathStation</span></code> index of every point, and\n<code class=\"docutils literal notranslate\"><span class=\"pre\">connection</span></code> is the originating connection index of every segment\nbetween consecutive points.</p></dd>"}
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
