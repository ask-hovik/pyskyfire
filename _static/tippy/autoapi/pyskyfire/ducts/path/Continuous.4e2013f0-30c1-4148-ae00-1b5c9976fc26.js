selector_to_html = {"a[href=\"#pyskyfire.ducts.path.Continuous.__post_init__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.path.Continuous.__post_init__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__post_init__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.ducts.path.Continuous\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.path.Continuous\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.ducts.path.</span></span><span class=\"sig-name descname\"><span class=\"pre\">Continuous</span></span></dt><dd><p>A locally resolved continuous-curvature connection.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">radius</span></code> is the exact centerline radius when known. If omitted on a\nmanually supplied layout, <code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">Path</span></code> estimates it transparently from\nneighboring station triples using Menger curvature. Attached to a\n<code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">PathNode</span></code>, it rounds that corner with a circular arc of <code class=\"docutils literal notranslate\"><span class=\"pre\">radius</span></code>,\nwhich the generator emits as one connection carrying its exact arc\n<code class=\"docutils literal notranslate\"><span class=\"pre\">length</span></code>, <code class=\"docutils literal notranslate\"><span class=\"pre\">center</span></code> and <code class=\"docutils literal notranslate\"><span class=\"pre\">normal</span></code>. A constant-section arc needs no\nfurther stations; <code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">Pipe</span></code> adds them where its section varies.</p></dd>", "a[href=\"#pyskyfire-ducts-path-continuous\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.ducts.path.Continuous<a class=\"headerlink\" href=\"#pyskyfire-ducts-path-continuous\" title=\"Link to this heading\">\u00b6</a></h1>"}
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
