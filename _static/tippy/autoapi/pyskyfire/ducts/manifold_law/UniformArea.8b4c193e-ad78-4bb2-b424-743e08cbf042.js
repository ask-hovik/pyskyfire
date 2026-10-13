selector_to_html = {"a[href=\"#pyskyfire.ducts.manifold_law.UniformArea.__post_init__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.manifold_law.UniformArea.__post_init__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__post_init__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.ducts.manifold_law.UniformArea\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.manifold_law.UniformArea\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.ducts.manifold_law.</span></span><span class=\"sig-name descname\"><span class=\"pre\">UniformArea</span></span></dt><dd><p>Bases: <code class=\"xref py py-obj docutils literal notranslate\"><span class=\"pre\">ManifoldLaw</span></code></p><p>Hold the passage area constant all the way round.</p><p>The plain toroidal manifold. It carries more dynamic head near the port\nthan at the dead point, since the same area passes a growing share of the\nflow, so it relies on being oversized enough that the variation does not\nreach the channels. That is the usual rule for a coolant manifold, and it\nis why the parameter here is simply the area the designer chose.</p></dd>", "a[href=\"#pyskyfire.ducts.manifold_law.UniformArea.residual\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.manifold_law.UniformArea.residual\">\n<span class=\"sig-name descname\"><span class=\"pre\">residual</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">collected</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">area</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Return the signed miss at one station.</p><p>Negative means the trial section is too small. Must increase with\nsection area so that the solver can bracket a root.</p></dd>", "a[href=\"#pyskyfire-ducts-manifold-law-uniformarea\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.ducts.manifold_law.UniformArea<a class=\"headerlink\" href=\"#pyskyfire-ducts-manifold-law-uniformarea\" title=\"Link to this heading\">\u00b6</a></h1>"}
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
