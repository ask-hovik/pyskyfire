selector_to_html = {"a[href=\"#pyskyfire.viz.core._Adder.shape\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.core._Adder.shape\">\n<span class=\"sig-name descname\"><span class=\"pre\">shape</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"o\"><span class=\"pre\">**</span></span><span class=\"n\"><span class=\"pre\">kw</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.viz.core._Adder\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.core._Adder\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.core.</span></span><span class=\"sig-name descname\"><span class=\"pre\">_Adder</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">fig</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><span class=\"pre\">plotly.graph_objects.Figure</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.viz.core._Adder.vline\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.core._Adder.vline\">\n<span class=\"sig-name descname\"><span class=\"pre\">vline</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"o\"><span class=\"pre\">**</span></span><span class=\"n\"><span class=\"pre\">kw</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire-viz-core-adder\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.core._Adder<a class=\"headerlink\" href=\"#pyskyfire-viz-core-adder\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.viz.core._Adder.hline\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.core._Adder.hline\">\n<span class=\"sig-name descname\"><span class=\"pre\">hline</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"o\"><span class=\"pre\">**</span></span><span class=\"n\"><span class=\"pre\">kw</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.viz.core._Adder.annotation\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.core._Adder.annotation\">\n<span class=\"sig-name descname\"><span class=\"pre\">annotation</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"o\"><span class=\"pre\">**</span></span><span class=\"n\"><span class=\"pre\">kw</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>"}
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
