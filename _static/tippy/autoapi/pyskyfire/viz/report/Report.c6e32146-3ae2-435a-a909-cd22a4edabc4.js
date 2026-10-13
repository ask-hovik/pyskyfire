selector_to_html = {"a[href=\"#pyskyfire.viz.report.Report\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.report.Report\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.report.</span></span><span class=\"sig-name descname\"><span class=\"pre\">Report</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">title</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/stdtypes.html#str\" title=\"(in Python v3.15)\"><span class=\"pre\">str</span></a></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">'pyskyfire</span> <span class=\"pre\">Report'</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.viz.report.Report.save_html\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.report.Report.save_html\">\n<span class=\"sig-name descname\"><span class=\"pre\">save_html</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">path</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/stdtypes.html#str\" title=\"(in Python v3.15)\"><span class=\"pre\">str</span></a></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"Tab.html#pyskyfire.viz.report.Tab\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.report.Tab\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.report.</span></span><span class=\"sig-name descname\"><span class=\"pre\">Tab</span></span></dt><dd></dd>", "a[href=\"#pyskyfire.viz.report.Report.add_tab\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.report.Report.add_tab\">\n<span class=\"sig-name descname\"><span class=\"pre\">add_tab</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">title</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/stdtypes.html#str\" title=\"(in Python v3.15)\"><span class=\"pre\">str</span></a></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference internal\" href=\"Tab.html#pyskyfire.viz.report.Tab\" title=\"pyskyfire.viz.report.Tab\"><span class=\"pre\">Tab</span></a></span></span></dt><dd></dd>", "a[href=\"#pyskyfire-viz-report-report\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.report.Report<a class=\"headerlink\" href=\"#pyskyfire-viz-report-report\" title=\"Link to this heading\">\u00b6</a></h1>"}
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
