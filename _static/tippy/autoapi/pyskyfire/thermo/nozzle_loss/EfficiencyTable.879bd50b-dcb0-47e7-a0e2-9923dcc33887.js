selector_to_html = {"a[href=\"#pyskyfire-thermo-nozzle-loss-efficiencytable\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.thermo.nozzle_loss.EfficiencyTable<a class=\"headerlink\" href=\"#pyskyfire-thermo-nozzle-loss-efficiencytable\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.thermo.nozzle_loss.EfficiencyTable\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.thermo.nozzle_loss.EfficiencyTable\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.thermo.nozzle_loss.</span></span><span class=\"sig-name descname\"><span class=\"pre\">EfficiencyTable</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">mixture_ratio</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">chamber_pressure</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">values</span></span></em>, <em class=\"sig-param\"><span class=\"keyword-only-separator o\"><abbr title=\"Keyword-only parameters separator (PEP 3102)\"><span class=\"pre\">*</span></abbr></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">log_pressure</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/functions.html#bool\" title=\"(in Python v3.15)\"><span class=\"pre\">bool</span></a></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">True</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Efficiency tabulated against mixture ratio and chamber pressure.</p><p class=\"rubric\">Notes</p><p>Queries outside the tabulated ranges are <strong>clamped</strong> to the edge values\nrather than extrapolated. Efficiency tables of this kind are fits over a\nstated envelope and go badly wrong outside it \u2013 the reference report notes\nthat the vendor tables it compared against can return negative specific\nimpulse when pushed past their intended range.</p></dd>", "a[href=\"#pyskyfire.thermo.nozzle_loss.EfficiencyTable.__call__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.thermo.nozzle_loss.EfficiencyTable.__call__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__call__</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">p_c</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">MR</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Interpolate the table at <code class=\"docutils literal notranslate\"><span class=\"pre\">p_c</span></code> [Pa] and mixture ratio <code class=\"docutils literal notranslate\"><span class=\"pre\">MR</span></code> [-].</p><p>Scalars return a float; array inputs broadcast against each other and\nreturn an array of the broadcast shape.</p></dd>"}
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
