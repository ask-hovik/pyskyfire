selector_to_html = {"a[href=\"#pyskyfire.pump.diffuser.DiffuserResult.loss_breakdown\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.pump.diffuser.DiffuserResult.loss_breakdown\">\n<span class=\"sig-name descname\"><span class=\"pre\">loss_breakdown</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire-pump-diffuser-diffuserresult\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.pump.diffuser.DiffuserResult<a class=\"headerlink\" href=\"#pyskyfire-pump-diffuser-diffuserresult\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.pump.diffuser.DiffuserResult\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.pump.diffuser.DiffuserResult\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.pump.diffuser.</span></span><span class=\"sig-name descname\"><span class=\"pre\">DiffuserResult</span></span></dt><dd><p>Losses of one diffuser operating point, G\u00fclich Table 3.8 (2).</p><p>Coefficients are referred to u2\u00b2/2, as G\u00fclich writes them. <code class=\"docutils literal notranslate\"><span class=\"pre\">vaneless</span></code>\nis the wall friction between impeller and vanes (T3.8.23), <code class=\"docutils literal notranslate\"><span class=\"pre\">inlet</span></code> the\nfriction in the triangular inlet region (T3.8.18), <code class=\"docutils literal notranslate\"><span class=\"pre\">shock</span></code> the loss of\ndecelerating c2 onto the throat (the 0.3 term of T3.8.19) and <code class=\"docutils literal notranslate\"><span class=\"pre\">channel</span></code>\nthe diffusing channel, 1 - cp - 1/AR\u00b2 times the throat head. <code class=\"docutils literal notranslate\"><span class=\"pre\">mixing</span></code>\nis the channel jets filling the outlet ring: they leave through part of\nthe circumference only, the vanes\u2019 trailing edges taking the rest, and\nmixing out keeps their swirl but loses their excess radial velocity, a\nBorda-Carnot loss on that component. This is the part of G\u00fclich\u2019s overflow\nloss that belongs to the diffuser; the energy left in the mixed flow is not\ncharged, because the volute collects it and charges its own loss.</p></dd>", "a[href=\"#pyskyfire.pump.diffuser.DiffuserResult.total_pressure_loss\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.pump.diffuser.DiffuserResult.total_pressure_loss\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">property</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">total_pressure_loss</span></span></dt><dd></dd>"}
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
