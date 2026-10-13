selector_to_html = {"a[href=\"#pyskyfire-pump-diffuser-diffuseroutlet\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.pump.diffuser.DiffuserOutlet<a class=\"headerlink\" href=\"#pyskyfire-pump-diffuser-diffuseroutlet\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.pump.diffuser.DiffuserOutlet\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.pump.diffuser.DiffuserOutlet\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.pump.diffuser.</span></span><span class=\"sig-name descname\"><span class=\"pre\">DiffuserOutlet</span></span></dt><dd><p>What the diffuser hands the volute, circumferentially mixed.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">cm</span></code> is the mean radial velocity over the whole outlet ring and <code class=\"docutils literal notranslate\"><span class=\"pre\">cu</span></code>\nthe mass-averaged swirl of the channel jets, both in m/s. Together with\n<code class=\"docutils literal notranslate\"><span class=\"pre\">radius</span></code> and <code class=\"docutils literal notranslate\"><span class=\"pre\">width</span></code> they are a volute inlet.</p></dd>"}
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
