selector_to_html = {"a[href=\"#pyskyfire.pump.disk_friction.DiskFriction.power\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.pump.disk_friction.DiskFriction.power\">\n<span class=\"sig-name descname\"><span class=\"pre\">power</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"keyword-only-separator o\"><abbr title=\"Keyword-only parameters separator (PEP 3102)\"><span class=\"pre\">*</span></abbr></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">rho</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">mu</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">n</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Friction power of every face [W] at one fluid state and speed.</p></dd>", "a[href=\"#pyskyfire.pump.disk_friction.DiskFriction.__post_init__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.pump.disk_friction.DiskFriction.__post_init__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__post_init__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire-pump-disk-friction-diskfriction\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.pump.disk_friction.DiskFriction<a class=\"headerlink\" href=\"#pyskyfire-pump-disk-friction-diskfriction\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.pump.disk_friction.DiskFriction\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.pump.disk_friction.DiskFriction\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.pump.disk_friction.</span></span><span class=\"sig-name descname\"><span class=\"pre\">DiskFriction</span></span></dt><dd><p>The disk friction of one impeller, all its wetted faces together.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">inner_radii</span></code> holds one entry per face: a shrouded impeller has two, the\nfront shroud starting at its eye seal and the back at the hub. An open\nimpeller has only the back one.</p></dd>"}
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
