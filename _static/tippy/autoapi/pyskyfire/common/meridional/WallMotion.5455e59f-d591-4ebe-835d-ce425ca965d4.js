selector_to_html = {"a[href=\"#pyskyfire.common.meridional.WallMotion.rotating\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.meridional.WallMotion.rotating\">\n<span class=\"sig-name descname\"><span class=\"pre\">rotating</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">wall</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">fraction</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Whether <code class=\"docutils literal notranslate\"><span class=\"pre\">wall</span></code> turns at each station <code class=\"docutils literal notranslate\"><span class=\"pre\">fraction</span></code>. Ends count as inside.</p></dd>", "a[href=\"#pyskyfire.common.meridional.WallMotion.whole\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.meridional.WallMotion.whole\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">classmethod</span></span><span class=\"w\"> </span></span><span class=\"sig-name descname\"><span class=\"pre\">whole</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"o\"><span class=\"pre\">*</span></span><span class=\"n\"><span class=\"pre\">walls</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>The named walls rotate over their whole length, the rest stand still.</p></dd>", "a[href=\"#pyskyfire.common.meridional.WallMotion\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.meridional.WallMotion\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.common.meridional.</span></span><span class=\"sig-name descname\"><span class=\"pre\">WallMotion</span></span></dt><dd><p>Which stretches of each wall rotate with the shaft.</p><p>Each wall lists the <code class=\"docutils literal notranslate\"><span class=\"pre\">(start,</span> <span class=\"pre\">end)</span></code> station fractions over which it\nturns; everywhere else it stands still. The default is a passage that\nstands still throughout, as a duct does.</p></dd>", "a[href=\"#pyskyfire.common.meridional.WallMotion.__post_init__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.meridional.WallMotion.__post_init__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__post_init__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire-common-meridional-wallmotion\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.common.meridional.WallMotion<a class=\"headerlink\" href=\"#pyskyfire-common-meridional-wallmotion\" title=\"Link to this heading\">\u00b6</a></h1>"}
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
