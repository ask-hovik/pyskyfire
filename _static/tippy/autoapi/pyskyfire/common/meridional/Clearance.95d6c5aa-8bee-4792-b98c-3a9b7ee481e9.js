selector_to_html = {"a[href=\"#pyskyfire-common-meridional-clearance\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.common.meridional.Clearance<a class=\"headerlink\" href=\"#pyskyfire-common-meridional-clearance\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.common.meridional.Clearance.__post_init__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.meridional.Clearance.__post_init__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__post_init__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.common.meridional.Clearance\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.common.meridional.Clearance\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.common.meridional.</span></span><span class=\"sig-name descname\"><span class=\"pre\">Clearance</span></span></dt><dd><p>A blade edge that stops short of a passage wall.</p></dd>"}
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
