selector_to_html = {"a[href=\"#pyskyfire-injector-element-hole\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.injector.element.Hole<a class=\"headerlink\" href=\"#pyskyfire-injector-element-hole\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.injector.element.Hole\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.injector.element.Hole\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.injector.element.</span></span><span class=\"sig-name descname\"><span class=\"pre\">Hole</span></span></dt><dd><p>One hole of an element, placed in the element\u2019s own frame.</p><p>The frame has its origin at the element centre on the face, +x along the\nchamber axis into the chamber, and +y along the element\u2019s orientation.\n<code class=\"docutils literal notranslate\"><span class=\"pre\">exit</span></code> is where the hole meets the face [m], <code class=\"docutils literal notranslate\"><span class=\"pre\">direction</span></code> a unit\nvector along the discharge, and <code class=\"docutils literal notranslate\"><span class=\"pre\">target</span></code> where the jet is drawn to and\nchecked against the wall: the impingement point for an impinging element.</p></dd>"}
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
