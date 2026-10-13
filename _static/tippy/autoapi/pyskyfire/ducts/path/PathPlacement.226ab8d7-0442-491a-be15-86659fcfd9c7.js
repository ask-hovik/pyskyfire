selector_to_html = {"a[href=\"#pyskyfire-ducts-path-pathplacement\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.ducts.path.PathPlacement<a class=\"headerlink\" href=\"#pyskyfire-ducts-path-pathplacement\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.ducts.path.PathPlacement.frames\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.path.PathPlacement.frames\">\n<span class=\"sig-name descname\"><span class=\"pre\">frames</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">path</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.ducts.path.PathPlacement\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.path.PathPlacement\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.ducts.path.</span></span><span class=\"sig-name descname\"><span class=\"pre\">PathPlacement</span></span></dt><dd><p>Place the first path station and establish its section roll.</p><p>The rotation maps template coordinates into assembly coordinates. The\ninitial floor direction is expressed in assembly coordinates, points from\nthe passage centroid toward its floor, and is projected perpendicular to\nthe starting tangent. Frame columns are (floor-to-roof, lateral, tangent).\nAll subsequent section orientations use minimal-rotation transport.</p></dd>", "a[href=\"#pyskyfire.ducts.path.PathPlacement.place\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.path.PathPlacement.place\">\n<span class=\"sig-name descname\"><span class=\"pre\">place</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">path</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.ducts.path.PathPlacement.__post_init__\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.path.PathPlacement.__post_init__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__post_init__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>"}
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
