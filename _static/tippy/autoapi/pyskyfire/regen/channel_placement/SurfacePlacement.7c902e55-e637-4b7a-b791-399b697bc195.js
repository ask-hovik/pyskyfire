selector_to_html = {"a[href=\"#pyskyfire.regen.channel_placement.SurfacePlacement.dtheta_dx\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.channel_placement.SurfacePlacement.dtheta_dx\">\n<span class=\"sig-name descname\"><span class=\"pre\">dtheta_dx</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">x</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">contour</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">circuit</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">None</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Return d\u03b8/dx [rad/m] for the representative lane centerline.\nDefault: 0 (axial/vertical channels).</p></dd>", "a[href=\"#pyskyfire-regen-channel-placement-surfaceplacement\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.regen.channel_placement.SurfacePlacement<a class=\"headerlink\" href=\"#pyskyfire-regen-channel-placement-surfaceplacement\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.regen.channel_placement.SurfacePlacement\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.channel_placement.SurfacePlacement\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.regen.channel_placement.</span></span><span class=\"sig-name descname\"><span class=\"pre\">SurfacePlacement</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">n_channel_positions</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/functions.html#int\" title=\"(in Python v3.15)\"><span class=\"pre\">int</span></a></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">helix_angle</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">0.0</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">channel_angular_width</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">None</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Bases: <code class=\"xref py py-obj docutils literal notranslate\"><span class=\"pre\">ChannelPlacement</span></code></p><p>Placement model for surface-mounted cooling channels.</p><p>Channels are positioned just outside the main wall stack,\noffset by total thickness and corrected for local contour angle.</p></dd>"}
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
