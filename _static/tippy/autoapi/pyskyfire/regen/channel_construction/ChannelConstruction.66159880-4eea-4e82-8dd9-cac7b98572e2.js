selector_to_html = {"a[href=\"#pyskyfire.regen.channel_construction.ChannelConstruction.centroid_path\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.channel_construction.ChannelConstruction.centroid_path\">\n<span class=\"sig-name descname\"><span class=\"pre\">centroid_path</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">cross_section</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Return the section centroid in world cylindrical coordinates.</p></dd>", "a[href=\"#pyskyfire.regen.channel_construction.ChannelConstruction\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.channel_construction.ChannelConstruction\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.regen.channel_construction.</span></span><span class=\"sig-name descname\"><span class=\"pre\">ChannelConstruction</span></span></dt><dd><p>One circuit\u2019s resolved section geometry along a nominal winding.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">planes</span></code> are what the section boundary is drawn in \u2013 transverse for a\nsquared section, flow-normal for a rounded one \u2013 while <code class=\"docutils literal notranslate\"><span class=\"pre\">flow_planes</span></code> are\nalways flow-normal and are what the centroid is measured in. <code class=\"docutils literal notranslate\"><span class=\"pre\">profiles</span></code>\ncarries the resolved dimensions, including the base-point path the section\nis anchored on, and <code class=\"docutils literal notranslate\"><span class=\"pre\">gas_area</span></code> the combustion-gas area each channel\u2019s\nsector keeps.</p></dd>", "a[href=\"#pyskyfire-regen-channel-construction-channelconstruction\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.regen.channel_construction.ChannelConstruction<a class=\"headerlink\" href=\"#pyskyfire-regen-channel-construction-channelconstruction\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.regen.channel_construction.ChannelConstruction.polygon\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.channel_construction.ChannelConstruction.polygon\">\n<span class=\"sig-name descname\"><span class=\"pre\">polygon</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">index</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">n_arc_points</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">8</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>"}
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
