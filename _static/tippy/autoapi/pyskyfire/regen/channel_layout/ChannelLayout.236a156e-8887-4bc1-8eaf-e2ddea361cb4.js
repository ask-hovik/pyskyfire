selector_to_html = {"a[href=\"#pyskyfire-regen-channel-layout-channellayout\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.regen.channel_layout.ChannelLayout<a class=\"headerlink\" href=\"#pyskyfire-regen-channel-layout-channellayout\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.regen.channel_layout.ChannelLayout\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.channel_layout.ChannelLayout\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.regen.channel_layout.</span></span><span class=\"sig-name descname\"><span class=\"pre\">ChannelLayout</span></span></dt><dd><p>A consumer-neutral collection of unique, circularly patternable seeds.</p></dd>", "a[href=\"ResolvedChannel.html#pyskyfire.regen.channel_layout.ResolvedChannel\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.channel_layout.ResolvedChannel\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.regen.channel_layout.</span></span><span class=\"sig-name descname\"><span class=\"pre\">ResolvedChannel</span></span></dt><dd><p>One unique channel seed and the circular pattern that reproduces it.</p></dd>", "a[href=\"#pyskyfire.regen.channel_layout.ChannelLayout.for_circuits\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.channel_layout.ChannelLayout.for_circuits\">\n<span class=\"sig-name descname\"><span class=\"pre\">for_circuits</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">circuit_indices</span></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference internal\" href=\"#pyskyfire.regen.channel_layout.ChannelLayout\" title=\"pyskyfire.regen.channel_layout.ChannelLayout\"><span class=\"pre\">ChannelLayout</span></a></span></span></dt><dd><p>Return a layout containing only the requested source circuits.</p></dd>", "a[href=\"#pyskyfire.regen.channel_layout.ChannelLayout.expanded_channels\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.channel_layout.ChannelLayout.expanded_channels\">\n<span class=\"sig-name descname\"><span class=\"pre\">expanded_channels</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://docs.python.org/3/builtins/stdtypes.html#tuple\" title=\"(in Python v3.15)\"><span class=\"pre\">tuple</span></a><span class=\"p\"><span class=\"pre\">[</span></span><a class=\"reference internal\" href=\"ResolvedChannel.html#pyskyfire.regen.channel_layout.ResolvedChannel\" title=\"pyskyfire.regen.channel_layout.ResolvedChannel\"><span class=\"pre\">ResolvedChannel</span></a><span class=\"p\"><span class=\"pre\">,</span></span><span class=\"w\"> </span><span class=\"p\"><span class=\"pre\">...</span></span><span class=\"p\"><span class=\"pre\">]</span></span></span></span></dt><dd><p>Expand every seed into its complete physical circular pattern.</p></dd>"}
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
