selector_to_html = {"a[href=\"#pyskyfire.ducts.port.ResolvedRadialPort._start_frame\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.port.ResolvedRadialPort._start_frame\">\n<span class=\"sig-name descname\"><span class=\"pre\">_start_frame</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd><p>Across the section plane first, then out of it.</p><p>A stub\u2019s cross-section splits naturally into the direction the dome\ncurves away in, which lies in the section plane, and the azimuthal one\nthe swept passage carries straight on through. Trimming the tube to the\npassage only moves it along the first.</p></dd>", "a[href=\"#pyskyfire-ducts-port-resolvedradialport\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.ducts.port.ResolvedRadialPort<a class=\"headerlink\" href=\"#pyskyfire-ducts-port-resolvedradialport\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.ducts.port.ResolvedRadialPort.sweep\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.port.ResolvedRadialPort.sweep\">\n<span class=\"sig-name descname\"><span class=\"pre\">sweep</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">n_points</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">64</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">n_stations</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">9</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>The stub\u2019s wall, trimmed back onto the dome it leaves [m, (N, M, 3)].</p><p>The first ring follows the dome rather than lying flat against it, which\nis what stops the tube floating a chord above the passage. It is sunk a\nlittle further than this station\u2019s own ellipse asks for, because the\nsection shrinks either side of the port and a tube buried slightly in\nmetal reads better than one with daylight under it.</p></dd>", "a[href=\"#pyskyfire.ducts.port.ResolvedRadialPort._drawing_margin\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.port.ResolvedRadialPort._drawing_margin\">\n<span class=\"sig-name descname\"><span class=\"pre\">_drawing_margin</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd><p>How much further than its own ellipse to sink the first ring [m].</p><p>The trim is taken against the port station\u2019s ellipse alone, which is\ncheap but leaves a gap: the passage is swept round a circle and its\nsection shrinks away from the port, so both effects pull the real\nsurface inboard of that ellipse across the width the flange covers.</p></dd>", "a[href=\"#pyskyfire.ducts.port.ResolvedRadialPort\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.ducts.port.ResolvedRadialPort\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.ducts.port.</span></span><span class=\"sig-name descname\"><span class=\"pre\">ResolvedRadialPort</span></span></dt><dd><p>Bases: <code class=\"xref py py-obj docutils literal notranslate\"><span class=\"pre\">ResolvedPort</span></code></p><p>A stub placed in machine coordinates, with the dome it leaves through.</p></dd>"}
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
