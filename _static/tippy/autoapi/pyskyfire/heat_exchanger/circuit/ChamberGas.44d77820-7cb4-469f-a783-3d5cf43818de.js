selector_to_html = {"a[href=\"#pyskyfire.heat_exchanger.circuit.ChamberGas\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.heat_exchanger.circuit.ChamberGas\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.heat_exchanger.circuit.</span></span><span class=\"sig-name descname\"><span class=\"pre\">ChamberGas</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">chamber</span></span></em>, <em class=\"sig-param\"><span class=\"keyword-only-separator o\"><abbr title=\"Keyword-only parameters separator (PEP 3102)\"><span class=\"pre\">*</span></abbr></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">position</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">(0.0,</span> <span class=\"pre\">0.0,</span> <span class=\"pre\">0.0)</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">rotation</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">None</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">coefficient</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">None</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Adapter to an unchanged thrust-chamber aerothermodynamics solution.</p><p>Assembly coordinates map into chamber coordinates by rotation and position.\nA supplied coefficient [W/m\u00b2/K] or callback(position, wall_T) is preferred\nfor a specific baffle geometry. The default uses the existing Bartz-style\nreference-temperature correlation as an engineering approximation; it is\nnot a resolved flow-around-baffles model.</p></dd>", "a[href=\"#pyskyfire.heat_exchanger.circuit.ChamberGas.world_position\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.heat_exchanger.circuit.ChamberGas.world_position\">\n<span class=\"sig-name descname\"><span class=\"pre\">world_position</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">point</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire-heat-exchanger-circuit-chambergas\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.heat_exchanger.circuit.ChamberGas<a class=\"headerlink\" href=\"#pyskyfire-heat-exchanger-circuit-chambergas\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.heat_exchanger.circuit.ChamberGas.sample\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.heat_exchanger.circuit.ChamberGas.sample\">\n<span class=\"sig-name descname\"><span class=\"pre\">sample</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">point</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">wall_temperature</span></span></em><span class=\"sig-paren\">)</span></dt><dd></dd>"}
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
