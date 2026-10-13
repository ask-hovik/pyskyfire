selector_to_html = {"a[href=\"#pyskyfire-viz-heat-exchanger-component-heatexchangercomponent\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.heat_exchanger_component.HeatExchangerComponent<a class=\"headerlink\" href=\"#pyskyfire-viz-heat-exchanger-component-heatexchangercomponent\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.viz.heat_exchanger_component.HeatExchangerComponent.build\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.heat_exchanger_component.HeatExchangerComponent.build\">\n<span class=\"sig-name descname\"><span class=\"pre\">build</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd></dd>", "a[href=\"#pyskyfire.viz.heat_exchanger_component.HeatExchangerComponent\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.viz.heat_exchanger_component.HeatExchangerComponent\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.viz.heat_exchanger_component.</span></span><span class=\"sig-name descname\"><span class=\"pre\">HeatExchangerComponent</span></span></dt><dd><p>Display an exchanger or chamber installation in the shared Scene.</p><p>Each channel is drawn as its wetted passage, the surface the coolant\ntouches, not the outer envelope of its walls. Raw exchangers use assembly\ncoordinates. Installations carry their simulation transform, relative to\nthe parent chamber; explicit Scene.add placement\noverrides it. Each channel has inlet/outlet anchors. Root anchors\ncircuit_N.inlet/outlet refer to the first copy; origin is the assembly datum.</p><p>Optional results provide T_static, p_static, velocity, coolant_enthalpy\nand T_coolant_wall, the wall temperature at the wetted interface. Wall\ntemperatures are area averages per exposed surface and segment. Inactive\nportions remain NaN. The opposite (gas-side) interface is not drawn.\nWith results, each circuit records its CoolProp fluid name as\n<code class=\"docutils literal notranslate\"><span class=\"pre\">metadata[\"fluid\"]</span></code>. <code class=\"docutils literal notranslate\"><span class=\"pre\">circuit_colors</span></code> paints each circuit apart;\nwithout it every channel takes the default colour.</p></dd>"}
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
