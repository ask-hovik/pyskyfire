selector_to_html = {"a[href=\"#pyskyfire-regen-coupled-solver-boundaryconditions\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.regen.coupled_solver.BoundaryConditions<a class=\"headerlink\" href=\"#pyskyfire-regen-coupled-solver-boundaryconditions\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.regen.coupled_solver.BoundaryConditions\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.coupled_solver.BoundaryConditions\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.regen.coupled_solver.</span></span><span class=\"sig-name descname\"><span class=\"pre\">BoundaryConditions</span></span></dt><dd><p>Boundary conditions for the coolant-side inlet.</p></dd>"}
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
