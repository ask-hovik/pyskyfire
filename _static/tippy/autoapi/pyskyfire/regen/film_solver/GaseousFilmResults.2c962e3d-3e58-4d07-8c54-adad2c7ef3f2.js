selector_to_html = {"a[href=\"#pyskyfire-regen-film-solver-gaseousfilmresults\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.regen.film_solver.GaseousFilmResults<a class=\"headerlink\" href=\"#pyskyfire-regen-film-solver-gaseousfilmresults\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.regen.film_solver.GaseousFilmResults\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.regen.film_solver.GaseousFilmResults\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.regen.film_solver.</span></span><span class=\"sig-name descname\"><span class=\"pre\">GaseousFilmResults</span></span></dt><dd></dd>"}
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
