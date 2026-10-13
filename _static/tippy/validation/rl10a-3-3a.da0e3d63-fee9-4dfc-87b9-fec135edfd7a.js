selector_to_html = {"a[href=\"#rl10a-3-3a\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">RL10A-3-3A<a class=\"headerlink\" href=\"#rl10a-3-3a\" title=\"Link to this heading\">\u00b6</a></h1><p>The HTML documentation build replaces this page with the standalone\nRL10A-3-3A validation report. This source page exists so that Sphinx can include\nthe report in its navigation and so that non-HTML documentation builders have a\nuseful fallback.</p><p>The generated report is maintained at\n<code class=\"docutils literal notranslate\"><span class=\"pre\">validation/RL10/RL10A-3-3A_Report.html</span></code> in the repository.</p>"}
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
