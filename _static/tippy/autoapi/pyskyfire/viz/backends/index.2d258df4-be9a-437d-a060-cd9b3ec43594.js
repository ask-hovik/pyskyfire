selector_to_html = {"a[href=\"#submodules\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Submodules<a class=\"headerlink\" href=\"#submodules\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"common/index.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.backends.common<a class=\"headerlink\" href=\"#module-pyskyfire.viz.backends.common\" title=\"Link to this heading\">\u00b6</a></h1><p>Scalar display validation shared by desktop and HTML backends.</p>", "a[href=\"threejs/index.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.backends.threejs<a class=\"headerlink\" href=\"#module-pyskyfire.viz.backends.threejs\" title=\"Link to this heading\">\u00b6</a></h1><p>Offline HTML rendering with mesh references and original scalar values.</p><p>The scene is serialized directly instead of exporting a VTK rendering window\nthrough glTF, so fields and node identities survive the browser boundary.</p>", "a[href=\"pyvista/index.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.backends.pyvista<a class=\"headerlink\" href=\"#module-pyskyfire.viz.backends.pyvista\" title=\"Link to this heading\">\u00b6</a></h1><p>Desktop rendering of a scene, retaining named scalar arrays and node IDs.</p>", "a[href=\"#module-pyskyfire.viz.backends\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.viz.backends<a class=\"headerlink\" href=\"#module-pyskyfire.viz.backends\" title=\"Link to this heading\">\u00b6</a></h1><p>Rendering adapters for the common scene model.</p>"}
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
