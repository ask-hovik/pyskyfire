selector_to_html = {"a[href=\"#pyskyfire-pump-diffuser-diffuserchecks\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pyskyfire.pump.diffuser.DiffuserChecks<a class=\"headerlink\" href=\"#pyskyfire-pump-diffuser-diffuserchecks\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#pyskyfire.pump.diffuser.DiffuserChecks\"]": "<dt class=\"sig sig-object py\" id=\"pyskyfire.pump.diffuser.DiffuserChecks\">\n<span class=\"property\"><span class=\"k\"><span class=\"pre\">class</span></span><span class=\"w\"> </span></span><span class=\"sig-prename descclassname\"><span class=\"pre\">pyskyfire.pump.diffuser.</span></span><span class=\"sig-name descname\"><span class=\"pre\">DiffuserChecks</span></span></dt><dd><p>The finished diffuser against G\u00fclich\u2019s criteria (7.9.1 item 10). Reported only.</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">gap_ratio</span></code> is d3/d2 against <code class=\"docutils literal notranslate\"><span class=\"pre\">minimum_gap_ratio</span></code> (Table 10.2).\n<code class=\"docutils literal notranslate\"><span class=\"pre\">inlet_width_ratio</span></code> b3/b2 should be 1.05 to 1.3. <code class=\"docutils literal notranslate\"><span class=\"pre\">throat_factor</span></code> is\nf_a3, normally 1.1 to 1.3. <code class=\"docutils literal notranslate\"><span class=\"pre\">deceleration</span></code> is c3q/c2: 0.7 to 0.85,\nhigher allowed below nq 12 to 15, and Table 10.2 asks for not much below\n0.8. <code class=\"docutils literal notranslate\"><span class=\"pre\">length_ratio</span></code> L/a3 should lie between 2.5 and 6 and\n<code class=\"docutils literal notranslate\"><span class=\"pre\">aspect_ratio</span></code> b3/a3 between 0.8 and 2. <code class=\"docutils literal notranslate\"><span class=\"pre\">exit_energy</span></code> is c4\u00b2/(2gH),\n0.02 to 0.04 for multistage pumps. <code class=\"docutils literal notranslate\"><span class=\"pre\">opening_angle_deg</span></code> is the equivalent\nconical half-angle (Eq. 1.45), which must stay below\n<code class=\"docutils literal notranslate\"><span class=\"pre\">allowed_opening_angle_deg</span></code> (Eq. 1.44) or the channel stalls.\n<code class=\"docutils literal notranslate\"><span class=\"pre\">vane_count_conflicts</span></code> lists each (nu2, nu3, m) with\n<code class=\"docutils literal notranslate\"><span class=\"pre\">|m|</span> <span class=\"pre\">=</span> <span class=\"pre\">|nu3</span> <span class=\"pre\">z_Le</span> <span class=\"pre\">-</span> <span class=\"pre\">nu2</span> <span class=\"pre\">z_La|</span></code> of 0 or 1 for nu up to 3, which G\u00fclich says\nto avoid (7.9.1 item 5, Table 10.2). <code class=\"docutils literal notranslate\"><span class=\"pre\">incidence_deg</span></code> is the vane\u2019s\nleading spiral against the approach flow, within +-3 degrees by his rule;\nit follows from the throat. <code class=\"docutils literal notranslate\"><span class=\"pre\">deflection_deg</span></code> is how far the channel\nturns the flow from throat to exit; G\u00fclich applies straight-diffuser data\nup to about 25 degrees. <code class=\"docutils literal notranslate\"><span class=\"pre\">drawn_throat_ratio</span></code> and <code class=\"docutils literal notranslate\"><span class=\"pre\">drawn_exit_ratio</span></code>\nare the throat and exit widths measured on the drawn vanes over the design\nvalues, which the camber construction meets only approximately.</p></dd>"}
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
