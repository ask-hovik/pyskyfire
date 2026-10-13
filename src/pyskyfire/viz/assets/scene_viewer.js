/* Offline scene display. Clicking a part selects the top-level component that
 * owns it, so one cooling channel identifies its whole thrust chamber, and
 * reports the placement arguments that reproduce its position in Python.
 * Shift+click gathers several components; right-click lists every part under
 * the cursor, down to circuits buried inside other components. A selection can
 * be saved as a named group, and any group can set the colour range on its own
 * while the rest of the scene fades back. */
(function () {
  "use strict";
  const hud = document.getElementById("hud");
  try {
    const data = JSON.parse(document.getElementById("psf-scene").textContent);
    const renderer = new THREE.WebGLRenderer({canvas: document.getElementById("c"), antialias: true});
    renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    renderer.setSize(window.innerWidth, window.innerHeight, false);
    const scene = new THREE.Scene();
    scene.background = new THREE.Color("white");
    const camera = new THREE.PerspectiveCamera(45, window.innerWidth/window.innerHeight, 0.001, 1000);
    // Trackball rather than orbit: rotation is unconstrained, so the view can be
    // spun continuously through the poles like a CAD viewer.
    const controls = new THREE.TrackballControls(camera, renderer.domElement);
    controls.rotateSpeed = 2.5;
    controls.zoomSpeed = 1.2;
    controls.panSpeed = 0.6;
    controls.dynamicDampingFactor = 0.25;
    // Ctrl/Cmd + left-drag pans, for trackpads without a comfortable right-drag.
    // The controls pick an action by matching event.button against these slots,
    // so pointing the pan slot at the left button for this press is enough. The
    // listener captures so it runs before the controls' own pointerdown.
    const buttons = {...controls.mouseButtons};
    renderer.domElement.addEventListener("pointerdown", event => {
      const pan = event.button === 0 && (event.ctrlKey || event.metaKey);
      Object.assign(controls.mouseButtons, pan ? {LEFT: -1, RIGHT: 0} : buttons);
    }, {capture: true});
    scene.add(new THREE.AmbientLight(0xffffff, 0.65));
    const light = new THREE.DirectionalLight(0xffffff, 0.85);
    light.position.set(5, 8, 6);
    scene.add(light);
    const geometry = data.meshes.map(m => {
      const g = new THREE.BufferGeometry();
      g.setAttribute("position", new THREE.Float32BufferAttribute(m.points, 3));
      g.setIndex(m.triangles);
      g.computeVertexNormals();
      return g;
    });
    function color(value, limits) {
      if (value === null) return [0.5, 0.5, 0.5];
      const [lo, hi] = limits;
      const index = Math.max(0, Math.min(255, Math.floor(256*(value-lo)/(hi-lo))));
      return data.color_table[index].map(x => x/255);
    }
    // Meshes are kept flat alongside their node, so the selected color field can
    // be swapped without rebuilding the hierarchy. Every node is also found by
    // its path ID, which is how selections and groups refer to parts.
    const parts = [];
    const nodes = {};
    function add(node, parent) {
      const group = new THREE.Group();
      group.name = node.id;
      group.userData = {node: node, fields: node.fields, metadata: node.metadata, anchors: node.anchors};
      group.matrix.fromArray(node.transform);
      group.matrixAutoUpdate = false;
      group.visible = node.visible;
      parent.add(group);
      nodes[node.id] = group;
      // Field values arrive by reference into a shared table, as meshes do.
      Object.values(node.fields).forEach(field => { field.values = data.values[field.values]; });
      if (node.mesh !== null) {
        const material = new THREE.MeshPhongMaterial({
          side: THREE.DoubleSide, transparent: node.opacity < 1, opacity: node.opacity
        });
        const mesh = new THREE.Mesh(geometry[node.mesh], material);
        mesh.name = node.id;
        group.add(mesh);
        // Seed geometry is shared between repeated channels, so colored copies
        // are created per part and cached by association rather than mutated.
        const part = {node: node, mesh: mesh, base: geometry[node.mesh], variants: {}, surface: node.color};
        mesh.userData.part = part;
        parts.push(part);
      }
      node.children.forEach(child => add(child, group));
    }
    const root = new THREE.Group();
    scene.add(root);
    data.nodes.forEach(node => add(node, root));
    root.updateMatrixWorld(true);

    function meshesOf(ids) {
      const meshes = new Set();
      ids.forEach(id => nodes[id].traverse(o => { if (o.isMesh) meshes.add(o); }));
      return meshes;
    }
    function shown(object) {
      for (; object; object = object.parent) if (!object.visible) return false;
      return true;
    }
    function topOf(object) {
      while (object.parent && object.parent !== root) object = object.parent;
      return object;
    }
    function title(node) {
      return node.metadata && node.metadata.label ? node.name + " (" + node.metadata.label + ")" : node.name;
    }
    function pathName(id) {
      const names = id.split("/");
      names[names.length - 1] = title(nodes[id].userData.node);
      return names.join(" › ");
    }
    function fluidOf(group) {
      for (; group && group !== root; group = group.parent) {
        const fluid = group.userData.metadata && group.userData.metadata.fluid;
        if (fluid) return fluid;
      }
      return null;
    }
    function element(tag, className, text) {
      const e = document.createElement(tag);
      if (className) e.className = className;
      if (text !== undefined) e.textContent = text;
      return e;
    }

    // Groups: one per fluid a component recorded from its results, one per tag
    // written in Python, and those saved from a selection here. A node already
    // covered by an ancestor with the same fluid or tag is not listed again.
    function collect(test) {
      const ids = [];
      (function visit(list, covered) {
        list.forEach(node => {
          const hit = !covered && test(node);
          if (hit) ids.push(node.id);
          visit(node.children, covered || hit);
        });
      })(data.nodes, false);
      return ids;
    }
    const fluids = new Set(), tags = new Set();
    (function visit(list) {
      list.forEach(node => {
        if (node.metadata && node.metadata.fluid) fluids.add(node.metadata.fluid);
        node.tags.forEach(tag => tags.add(tag));
        visit(node.children);
      });
    })(data.nodes);
    const builtin = [
      ...[...fluids].map(fluid => ({kind: "Fluid", name: fluid,
                                   members: collect(node => node.metadata && node.metadata.fluid === fluid)})),
      ...[...tags].map(tag => ({kind: "Tag", name: tag, members: collect(node => node.tags.includes(tag))})),
    ];
    // Saved groups outlive the page where the browser allows local storage
    // (not in the sandboxed iframes a report embeds), keyed by the scene title.
    // Members missing from this scene are kept in storage but not shown.
    const storageKey = "pyskyfire-groups:" + data.title;
    let saved = [];
    try { saved = JSON.parse(window.localStorage.getItem(storageKey)) || []; } catch (error) { saved = []; }
    function store() {
      try { window.localStorage.setItem(storageKey, JSON.stringify(saved)); } catch (error) { /* unavailable */ }
    }
    let groups = [];
    function rebuildGroups() {
      groups = builtin.concat(saved.map(g => ({kind: "Saved", name: g.name,
                                                members: g.members.filter(id => id in nodes)})))
        .filter(g => g.members.length);
      groups.forEach(g => { g.key = g.kind + ":" + g.name; });
    }
    rebuildGroups();

    // Selected parts are tinted a calm blue, as CAD tools mark a selection, and
    // a part hovered in a list previews it more lightly. The blue is mixed in
    // rather than painted over, so a field's colours still show through: the
    // lit surface is scaled down and the blue added, unlit, in its place.
    const tint = new THREE.Color(0x4a90d9);
    let chosen = new Set(), hovered = new Set();
    function shade(part) {
      const mix = chosen.has(part.mesh) ? 0.55 : hovered.has(part.mesh) ? 0.35 : 0;
      part.mesh.material.color.setRGB(...part.surface).multiplyScalar(1 - mix);
      part.mesh.material.emissive.copy(tint).multiplyScalar(mix);
    }
    function reshade(before, after) {
      new Set([...before, ...after]).forEach(mesh => shade(mesh.userData.part));
    }
    function highlight(ids) {
      const before = hovered;
      hovered = meshesOf(ids);
      reshade(before, hovered);
    }

    const legend = document.getElementById("legend");
    const caption = document.createElement("div");
    const gradient = document.createElement("div");
    gradient.id = "gradient";
    gradient.style.background = "linear-gradient(to right," + data.color_table.map(c => "rgb("+c.join(",")+")").join(",") + ")";
    const limits = document.createElement("div");
    limits.id = "limits";
    const low = document.createElement("span");
    const high = document.createElement("span");
    limits.append(low, high);
    const missing = document.createElement("div");
    missing.textContent = "Grey: unavailable values";
    const selector = document.createElement("select");
    selector.id = "field";
    const names = Object.keys(data.styles);
    [["", "None (component colors)"]].concat(names.map(name => [name, data.styles[name].label]))
      .forEach(([value, text]) => {
        const option = document.createElement("option");
        option.value = value;
        option.textContent = text;
        selector.appendChild(option);
      });
    selector.value = data.style ? data.style.name : "";
    // The range a group sets covers only its members; everything else fades so
    // the parts the scale describes stand out.
    const range = document.createElement("select");
    range.id = "range";
    range.title = "Colour range and focus";
    function fillRange() {
      const current = range.value;
      range.textContent = "";
      [["", "Range: whole scene"]].concat(groups.map(g => [g.key, "Range: " + g.kind + " · " + g.name]))
        .forEach(([value, text]) => {
          const option = document.createElement("option");
          option.value = value;
          option.textContent = text;
          range.appendChild(option);
        });
      range.value = groups.some(g => g.key === current) ? current : "";
    }
    fillRange();
    legend.append(selector, range, caption, gradient, limits, missing);
    function showLegend() {
      legend.hidden = names.length === 0 && groups.length === 0;
      range.hidden = groups.length === 0;
    }
    showLegend();

    function applyColors() {
      const name = selector.value;
      const style = name ? data.styles[name] : null;
      const focus = groups.find(g => g.key === range.value) || null;
      const members = focus ? meshesOf(focus.members) : null;
      let span = style ? style.limits : null;
      let found = Boolean(style && style.has_values);
      if (style && focus) {
        let lo = Infinity, hi = -Infinity;
        parts.forEach(part => {
          const field = members.has(part.mesh) && shown(part.mesh) && part.node.fields[name];
          if (field) field.values.forEach(v => { if (v !== null) { lo = Math.min(lo, v); hi = Math.max(hi, v); } });
        });
        found = lo <= hi;
        if (found && lo === hi) {
          const pad = Math.max(Math.abs(lo)*1e-6, 1e-12);
          lo -= pad;
          hi += pad;
        }
        span = found ? [lo, hi] : [0, 1];
      }
      parts.forEach(part => {
        const inFocus = !members || members.has(part.mesh);
        const field = inFocus && style && part.node.fields[name];
        const material = part.mesh.material;
        if (!field) {
          part.mesh.geometry = part.base;
          material.vertexColors = false;
          part.surface = part.node.color;
        } else {
          let g = part.variants[field.association];
          if (!g) {
            g = field.association === "cell" ? part.base.toNonIndexed() : part.base.clone();
            part.variants[field.association] = g;
          }
          const colors = [];
          field.values.forEach(value => {
            const rgb = color(value, span);
            for (let i = 0; i < (field.association === "cell" ? 3 : 1); ++i) colors.push(...rgb);
          });
          g.setAttribute("color", new THREE.Float32BufferAttribute(colors, 3));
          part.mesh.geometry = g;
          material.vertexColors = true;
          part.surface = [1, 1, 1];
        }
        shade(part);
        const opacity = inFocus ? part.node.opacity : Math.min(part.node.opacity, 0.12);
        material.opacity = opacity;
        material.transparent = opacity < 1;
        // Faded parts must not hide the focused ones behind them.
        material.depthWrite = inFocus;
        material.needsUpdate = true;
      });
      caption.textContent = (style ? style.label + " [" + style.units + "]" : "Component colors")
        + (focus ? " · " + focus.name : "")
        + (style && !found ? " — no finite values" : "");
      const scale = Boolean(style && found);
      gradient.hidden = !scale;
      limits.hidden = !scale;
      missing.hidden = !scale;
      if (scale) {
        low.textContent = span[0].toPrecision(5);
        high.textContent = span[1].toPrecision(5);
      }
    }
    selector.addEventListener("change", applyColors);
    range.addEventListener("change", applyColors);
    applyColors();

    const box = new THREE.Box3().setFromObject(root);
    if (!box.isEmpty()) {
      const center = box.getCenter(new THREE.Vector3());
      const size = box.getSize(new THREE.Vector3());
      const halfFov = Math.min(camera.fov*Math.PI/360, Math.atan(Math.tan(camera.fov*Math.PI/360)*camera.aspect));
      const dist = Math.max(size.length()/2, 0.001)/Math.sin(halfFov)*1.1;
      camera.position.copy(center).add(new THREE.Vector3(1, 0.5, 1).normalize().multiplyScalar(dist));
      camera.near = dist/1000;
      camera.far = dist*1000;
      camera.updateProjectionMatrix();
      controls.target.copy(center);
    } else camera.position.set(2, 1, 2);
    controls.position0.copy(camera.position);
    controls.target0.copy(controls.target);
    controls.up0.copy(camera.up);

    // Corner orientation triad. It renders in its own viewport with the main
    // camera's rotation, so it never occludes geometry or scales with zoom.
    let gizmo = null, gizmoCamera = null;
    if (data.show_axes) {
      gizmo = new THREE.Scene();
      gizmo.add(new THREE.AxesHelper(1));
      gizmoCamera = new THREE.OrthographicCamera(-1.7, 1.7, 1.7, -1.7, 0.1, 10);
      ["x", "y", "z"].forEach((label, i) => {
        const canvas = document.createElement("canvas");
        canvas.width = canvas.height = 64;
        const ctx = canvas.getContext("2d");
        ctx.font = "bold 44px system-ui";
        ctx.textAlign = "center";
        ctx.textBaseline = "middle";
        ctx.fillStyle = ["#d00000", "#00a000", "#0000d0"][i];
        ctx.fillText(label, 32, 32);
        const sprite = new THREE.Sprite(new THREE.SpriteMaterial({
          map: new THREE.CanvasTexture(canvas), transparent: true, depthTest: false
        }));
        sprite.scale.setScalar(0.6);
        sprite.position.setComponent(i, 1.3);
        gizmo.add(sprite);
      });
    }
    function render() {
      const width = window.innerWidth, height = window.innerHeight;
      renderer.autoClear = true;
      renderer.setScissorTest(false);
      renderer.setViewport(0, 0, width, height);
      renderer.render(scene, camera);
      if (!gizmo) return;
      const size = Math.max(70, Math.min(130, Math.min(width, height)*0.2));
      const x = width - size - 12;
      renderer.autoClear = false;
      renderer.setViewport(x, 12, size, size);
      renderer.setScissor(x, 12, size, size);
      renderer.setScissorTest(true);
      renderer.clearDepth();
      gizmoCamera.position.subVectors(camera.position, controls.target).normalize().multiplyScalar(4);
      gizmoCamera.up.copy(camera.up);
      gizmoCamera.lookAt(0, 0, 0);
      renderer.render(gizmo, gizmoCamera);
    }

    const info = document.getElementById("info");
    const raycaster = new THREE.Raycaster();
    const pointer = new THREE.Vector2();

    function number(value) {
      // Round-trippable but not noisy: Python accepts either form verbatim.
      return String(Number(value.toPrecision(9)));
    }
    function triple(values) {
      return "(" + values.map(number).join(", ") + ")";
    }
    function row(label, text, mono) {
      return element("div", mono ? "mono" : "", label ? label + " " + text : text);
    }

    // The selection is a list of node IDs in the order they were picked.
    let selected = [];
    function select(ids, at) {
      selected = [...new Set(ids)];
      const before = chosen;
      chosen = meshesOf(selected);
      reshade(before, chosen);
      if (selected.length) document.body.dataset.selected = selected.join("\n");
      else delete document.body.dataset.selected;
      showInfo(at);
    }
    function toggle(ids, at) {
      const remove = ids.every(id => selected.includes(id));
      select(remove ? selected.filter(id => !ids.includes(id)) : selected.concat(ids), at);
    }
    function deselect() {
      select([]);
    }

    function details(group) {
      const node = group.userData.node;
      const body = [];
      const kind = node.metadata && node.metadata.kind;
      if (kind) body.push(row("", kind.replace(/_/g, " ")));
      const fluid = fluidOf(group);
      if (fluid) body.push(row("Fluid:", fluid));
      if (node.tags.length) body.push(row("Tags:", node.tags.join(", ")));
      if (node.placement) {
        body.push(row("position=", triple(node.placement.position), true));
        body.push(row("rotation=", triple(node.placement.rotation), true));
      } else body.push(row("Part of", pathName(group.parent.name)));
      // Fields usually sit on the channels or blades below a component.
      const fields = new Set();
      (function visit(n) { Object.keys(n.fields).forEach(f => fields.add(f)); n.children.forEach(visit); })(node);
      if (fields.size) body.push(row("Fields:", [...fields].join(", ")));
      const anchors = Object.keys(node.anchors);
      if (anchors.length) {
        body.push(row("", "Anchors (scene coordinates)"));
        anchors.forEach(name => {
          // The anchor's own translation carried through the node's placement.
          const local = new THREE.Matrix4().fromArray(node.anchors[name]);
          const at = new THREE.Vector3().setFromMatrixPosition(local).applyMatrix4(group.matrixWorld);
          body.push(row(name, triple([at.x, at.y, at.z]), true));
        });
      }
      return body;
    }
    function listing() {
      return selected.map(id => {
        const line = element("div", "item");
        line.append(element("span", "", pathName(id)));
        const drop = element("button", "drop", "×");
        drop.title = "Remove from selection";
        drop.addEventListener("click", () => { highlight([]); toggle([id]); });
        line.append(drop);
        line.addEventListener("mouseenter", () => highlight([id]));
        line.addEventListener("mouseleave", () => highlight([]));
        return line;
      });
    }
    function creator() {
      const form = element("form", "");
      form.id = "create";
      const input = element("input", "");
      input.placeholder = "Group name";
      input.required = true;
      const button = element("button", "", "Create group");
      button.type = "submit";
      form.append(input, button);
      form.addEventListener("submit", event => {
        event.preventDefault();
        const name = input.value.trim();
        if (!name) return;
        // Saving under an existing name replaces that group.
        saved = saved.filter(g => g.name !== name).concat([{name: name, members: selected.slice()}]);
        store();
        rebuildGroups();
        fillRange();
        showLegend();
        renderGroups();
        input.value = "";
      });
      return form;
    }
    // Place the box beside the event that opened it. Once it is showing it
    // stays put, whether it was dragged there or not, so adding parts to the
    // selection does not chase the cursor.
    function showInfo(at) {
      if (!selected.length) {
        info.hidden = true;
        return;
      }
      const opening = info.hidden;
      info.textContent = "";
      const head = element("div", "");
      head.id = "head";
      const heading = element("span", "");
      heading.id = "title";
      heading.textContent = selected.length === 1 ? title(nodes[selected[0]].userData.node)
        : selected.length + " parts selected";
      const close = element("button", "", "×");
      close.id = "close";
      close.title = "Close and deselect";
      close.addEventListener("click", deselect);
      head.append(heading, close);
      const body = element("div", "");
      body.id = "body";
      body.append(...(selected.length === 1 ? details(nodes[selected[0]]) : listing()), creator());
      info.append(head, body);
      info.hidden = false;
      const margin = 12;
      if (opening) {
        info.style.left = "0px";
        info.style.top = "0px";
        const size = info.getBoundingClientRect();
        const x = at ? at.clientX + 14 : margin;
        const y = at ? at.clientY + 14 : 48;
        info.style.left = Math.max(margin, Math.min(x, window.innerWidth - size.width - margin)) + "px";
        info.style.top = Math.max(margin, Math.min(y, window.innerHeight - size.height - margin)) + "px";
      } else keepInside();
    }
    function keepInside() {
      const size = info.getBoundingClientRect();
      info.style.left = Math.max(0, Math.min(size.left, window.innerWidth - size.width)) + "px";
      info.style.top = Math.max(0, Math.min(size.top, window.innerHeight - size.height)) + "px";
    }
    // Drag the box by its heading.
    info.addEventListener("pointerdown", event => {
      const head = event.target.closest("#head");
      if (!head || event.target.closest("button") || event.button !== 0) return;
      const start = info.getBoundingClientRect();
      const dx = event.clientX - start.left, dy = event.clientY - start.top;
      head.setPointerCapture(event.pointerId);
      const move = e => {
        info.style.left = (e.clientX - dx) + "px";
        info.style.top = (e.clientY - dy) + "px";
        keepInside();
      };
      const end = () => {
        head.removeEventListener("pointermove", move);
        head.removeEventListener("pointerup", end);
      };
      head.addEventListener("pointermove", move);
      head.addEventListener("pointerup", end);
      event.preventDefault();
    });

    function hits(event) {
      const rect = renderer.domElement.getBoundingClientRect();
      pointer.x = ((event.clientX - rect.left)/rect.width)*2 - 1;
      pointer.y = -((event.clientY - rect.top)/rect.height)*2 + 1;
      raycaster.setFromCamera(pointer, camera);
      return raycaster.intersectObjects(root.children, true).filter(hit => shown(hit.object));
    }
    function click(event) {
      const hit = hits(event)[0];
      if (!hit) {
        if (!event.shiftKey) deselect();
        return;
      }
      // Climb to the child of the root: the node a single Scene.add() created.
      const id = topOf(hit.object).name;
      if (event.shiftKey) toggle([id], event);
      else select([id], event);
    }

    // Right-click lists every part under the cursor: each component hit, front
    // to back, with its circuits and other named parts. Repeated instances
    // (channels, blades, vanes) are listed only where the cursor meets them.
    const pick = document.getElementById("pick");
    function closePick() {
      pick.hidden = true;
      highlight([]);
    }
    function openPick(event) {
      const found = hits(event);
      if (!found.length) return closePick();
      const reached = new Set();
      found.forEach(hit => { for (let o = hit.object; o && o !== root; o = o.parent) reached.add(o.name); });
      const tops = [...new Set(found.map(hit => topOf(hit.object)))];
      pick.textContent = "";
      pick.append(element("div", "heading", "Parts under the cursor"));
      function visit(group, depth) {
        const node = group.userData.node;
        if (!shown(group)) return;
        if (depth && node.metadata && "instance_index" in node.metadata && !reached.has(node.id)) return;
        const line = element("div", "item", title(node));
        line.style.paddingLeft = (6 + 14*depth) + "px";
        const fluid = node.metadata && node.metadata.fluid;
        if (fluid) line.append(element("span", "note", " " + fluid));
        line.addEventListener("mouseenter", () => highlight([node.id]));
        line.addEventListener("mouseleave", () => highlight([]));
        line.addEventListener("click", e => {
          closePick();
          if (e.shiftKey) toggle([node.id], event);
          else select([node.id], event);
        });
        pick.append(line);
        group.children.filter(child => child.isGroup).forEach(child => visit(child, depth + 1));
      }
      tops.forEach(top => visit(top, 0));
      pick.hidden = false;
      pick.style.left = "0px";
      pick.style.top = "0px";
      const size = pick.getBoundingClientRect();
      pick.style.left = Math.max(4, Math.min(event.clientX + 4, window.innerWidth - size.width - 4)) + "px";
      pick.style.top = Math.max(4, Math.min(event.clientY + 4, window.innerHeight - size.height - 4)) + "px";
    }
    window.addEventListener("pointerdown", event => {
      if (!pick.hidden && !pick.contains(event.target)) closePick();
    }, {capture: true});

    // Groups panel: clicking a group selects its members, Shift+click adds them.
    const panel = document.getElementById("groups");
    function renderGroups() {
      panel.textContent = "";
      panel.hidden = groups.length === 0;
      panel.append(element("div", "heading", "Groups"));
      groups.forEach(g => {
        const line = element("div", "item");
        line.append(element("span", "note", g.kind + " "), element("span", "", g.name),
                    element("span", "note", " (" + g.members.length + ")"));
        line.title = g.members.map(pathName).join("\n");
        line.addEventListener("mouseenter", () => highlight(g.members));
        line.addEventListener("mouseleave", () => highlight([]));
        line.addEventListener("click", event => {
          if (event.target.closest("button")) return;
          if (event.shiftKey) select(selected.concat(g.members));
          else select(g.members);
        });
        if (g.kind === "Saved") {
          const drop = element("button", "drop", "×");
          drop.title = "Delete group";
          drop.addEventListener("click", () => {
            saved = saved.filter(s => s.name !== g.name);
            store();
            highlight([]);
            rebuildGroups();
            fillRange();
            showLegend();
            renderGroups();
            applyColors();
          });
          line.append(drop);
        }
        panel.append(line);
      });
    }
    renderGroups();

    // Distinguish clicks from trackball drags, which share the same buttons.
    let down = null;
    renderer.domElement.addEventListener("pointerdown", event => {
      down = event.button === 0 || event.button === 2 ? {x: event.clientX, y: event.clientY, button: event.button} : null;
    });
    renderer.domElement.addEventListener("pointerup", event => {
      const still = down && event.button === down.button && Math.hypot(event.clientX - down.x, event.clientY - down.y) < 5;
      if (still && event.button === 0) click(event);
      if (still && event.button === 2) openPick(event);
      down = null;
    });

    controls.update();
    render();
    hud.textContent = "Drag: rotate • Right-drag or Ctrl+drag: pan • Wheel: zoom • Click: inspect • "
      + "Shift+click: add • Right-click: parts under cursor • R: reset view";
    window.addEventListener("keydown", event => {
      if (event.target.closest && event.target.closest("input, select, textarea")) return;
      if (event.key === "r" || event.key === "R") controls.reset();
      if (event.key === "Escape") {
        if (!pick.hidden) closePick();
        else deselect();
      }
    });
    document.body.dataset.sceneReady = "true";
    window.addEventListener("resize", () => {
      renderer.setSize(window.innerWidth, window.innerHeight, false);
      camera.aspect = window.innerWidth/window.innerHeight;
      camera.updateProjectionMatrix();
      controls.handleResize();
      if (!info.hidden) keepInside();
    });
    function animate() {
      controls.update();
      render();
      requestAnimationFrame(animate);
    }
    animate();
  } catch (error) {
    hud.textContent = "Failed to display scene: " + error.message;
    document.body.dataset.sceneError = error.message;
    console.error(error);
  }
})();
