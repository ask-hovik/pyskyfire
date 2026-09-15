/* Offline scene display. Clicking a part selects the top-level component that
 * owns it, so one cooling channel identifies its whole thrust chamber, and
 * reports the placement arguments that reproduce its position in Python. */
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
    // be swapped without rebuilding the hierarchy.
    const parts = [];
    function add(node, parent) {
      const group = new THREE.Group();
      group.name = node.id;
      group.userData = {node: node, fields: node.fields, metadata: node.metadata, anchors: node.anchors};
      group.matrix.fromArray(node.transform);
      group.matrixAutoUpdate = false;
      group.visible = node.visible;
      parent.add(group);
      if (node.mesh !== null) {
        const material = new THREE.MeshPhongMaterial({
          side: THREE.DoubleSide, transparent: node.opacity < 1, opacity: node.opacity
        });
        const mesh = new THREE.Mesh(geometry[node.mesh], material);
        mesh.name = node.id;
        group.add(mesh);
        // Seed geometry is shared between repeated channels, so colored copies
        // are created per part and cached by association rather than mutated.
        parts.push({node: node, mesh: mesh, base: geometry[node.mesh], variants: {}});
      }
      node.children.forEach(child => add(child, group));
    }
    const root = new THREE.Group();
    scene.add(root);
    data.nodes.forEach(node => add(node, root));
    root.updateMatrixWorld(true);

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
    legend.append(selector, caption, gradient, limits, missing);
    legend.hidden = names.length === 0;

    function applyField(name) {
      const style = name ? data.styles[name] : null;
      parts.forEach(part => {
        const field = style && part.node.fields[name];
        const material = part.mesh.material;
        if (!field) {
          part.mesh.geometry = part.base;
          material.vertexColors = false;
          material.color.setRGB(...part.node.color);
        } else {
          let g = part.variants[field.association];
          if (!g) {
            g = field.association === "cell" ? part.base.toNonIndexed() : part.base.clone();
            part.variants[field.association] = g;
          }
          const colors = [];
          field.values.forEach(value => {
            const rgb = color(value, style.limits);
            for (let i = 0; i < (field.association === "cell" ? 3 : 1); ++i) colors.push(...rgb);
          });
          g.setAttribute("color", new THREE.Float32BufferAttribute(colors, 3));
          part.mesh.geometry = g;
          material.vertexColors = true;
          material.color.setRGB(1, 1, 1);
        }
        material.needsUpdate = true;
      });
      caption.textContent = style
        ? style.label + " [" + style.units + "]" + (style.has_values ? "" : " — no finite values")
        : "Component colors";
      const scale = Boolean(style && style.has_values);
      gradient.hidden = !scale;
      limits.hidden = !scale;
      missing.hidden = !scale;
      if (scale) {
        low.textContent = style.limits[0].toPrecision(5);
        high.textContent = style.limits[1].toPrecision(5);
      }
    }
    selector.addEventListener("change", () => applyField(selector.value));
    applyField(selector.value);

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
    const outline = new THREE.Box3Helper(new THREE.Box3(), new THREE.Color(0x000000));
    outline.visible = false;
    scene.add(outline);
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
      const line = document.createElement("div");
      line.className = mono ? "mono" : "";
      line.textContent = label ? label + " " + text : text;
      return line;
    }
    function deselect() {
      info.hidden = true;
      outline.visible = false;
      delete document.body.dataset.selected;
    }
    function describe(group, event) {
      const node = group.userData.node;
      info.textContent = "";
      const close = document.createElement("button");
      close.id = "close";
      close.textContent = "×";
      close.title = "Close";
      close.addEventListener("click", deselect);
      info.appendChild(close);
      const heading = document.createElement("div");
      heading.id = "title";
      heading.textContent = node.name;
      info.appendChild(heading);
      const kind = node.metadata && node.metadata.kind;
      if (kind) info.appendChild(row("", kind.replace(/_/g, " ")));
      info.appendChild(row("position=", triple(node.placement.position), true));
      info.appendChild(row("rotation=", triple(node.placement.rotation), true));
      const fields = Object.keys(node.fields);
      if (fields.length) info.appendChild(row("Fields:", fields.join(", ")));
      const anchors = Object.keys(node.anchors);
      if (anchors.length) {
        info.appendChild(row("", "Anchors (scene coordinates)"));
        anchors.forEach(name => {
          // The anchor's own translation carried through the node's placement.
          const local = new THREE.Matrix4().fromArray(node.anchors[name]);
          const at = new THREE.Vector3().setFromMatrixPosition(local).applyMatrix4(group.matrixWorld);
          info.appendChild(row(name, triple([at.x, at.y, at.z]), true));
        });
      }
      info.hidden = false;
      // Place beside the click, then pull back inside the window if it overflows.
      const margin = 12;
      info.style.left = "0px";
      info.style.top = "0px";
      const size = info.getBoundingClientRect();
      const x = Math.max(margin, Math.min(event.clientX + 14, window.innerWidth - size.width - margin));
      const y = Math.max(margin, Math.min(event.clientY + 14, window.innerHeight - size.height - margin));
      info.style.left = x + "px";
      info.style.top = y + "px";
      outline.box.setFromObject(group);
      outline.visible = !outline.box.isEmpty();
      document.body.dataset.selected = node.id;
    }
    function select(event) {
      const rect = renderer.domElement.getBoundingClientRect();
      pointer.x = ((event.clientX - rect.left)/rect.width)*2 - 1;
      pointer.y = -((event.clientY - rect.top)/rect.height)*2 + 1;
      raycaster.setFromCamera(pointer, camera);
      const hit = raycaster.intersectObjects(root.children, true)[0];
      if (!hit) return deselect();
      // Climb to the child of the root: the node a single Scene.add() created.
      let group = hit.object;
      while (group.parent && group.parent !== root) group = group.parent;
      describe(group, event);
    }
    // Distinguish a click from a trackball drag, which shares the same button.
    let down = null;
    renderer.domElement.addEventListener("pointerdown", event => {
      down = event.button === 0 ? {x: event.clientX, y: event.clientY} : null;
    });
    renderer.domElement.addEventListener("pointerup", event => {
      if (down && Math.hypot(event.clientX - down.x, event.clientY - down.y) < 5) select(event);
      down = null;
    });

    controls.update();
    render();
    hud.textContent = "Drag: rotate • Right-drag: pan • Wheel: zoom • Click: inspect • R: reset view";
    window.addEventListener("keydown", event => {
      if (event.key === "r" || event.key === "R") controls.reset();
      if (event.key === "Escape") deselect();
    });
    document.body.dataset.sceneReady = "true";
    window.addEventListener("resize", () => {
      renderer.setSize(window.innerWidth, window.innerHeight, false);
      camera.aspect = window.innerWidth/window.innerHeight;
      camera.updateProjectionMatrix();
      controls.handleResize();
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
