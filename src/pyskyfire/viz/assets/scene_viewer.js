/* Offline scene display. Node identity and scalar values remain available for
 * future selection controls; this viewer intentionally has no picking UI. */
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
    const controls = new THREE.OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;
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
    function color(value) {
      if (value === null) return [0.5, 0.5, 0.5];
      const [lo, hi] = data.style.limits;
      const index = Math.max(0, Math.min(255, Math.floor(256*(value-lo)/(hi-lo))));
      return data.color_table[index].map(x => x/255);
    }
    function add(node, parent) {
      const group = new THREE.Group();
      group.name = node.id;
      group.userData = {fields: node.fields, metadata: node.metadata, anchors: node.anchors};
      group.matrix.fromArray(node.transform);
      group.matrixAutoUpdate = false;
      group.visible = node.visible;
      parent.add(group);
      if (node.mesh !== null) {
        let g = geometry[node.mesh];
        const field = data.style && node.fields[data.style.name];
        if (field) {
          g = field.association === "cell" ? g.toNonIndexed() : g.clone();
          const colors = [];
          field.values.forEach(value => {
            const rgb = color(value);
            for (let i=0; i<(field.association === "cell" ? 3 : 1); ++i) colors.push(...rgb);
          });
          g.setAttribute("color", new THREE.Float32BufferAttribute(colors, 3));
        }
        const material = new THREE.MeshPhongMaterial({
          color: field ? 0xffffff : new THREE.Color(...node.color),
          vertexColors: Boolean(field), side: THREE.DoubleSide,
          transparent: node.opacity < 1, opacity: node.opacity
        });
        const mesh = new THREE.Mesh(g, material);
        mesh.name = node.id;
        group.add(mesh);
      }
      node.children.forEach(child => add(child, group));
    }
    const root = new THREE.Group();
    scene.add(root);
    data.nodes.forEach(node => add(node, root));
    root.updateMatrixWorld(true);
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
    if (data.style) {
      const legend = document.getElementById("legend");
      legend.hidden = false;
      legend.textContent = data.style.label + " [" + data.style.units + "]" +
        (data.style.has_values ? "" : " — no finite values");
      const gradient = document.createElement("div");
      gradient.id = "gradient";
      gradient.style.background = "linear-gradient(to right," + data.color_table.map(c => "rgb("+c.join(",")+")").join(",") + ")";
      legend.appendChild(gradient);
      const limits = document.createElement("div");
      limits.id = "limits";
      data.style.limits.forEach(value => {
        const span = document.createElement("span");
        span.textContent = value.toPrecision(5);
        limits.appendChild(span);
      });
      legend.appendChild(limits);
      if (!data.style.has_values) {
        gradient.hidden = true;
        limits.hidden = true;
      }
      const missing = document.createElement("div");
      missing.textContent = "Grey: unavailable values";
      legend.appendChild(missing);
    }
    controls.update();
    renderer.render(scene, camera);
    hud.textContent = "Drag: rotate • Right-drag: pan • Wheel: zoom";
    document.body.dataset.sceneReady = "true";
    window.addEventListener("resize", () => {
      renderer.setSize(window.innerWidth, window.innerHeight, false);
      camera.aspect = window.innerWidth/window.innerHeight;
      camera.updateProjectionMatrix();
    });
    function animate() {
      controls.update();
      renderer.render(scene, camera);
      requestAnimationFrame(animate);
    }
    animate();
  } catch (error) {
    hud.textContent = "Failed to display scene: " + error.message;
    document.body.dataset.sceneError = error.message;
    console.error(error);
  }
})();
