
import numpy as np
from vispy import app, scene
from vispy.scene import visuals, SceneCanvas
from vispy.visuals.filters.clipping_planes import PlanesClipper
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import json

    # Professional style + Latin Modern Roman 10pt
plt.style.use('ggplot')
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Latin Modern Roman'],
    'font.size': 10,
    'mathtext.fontset': 'cm',                 # ← built-in Computer Modern
})

def plot_engine_3D(
    thrust_chamber,
    contour=None,
    fraction=1.0,
    clip_plane=None,
    draw_section_line=False,
    line_color=(1.0, 0.0, 0.0, 1.0),
    draw_only_section_line=False, 
    show_axis=True
):
    """
    Visualizes the cooling channels defined in the thrust chamber using VisPy,
    optionally applying an arbitrary clipping plane, drawing the cross-section line,
    customizing line color, drawing only the intersection line, and overlaying
    a contour line for inspection. Supports multiple walls per channel if
    compute_point_cloud returns them.

    Parameters
    ----------
    thrust_chamber : ThrustChamber
        A thrust chamber object with cooling circuit geometry.
    contour : Contour or None
        If provided, an instance of Contour whose inner wall will be drawn at θ=0.
    fraction : float, optional
        Fraction of channels to plot per circuit (0 to 1).
    clip_plane : tuple or None
        If provided, should be a tuple ((px, py, pz), (nx, ny, nz)), where:
        - (px, py, pz) is a point on the plane
        - (nx, ny, nz) is the plane normal vector
    draw_section_line : bool, optional
        If True, compute and draw the intersection line of the mesh
        with the clipping plane.
    line_color : tuple, optional
        RGBA color for the section line, e.g. (1,0,0,1).
    draw_only_section_line : bool, optional
        If True, skip drawing the full mesh and render only the intersection line
        with default camera facing perpendicular to slice and orthographic projection.

    Returns
    -------
    canvas : vispy.scene.SceneCanvas
        The VisPy canvas for interactive viewing.
    """
    # Build channel geometry
    thrust_chamber.build_channel_centerlines(mode="plot")
    for circuit in thrust_chamber.cooling_circuit_group.circuits:
        circuit.finalize()
        circuit.compute_geometry()

    # Prepare clip plane
    pt = None
    norm_arr = None
    if clip_plane is not None:
        pt, norm = clip_plane
        norm_arr = np.array(norm, dtype=np.float32)

    # Set up VisPy scene
    canvas = SceneCanvas(keys="interactive", show=True, bgcolor="white", size=(4000, 5000))
    view = canvas.central_widget.add_view()
    if show_axis:
        visuals.XYZAxis(parent=view.scene)


    # Draw contour line at θ=0 for inspection
    if contour is not None:
        xs = np.array(contour.xs, dtype=np.float32)
        rs = np.array(contour.rs, dtype=np.float32)
        contour_pts = np.column_stack((xs, rs, np.zeros_like(xs)))
        visuals.Line(
            pos=contour_pts,
            connect='strip',
            color=(1.0, 0.0, 0.0, 1.0),
            width=2,
            parent=view.scene
        )

    # Draw clip plane normal arrow
    if clip_plane is not None:
        start = np.array(pt, dtype=np.float32)
        end = start + norm_arr
        visuals.Arrow(
            pos=np.vstack([start, end]),
            arrow_size=10,
            color=(0.0, 0.0, 0.0, 1.0),
            parent=view.scene
        )

    # Choose camera
    if draw_only_section_line and clip_plane is not None:
        norm_u = norm_arr / np.linalg.norm(norm_arr)
        nx, ny, nz = norm_u
        elevation = np.degrees(np.arcsin(nz))
        azimuth = -np.degrees(np.arctan2(nx, ny)) + 180
        cam = scene.cameras.TurntableCamera(
            elevation=elevation,
            azimuth=azimuth,
            fov=0.0,
            up="z"
        )
    else:
        cam = scene.cameras.TurntableCamera(
            fov=45.0,
            azimuth=30.0,
            elevation=30.0
        )
    view.camera = cam

    # Collect vertices/faces, supporting multiple walls per channel
    all_verts = []
    all_faces = []
    vert_offset = 0
    for circuit in thrust_chamber.cooling_circuit_group.circuits:
        n_channels = circuit.placement.n_channel_positions * circuit.placement.n_channels_per_leaf
        n_plot = int(round(fraction * n_channels))
        for ch_idx in range(max(0, min(n_plot, n_channels))):
            cyl = circuit.centerlines[ch_idx]
            pc = circuit.point_cloud[ch_idx]
            # Detect multiple wall point-clouds (tuple/list or array with ndim=4)
            if isinstance(pc, (tuple, list)) or (isinstance(pc, np.ndarray) and pc.ndim == 4):
                pcs = pc
            else:
                pcs = [pc]

            x, r, theta = cyl.T
            center_cart = np.column_stack((x, r*np.cos(theta), r*np.sin(theta)))

            for wall_pc in pcs:
                N, M, _ = wall_pc.shape
                # build sections
                sections = np.empty((N, M, 3), dtype=np.float32)
                for i in range(N):
                    # flip sign for outward
                    sections[i] = center_cart[i] - wall_pc[i]
                # close loop
                sec_closed = np.concatenate([sections, sections[:, :1, :]], axis=1)
                Mcl = M + 1
                verts = sec_closed.reshape(-1, 3)
                faces = []
                for ii in range(N - 1):
                    for jj in range(Mcl - 1):
                        v0 = vert_offset + ii*Mcl + jj
                        v1 = v0 + 1
                        v2 = vert_offset + (ii+1)*Mcl + jj
                        v3 = v2 + 1
                        faces.append([v0, v1, v2])
                        faces.append([v2, v1, v3])
                all_verts.append(verts)
                all_faces.append(np.array(faces, dtype=np.uint32))
                vert_offset += verts.shape[0]

    if not all_verts:
        raise RuntimeError("No cooling channels to visualize.")
    verts_all = np.vstack(all_verts)
    faces_all = np.vstack(all_faces)

    # Compute section segments if requested
    section_segments = None
    if clip_plane is not None and draw_section_line:
        dists = np.dot(verts_all - pt, norm_arr)
        segments = []
        for face in faces_all:
            ds = dists[face]
            vs = verts_all[face]
            ips = []
            for i in range(3):
                d0, d1 = ds[i], ds[(i+1) % 3]
                if d0 * d1 < 0:
                    t = d0 / (d0 - d1)
                    p = vs[i] + t * (vs[(i+1) % 3] - vs[i])
                    ips.append(p)
            if len(ips) == 2:
                segments.extend(ips)
        if segments:
            section_segments = np.array(segments, dtype=np.float32)

    if draw_only_section_line:
        if section_segments is None:
            raise RuntimeError("No intersection line: check clip_plane and draw_section_line flags.")
        visuals.Line(
            pos=section_segments,
            connect='segments',
            color=line_color,
            width=2,
            parent=view.scene
        )
        app.run()
        return canvas

    mesh = visuals.Mesh(
        vertices=verts_all,
        faces=faces_all,
        color=(0.5, 0.8, 0.9, 1.0),
        shading="smooth",
        parent=view.scene
    )
    if clip_plane is not None:
        clipper = PlanesClipper()
        mesh.attach(clipper)
        clipper.clipping_planes = np.array([[pt, norm_arr]], dtype=np.float32)

    if section_segments is not None:
        visuals.Line(
            pos=section_segments,
            connect='segments',
            color=line_color,
            width=2,
            parent=view.scene
        )

    app.run()
    return canvas, view

from vispy.io import write_png


def save_canvas_png(canvas,
                    filename="engine_4k.png",
                    out_px=(4000, 3000),   # final physical resolution
                    bg=(1, 1, 1, 1)):
    """
    Export *canvas* to a PNG that really fills the whole image.
    Handles Hi-DPI screens automatically.
    """
    # --- figure out how many *logical* pixels we need -------------
    #scale   = canvas.pixel_scale or 1          # 2 on Retina / HiDPI, else 1
    #logical = (int(out_px[0] / scale),
    #           int(out_px[1] / scale))

    # --- temporarily resize the canvas ---------------------------
    #orig_size = canvas.size          # save to restore later
    #canvas.size = logical            # triggers a full resize event
    canvas.update()
    canvas.app.process_events()      # make sure transforms catch up

    # --- off-screen render & write -------------------------------
    img = canvas.render( bgcolor=bg)   # fills whole FBO now
    write_png(filename, img)

    # --- restore original on-screen size (optional) --------------
    #canvas.size = orig_size



def plot_precomputed_properties(thrust_chamber, circuit_index=0):
    """
    Plot all precomputed thermal properties for a chosen cooling circuit.

    Parameters:
        thrust_chamber: ThrustChamber3 object that has precomputed thermal properties.
        circuit_index (int): Index of the cooling circuit to plot (default: 0).

    This function creates a 2x3 grid of subplots showing:
        - dA/dx for the thermal exhaust side,
        - dA/dx for the thermal coolant side,
        - Cooling channel cross-sectional area (A_coolant),
        - dA/dx of the coolant channel area,
        - Coolant hydraulic diameter (Dh_coolant), and
        - Local radius of curvature.
    """
    # Get the chosen circuit
    circuit = thrust_chamber.cooling_circuit_group.circuits[circuit_index]
    # x-domain over which properties are computed
    x = circuit.x_domain

    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    
    # dA/dx Thermal Exhaust
    axs[0, 0].plot(x, circuit.dA_dx_thermal_exhaust_vals, 'b-')
    axs[0, 0].set_title("dA/dx Thermal Hot Gas")
    axs[0, 0].set_xlabel("Axial Position (x)")
    axs[0, 0].set_ylabel("dA/dx")

    # dA/dx Thermal Coolant
    axs[0, 1].plot(x, circuit.dA_dx_thermal_coolant_vals, 'g-')
    axs[0, 1].set_title("dA/dx Thermal Coolant")
    axs[0, 1].set_xlabel("Axial Position (x)")
    axs[0, 1].set_ylabel("dA/dx")

    # Coolant Cross-Sectional Area (A_coolant)
    axs[0, 2].plot(x, circuit.A_coolant_vals, 'r-')
    axs[0, 2].set_title("Coolant Cross-Sectional Area")
    axs[0, 2].set_xlabel("Axial Position (x)")
    axs[0, 2].set_ylabel(r"A ($\mathrm{m^2}$)")

    # dA/dx Coolant (Derivative of A_coolant)
    axs[1, 0].plot(x, circuit.dA_dx_coolant_vals, 'm-')
    axs[1, 0].set_title("dA/dx Coolant")
    axs[1, 0].set_xlabel("Axial Position (x)")
    axs[1, 0].set_ylabel("dA/dx")

    # Hydraulic Diameter (Dh_coolant)
    axs[1, 1].plot(x, circuit.Dh_coolant_vals, 'c-')
    axs[1, 1].set_title("Coolant Hydraulic Diameter ")
    axs[1, 1].set_xlabel("Axial Position (x)")
    axs[1, 1].set_ylabel(r"$\mathrm{D_h}$ (m)")

    # Radius of Curvature
    axs[1, 2].plot(x, circuit.radius_of_curvature_vals, 'k-')
    axs[1, 2].set_title("Radius of Curvature")
    axs[1, 2].set_xlabel("Axial Position (x)")
    axs[1, 2].set_ylabel("R (m)")

    plt.tight_layout()



import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt


def plot_contour(contours,
                 colors=None,
                 linestyles=None,
                 show_labels=True,
                 title="Contour Profiles"):
    """
    Plot bell-nozzle or toroidal-aerospike contours, mirrored about the x-axis.

    Parameters
    ----------
    contours : sequence
        Each element can be either a classic `Contour` (xs/rs) or a
        `ContourToroidalAerospike` (xs_outer / rs_outer / xs_inner / rs_inner).
    colors, linestyles, show_labels, title
        Same meaning as the original routine.
    """
    # ------------------------------------------------------------------ set-up
    plt.style.use("ggplot")
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Latin Modern Roman"],
        "font.size": 10,
        "savefig.dpi": 600
    })

    n = len(contours)
    if colors is None:
        colors = [f"C{i}" for i in range(n)]
    if linestyles is None:
        linestyles = ["-"] + ["--"] * (n - 1)

    if len(colors) < n:
        raise ValueError(f"Need ≥ {n} colours, got {len(colors)}.")
    if len(linestyles) < n:
        raise ValueError(f"Need ≥ {n} linestyles, got {len(linestyles)}.")

    fig, ax = plt.subplots(figsize=(7, 5))

    # ------------------------------------------------------------ main loop
    for idx, contour in enumerate(contours):
        colour = colors[idx]
        lstyle = linestyles[idx]

        # ---------- Classic bell nozzle: xs / rs --------------------------
        if hasattr(contour, "xs") and hasattr(contour, "rs"):
            xs = contour.xs
            rs = contour.rs

            label = getattr(contour, "name", f"Contour {idx + 1}")
            ax.plot(xs,  rs,  color=colour, linestyle=lstyle, linewidth=2, label=label)
            ax.plot(xs, -np.asarray(rs), color=colour, linestyle=lstyle, linewidth=2, label="_nolegend_")

        # ---------- Aerospike: xs_outer / rs_outer (+ inner) --------------
        elif (hasattr(contour, "xs_outer") and hasattr(contour, "rs_outer")
              and hasattr(contour, "xs_inner") and hasattr(contour, "rs_inner")):

            xs_o, rs_o = contour.xs_outer, contour.rs_outer
            xs_i, rs_i = contour.xs_inner, contour.rs_inner

            label = getattr(contour, "name", f"Contour {idx + 1}")
            # outer wall -- in legend
            ax.plot(xs_o,  rs_o,  color=colour, linestyle=lstyle, linewidth=2, label=label)
            ax.plot(xs_o, -rs_o,  color=colour, linestyle=lstyle, linewidth=2, label="_nolegend_")
            # inner wall -- no legend, lighter weight
            ax.plot(xs_i,  rs_i,  color=colour, linestyle=lstyle, linewidth=1.0, label="_nolegend_")
            ax.plot(xs_i, -rs_i,  color=colour, linestyle=lstyle, linewidth=1.0, label="_nolegend_")

        else:
            raise TypeError(f"Object at index {idx} lacks the required contour attributes.")

    # ---------------------------------------------------------- cosmetics
    ax.set_aspect("equal", "box")
    ax.set_xlabel("Axial position, x (m)")
    ax.set_ylabel("Radius, r (m)")
    if title:
        ax.set_title(title, pad=12)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(width=1.2)

    if show_labels:
        lg = ax.legend(frameon=True, framealpha=0.9, edgecolor="black")
        lg.get_frame().set_linewidth(1.0)

    plt.tight_layout()
    return fig, ax



def plot_temperature_profile(results, thrust_chamber, circuit_index, x_query,
                             n_bl=1000):
    """
    Plot T(y) through gas BL, wall layers, coolant BL at axial position x_query.
    """
    


    # 1) find nearest node
    x_arr = np.array(results['x'])
    i = np.argmin(np.abs(x_arr - x_query))
    x0 = x_arr[i]

    # 2) extract temperatures and heat flux
    T_inf = thrust_chamber.combustion_transport.get_T(x0)
    T_hw = results['T'][i, -1]
    T_cw = results['T'][i, 1]
    T_c  = results['T_static'][i]
    qpp  = results['dQ_dA'][i]

    # 3) gas-side BL
    k_g = thrust_chamber.combustion_transport.get_k(x0)
    h_g = qpp / (T_inf - T_hw)
    delta_g = 7 * k_g / h_g
    y_g = np.linspace(-delta_g, 0.0, n_bl)
    T_g = np.where(
        np.abs(y_g) <= delta_g,
        T_inf + (T_hw - T_inf) * (1 - (np.abs(y_g) / delta_g)**(1/7)),
        T_inf
    )

    Ts_rev  = results['T'][i, 1:]      # coolant-side → hot-side
    Ts_wall = Ts_rev[::-1]             # hot-side → coolant-side

    walls      = thrust_chamber.wall_group.walls
    thicknesses = [w.thickness(x0) for w in walls]
    # y‐locations of each interface, starting at the hot wall face
    y_w = np.insert(np.cumsum(thicknesses), 0, 0.0)
    wall_thickness = y_w[-1]

    # 5) coolant-side BL
    p_static = results['p_static'][i]
    T_film = 0.5 * (T_c + T_cw)
    coolant = thrust_chamber.cooling_circuit_group.circuits[circuit_index].coolant_transport
    k_c = coolant.get_k(T_film, p_static)
    h_c = qpp / (T_cw - T_c)
    delta_c = 7 * k_c / h_c
    y_c = np.linspace(0.0, delta_c, n_bl)
    T_cBL = np.where(
        y_c <= delta_c,
        T_c + (T_cw - T_c) * (1 - (y_c / delta_c)**(1/7)),
        T_c
    )

    # Combine
    y_all = np.concatenate([y_g, y_w, y_c + wall_thickness])
    T_all = np.concatenate([T_g, Ts_wall, T_cBL])

    # Extents for freestreams (at least 3x wall thickness, ensure covers BL)
    x_min = -delta_g - 2*wall_thickness
    x_max = wall_thickness + delta_c + 2*wall_thickness

    fig, ax = plt.subplots(figsize=(7, 5))

    # Shaded regions
    ax.axvspan(-delta_g, 0, color='lightgray', alpha=0.7)
    ax.axvspan(0, wall_thickness, color='darkgray', alpha=0.8)
    ax.axvspan(wall_thickness, wall_thickness + delta_c, color='lightgray', alpha=0.7)

    # Horizontal freestream temperature lines
    ax.hlines(T_inf, x_min, -delta_g, colors='tab:red', linestyle='--', linewidth=2)
    ax.hlines(T_c, wall_thickness + delta_c, x_max, colors='tab:red', linestyle='--', linewidth=2)

    # Temperature profile in professional red
    ax.plot(y_all, T_all, color='tab:red', linestyle='--', linewidth=2, label='Temperature profile')

    # Region labels, rotated and centered vertically
    y_lim = ax.get_ylim()
    y_mid = 0.5 * (y_lim[0] + y_lim[1])
    ax.text(x_min - 0.03*x_min, y_mid, 'Freestream gas', va='center', ha='left', rotation=45)
    ax.text(-0.5*delta_g, y_mid, 'Gas BL', va='center', ha='center', rotation=45)
    ax.text(wall_thickness/2, y_mid, 'Wall', va='center', ha='center', rotation=45)
    ax.text(wall_thickness + delta_c, y_mid, 'Coolant BL',
            va='center', ha='center', rotation=45)
    ax.text(x_max, y_mid, 'Freestream coolant', va='center', ha='right', rotation=45)

    # Labels and formatting
    ax.set_xlim(x_min, x_max)
    ax.set_xlabel('Distance from hot-wall interface, y (m)')
    ax.set_ylabel('Temperature, T (K)')
    ax.set_title(f'Temperature profile at x = {x0:.3f} m, circuit {circuit_index}')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(width=1.2)
    legend = ax.legend(frameon=True, framealpha=0.9, edgecolor='black')
    legend.get_frame().set_linewidth(1.0)

    plt.tight_layout()

import matplotlib.cm as cm
import matplotlib.colors as mcolors

import os, json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import os, json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap, Normalize
import matplotlib.ticker as ticker      # put this near your other imports

def plot_theta_vs_epsilon():
    """
    θ vs ε plot with:
      • θᵉ curves  : green→blue
      • θⁿ curves  : orange→red
      • Computer-Modern math
      • Endpoint % tags + centred θₑ / θₙ (labels in black)
      • Per-percent fine-tune of tag vertical position
      • Extra x-margin so tags stay on the canvas
    """

    # ---------------- 1  font for math ----------------

    # ---------------- 2  colour ramps -----------------
    #cmap_e = mcolors.LinearSegmentedColormap.from_list(
    #    "green_blue", ["#00c46a", "#005dff"])
    #cmap_n = mcolors.LinearSegmentedColormap.from_list(
    #    "orange_red", ["#ff9b00", "#ff1a1a"])
    #norm = mcolors.Normalize(vmin=60, vmax=100)
    def muted_cmap(name, lo=0.40, hi=0.80, steps=256):
        base = cm.get_cmap(name)
        return ListedColormap(base(np.linspace(lo, hi, steps)), name=f"{name}_muted")

    cmap_e = muted_cmap('BuGn',   0.60, 0.80)     # bluish-green   for θᵉ
    cmap_n = muted_cmap('YlOrRd', 0.40, 0.80)     # orange-to-red  for θⁿ
    norm   = Normalize(vmin=60, vmax=100)         # 60 % → 0 … 100 % → 1

    def colour(pct, which):           # pct = 60…100
        return cmap_e(norm(pct)) if which == 'e' else cmap_n(norm(pct))

    # ---------------- 3  per-label y-nudges ----------
    #    edit the dict below if two labels touch;
    #    units are “points” in annotate xytext
    y_tweak_e = {
        60: 0,      # already used for θₑ/θₙ label
        70: 0,
        80:  3,
        90: -0,
        100: -2,
    }

    y_tweak_n = {
        60: 0,      # already used for θₑ/θₙ label
        70: 0,
        80:  4,
        90: 0,
        100: -6,
    }

    # ---------------- 4  helpers & data --------------
    def pct_from_key(k):
        try:
            return int(k.split('_')[2])
        except Exception:
            return 0

    base = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(base, 'data', 'theta_e.json')) as fe:
        data_e = json.load(fe)
    with open(os.path.join(base, 'data', 'theta_n.json')) as fn:
        data_n = json.load(fn)

    # ---------------- 5  figure/axes -----------------
    plt.figure(figsize=(5, 3.5))
    plt.xscale('log')
    ax = plt.gca()
    ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0))              # keep decade ticks
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    # optional: hide the minor tick labels (they're too dense to print nicely)
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    plt.xlabel(r'$\varepsilon$')
    plt.ylabel(r'$\theta\,(^{\circ})$')
    #plt.title(r'$\theta$ vs $\varepsilon$ curves')
    plt.grid(True, linewidth=0.6, alpha=0.5)

    # ---------------- 6  θᵉ curves -------------------
    for k, c in data_e.items():
        eps, th = c['epsilon'], c['theta']
        pct = pct_from_key(k)
        col = colour(pct, 'e')
        line, = plt.plot(eps, th, lw=2, color=col)

        # endpoint %
        plt.annotate(f"{pct}%",
                     xy=(eps[-1], th[-1]),
                     xytext=(5, y_tweak_e.get(pct, 0)),
                     textcoords="offset points",
                     ha='left', va='center',
                     fontsize=9, color=col)

        # centred θₑ on the 60 % curve
        if pct == 60:
            mid = len(eps) // 2
            plt.annotate(r'$\theta_e$', xy=(eps[mid], th[mid]),
                         xytext=(0, 10), textcoords="offset points",
                         ha='center', va='bottom',
                         fontsize=11, color='black')

    # ---------------- 7  θⁿ curves -------------------
    for k, c in data_n.items():
        eps, th = c['epsilon'], c['theta']
        pct = pct_from_key(k)
        col = colour(pct, 'n')
        line, = plt.plot(eps, th, lw=2, color=col)

        plt.annotate(f"{pct}%",
                     xy=(eps[-1], th[-1]),
                     xytext=(5, y_tweak_n.get(pct, 0)),
                     textcoords="offset points",
                     ha='left', va='center',
                     fontsize=9, color=col)

        if pct == 60:
            mid = len(eps) // 2
            plt.annotate(r'$\theta_n$', xy=(eps[mid], th[mid]),
                         xytext=(0, 10), textcoords="offset points",
                         ha='center', va='bottom',
                         fontsize=11, color='black')

    # ---------------- 8  give labels breathing room --
    ax = plt.gca()
    xmin, xmax = ax.get_xlim()
    ax.set_xlim(xmin, xmax * 1.3)   # 5 % extra space on the right

    plt.tight_layout()


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

def plot_wall_temperature(
    *cooling_data_dicts,
    circuits,
    label_interfaces=None,
    annotate_ha=None,
    annotate_va=None,
    title=None,
    plot_hot=True,
    plot_interfaces=False,
    plot_coolant_wall=False
):
    """
    Plot wall temperatures for one or more cooling circuits, with:
      - a blue→red gradient that covers coolant wall → internal interfaces → hot wall,
      - optional per‐circuit labels (Cold Wall / Interface n / Hot Wall),
      - a circuit‐name annotation placed at the midpoint of the topmost curve,
        with a small “margin” offset away from the point.

    Args:
        *cooling_data_dicts:  
            One dict per circuit, each with keys:
              'x' : 1D array of axial positions
              'T' : 2D array of shape (len(x), n_layers), where:
                     - column 0 is coolant bulk (ignored here),
                     - column 1 is the coolant‐side wall,
                     - columns 2..(n_layers−2) are internal interfaces,
                     - column n_layers−1 is the hot‐side wall.
        circuits (list):  
            A list of “circuit‐like” objects (length = N).  Each must have:
              .name (string) → used for the annotation.
        label_interfaces (list of bool, optional):  
            Length = N.  If True for circuit _i_, then each of its plotted curves
            will get a legend entry (“Cold Wall”, “Interface n”, “Hot Wall”).  
            If False, that circuit’s curves are drawn without individual labels.
            Default: all False.
        annotate_ha (list of str, optional):  
            Length = N.  Horizontal alignment for each circuit annotation: 
            one of {"left", "center", "right"}.  Default: all "center".
        annotate_va (list of str, optional):  
            Length = N.  Vertical alignment for each circuit annotation: 
            one of {"top", "center", "bottom"}.  Default: all "bottom".
        title (str, optional):  
            Plot title.
        plot_hot (bool, optional):  
            If True, plot each circuit’s hot‐wall (last column).
        plot_interfaces (bool, optional):  
            If True, plot each circuit’s internal interfaces (columns 2..n_layers−2).
        plot_coolant_wall (bool, optional):  
            If True, plot each circuit’s coolant‐side wall (column 1).

    Returns:
        fig, ax
    """
    N = len(cooling_data_dicts)
    if len(circuits) != N:
        raise ValueError("`circuits` must be the same length as `cooling_data_dicts`.")

    # Default: don’t label any interfaces if not provided
    if label_interfaces is None:
        label_interfaces = [False] * N
    if len(label_interfaces) != N:
        raise ValueError("`label_interfaces` must have length = number of circuits.")

    # Default alignments
    if annotate_ha is None:
        annotate_ha = ["center"] * N
    if annotate_va is None:
        annotate_va = ["bottom"] * N
    if len(annotate_ha) != N or len(annotate_va) != N:
        raise ValueError("`annotate_ha` and `annotate_va` must each have length = number of circuits.")

    # Build a blue→red colormap once
    cmap = LinearSegmentedColormap.from_list('blue_red', ['tab:blue', 'tab:red'])

    fig, ax = plt.subplots(figsize=(8, 5))

    # Will store (x_mid, y_top) for each circuit, or None if nothing plotted
    annotation_positions = []

    for i, data in enumerate(cooling_data_dicts):
        x = np.asarray(data['x'])
        T = np.asarray(data['T'])
        circuit_name = circuits[i].name
        do_labels = bool(label_interfaces[i])

        # 1) Choose which columns to plot for this circuit
        cols_to_plot = []
        if plot_coolant_wall and T.shape[1] > 1:
            cols_to_plot.append(1)
        if plot_interfaces and T.shape[1] > 2:
            cols_to_plot.extend(range(2, T.shape[1] - 1))
        if plot_hot:
            cols_to_plot.append(T.shape[1] - 1)

        if len(cols_to_plot) == 0:
            annotation_positions.append(None)
            continue

        # 2) Pick evenly spaced colors on the blue→red gradient
        if len(cols_to_plot) == 1:
            colors = [cmap(0.5)]
        else:
            colors = [cmap(j / (len(cols_to_plot) - 1)) for j in range(len(cols_to_plot))]

        # 3) Plot and collect “midpoint” temperatures
        x_mid = 0.5 * (x.min() + x.max())
        mid_idx = np.abs(x - x_mid).argmin()
        y_values_at_mid = []

        for idx_color, col in enumerate(cols_to_plot):
            y = T[:, col]
            y_values_at_mid.append(y[mid_idx])

            if do_labels:
                if col == 1 and plot_coolant_wall:
                    line_label = "Cold Wall"
                elif col == T.shape[1] - 1 and plot_hot:
                    line_label = "Hot Wall"
                else:
                    interface_number = col - 1
                    line_label = f"Interface {interface_number}"
            else:
                line_label = None

            linestyle = '--' if (plot_interfaces and 2 <= col <= T.shape[1] - 2) else '-'

            ax.plot(
                x, y,
                color=colors[idx_color],
                linestyle=linestyle,
                label=line_label
            )

        # Find the topmost curve at x_mid
        y_top = max(y_values_at_mid)
        annotation_positions.append((x_mid, y_top))

    # 4) Draw a legend if any circuit requested per-layer labels
    if any(label_interfaces) or (plot_hot or plot_interfaces or plot_coolant_wall):
        ax.legend()

    # 5) Add each circuit’s annotation with a small offset (“margin”) away from the point
    margin = 10  # in points
    for i, pos in enumerate(annotation_positions):
        if pos is None:
            continue
        x_annot, y_annot = pos
        ha = annotate_ha[i]
        va = annotate_va[i]

        # Determine offset in points based on alignment
        if ha == 'left':
            dx = -margin
        elif ha == 'right':
            dx = margin
        else:  # 'center'
            dx = 0

        if va == 'bottom':
            dy = margin
        elif va == 'top':
            dy = -margin
        else:  # 'center'
            dy = 0

        ax.annotate(
            circuits[i].name,
            xy=(x_annot, y_annot),
            xytext=(dx, dy),
            textcoords='offset points',
            ha=ha,
            va=va,
            fontsize=10,
            fontweight='bold'
        )

    ax.set_xlabel('Axial Position (m)')
    ax.set_ylabel('Temperature (K)')
    if title:
        ax.set_title(title)

    ax.grid(True)
    plt.tight_layout()
    return fig, ax




def plot_coolant_temperature(*cooling_data_dicts, labels=None, title=None):
    """
    Plot coolant static temperature for one or more cooling simulations.

    Parameters
    ----------
    *cooling_data_dicts : dict
        Must contain keys 'x' and 'T_static'.
    labels : sequence of str, optional
    title : str, optional
    """
    fig, ax = plt.subplots(figsize=(8, 5))
    for i, data in enumerate(cooling_data_dicts):
        x = np.asarray(data['x'])
        y = np.asarray(data['T_static'])
        lbl = labels[i] if labels is not None else None
        ax.plot(x, y, color='tab:red', label=lbl)

    ax.set_xlabel('Axial Position (m)')
    ax.set_ylabel('Coolant Temperature (K)')
    if title:
        ax.set_title(title)
    if labels:
        ax.legend()
    ax.grid(True)
    plt.tight_layout()
    return fig, ax


def plot_coolant_pressure(*cooling_data_dicts, static=True, stagnation=True, labels=None, title=None):
    """
    Plot coolant static pressure for one or more cooling simulations.

    Parameters
    ----------
    *cooling_data_dicts : dict
        Must contain keys 'x' and 'p_static'.
    labels : sequence of str, optional
    title : str, optional
    """
    fig, ax = plt.subplots(figsize=(8, 5))
    for i, data in enumerate(cooling_data_dicts):
        x = np.asarray(data['x'])
        y1 = np.asarray(data['p_static'])
        y2 = np.asarray(data['p_stagnation'])
        lbl = labels[i] if labels is not None else None
        if static:
            ax.plot(x, y1, color='tab:red', label="Static Pressure")
        if stagnation:
            ax.plot(x, y2, color='tab:purple', label="Stagnation Pressure")

    ax.set_xlabel('Axial Position (m)')
    ax.set_ylabel('Pressure (Pa)')
    if title:
        ax.set_title(title)
    if labels:
        ax.legend()
    ax.grid(True)
    plt.tight_layout()
    return fig, ax


def plot_heat_flux(*cooling_data_dicts, labels=None, title=None):
    """
    Plot heat flux (dQ/dA) for one or more cooling simulations.

    Parameters
    ----------
    *cooling_data_dicts : dict
        Must contain keys 'x' and 'dQ_dA'.
    labels : sequence of str, optional
    title : str, optional
    """
    fig, ax = plt.subplots(figsize=(8, 5))
    for i, data in enumerate(cooling_data_dicts):
        x = np.asarray(data['x'])
        y = np.asarray(data['dQ_dA'])
        lbl = labels[i] if labels is not None else None
        ax.plot(x, y, color='tab:red', label=lbl)

    ax.set_xlabel('Axial Position (m)')
    ax.set_ylabel('Heat Flux (W/m²)')
    if title:
        ax.set_title(title)
    if labels:
        ax.legend()
    ax.grid(True)
    plt.tight_layout()
    return fig, ax


def plot_velocity(*cooling_data_dicts, labels=None, title=None):
    """
    Plot coolant velocity for one or more cooling simulations.

    Parameters
    ----------
    *cooling_data_dicts : dict
        Must contain keys 'x' and 'velocity'.
    labels : sequence of str, optional
    title : str, optional
    """
    fig, ax = plt.subplots(figsize=(8, 5))
    for i, data in enumerate(cooling_data_dicts):
        x = np.asarray(data['x'])
        y = np.asarray(data['velocity'])
        lbl = labels[i] if labels is not None else None
        ax.plot(x, y, color='tab:red', label=lbl)

    ax.set_xlabel('Axial Position (m)')
    ax.set_ylabel('Velocity (m/s)')
    if title:
        ax.set_title(title)
    if labels:
        ax.legend()
    ax.grid(True)
    plt.tight_layout()
    return fig, ax
