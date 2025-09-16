'''import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Patch

def set_axes_equal(ax):
    """
    Sets equal scaling for a 3D plot by adjusting the axis limits.
    This ensures that a sphere appears as a sphere, not an ellipsoid.
    """
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    y_range = abs(y_limits[1] - y_limits[0])
    z_range = abs(z_limits[1] - z_limits[0])
    max_range = max([x_range, y_range, z_range])
    
    x_mid = np.mean(x_limits)
    y_mid = np.mean(y_limits)
    z_mid = np.mean(z_limits)
    
    ax.set_xlim3d([x_mid - max_range / 2, x_mid + max_range / 2])
    ax.set_ylim3d([y_mid - max_range / 2, y_mid + max_range / 2])
    ax.set_zlim3d([z_mid - max_range / 2, z_mid + max_range / 2])

def plot_impeller_3D(imp):
    """
    Plots the impeller geometry in 3D.
    
    This function does the following:
      1. Plots each meridional blade curve (streamline) in Cartesian coordinates.
      2. Creates semi-transparent surfaces for the lower and upper shrouds by revolving
         the first and last streamlines about the z-axis.
      3. Constructs a semi-transparent surface for the impeller blades by revolving the entire
         set of streamlines (which form the blade meridional geometry) for each blade. The blades
         are distributed evenly around the impeller (using the parameter z for the number of blades).
    
    Parameters:
      new_meridionals (list): A list of streamlines, each a list of (z, r, theta) tuples.
                              These curves represent the blade geometry in the meridional plane.
      z (int): The number of blades.
      title (str): Title for the 3D plot.
      show (bool): If True, displays the plot.
    """

    new_meridionals = imp.geometry["streamlines_3D"]
    z = imp.z
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # ---------------------------------------
    # 1. Plot the individual streamlines.
    # ---------------------------------------
    """
    for streamline in new_meridionals:
        # Unpack (z, r, theta) for this streamline.
        z_coords, r_coords, theta_coords = zip(*streamline)
        r_arr = np.array(r_coords)
        # Convert theta (in degrees) to radians.
        theta_arr = np.deg2rad(np.array(theta_coords))
        # Convert to Cartesian coordinates.
        x_arr = r_arr * np.cos(theta_arr)
        y_arr = r_arr * np.sin(theta_arr)
        ax.plot(x_arr, y_arr, z_coords, color='k', lw=1)
    """
    # ---------------------------------------------------------
    # 2. Plot shroud surfaces by revolving the first and last streamlines.
    # ---------------------------------------------------------
    lower_shroud = new_meridionals[0]
    upper_shroud = new_meridionals[-1]
    
    lower_shroud = np.array(lower_shroud)  # shape (N, 3)
    upper_shroud = np.array(upper_shroud)
    
    z_lower = lower_shroud[:, 0]
    r_lower = lower_shroud[:, 1]
    z_upper = upper_shroud[:, 0]
    r_upper = upper_shroud[:, 1]
    
    # Define a full revolution.
    theta_full = np.linspace(0, 2*np.pi, 100)
    
    # Lower shroud:
    Theta_lower, Z_lower_mesh = np.meshgrid(theta_full, z_lower)
    R_lower_mesh = np.tile(r_lower, (len(theta_full), 1)).T  # replicate along columns
    X_lower = R_lower_mesh * np.cos(Theta_lower)
    Y_lower = R_lower_mesh * np.sin(Theta_lower)
    ax.plot_surface(X_lower, Y_lower, Z_lower_mesh, color='blue', alpha=0.3,
                    rstride=4, cstride=4, linewidth=0)
    
    # Upper shroud:
    Theta_upper, Z_upper_mesh = np.meshgrid(theta_full, z_upper)
    R_upper_mesh = np.tile(r_upper, (len(theta_full), 1)).T
    X_upper = R_upper_mesh * np.cos(Theta_upper)
    Y_upper = R_upper_mesh * np.sin(Theta_upper)
    ax.plot_surface(X_upper, Y_upper, Z_upper_mesh, color='red', alpha=0.3,
                    rstride=4, cstride=4, linewidth=0)
    
    # ---------------------------------------------------------
    # 3. Plot the blade surfaces.
    # ---------------------------------------------------------
    # Combine the meridional streamlines into a 2D grid.
    # Assume new_meridionals is organized such that each streamline represents a row
    # and points along each streamline represent the columns.
    blade_surface = np.array(new_meridionals)  # shape: (n_streamlines, n_points, 3)
    
    # Extract the meridional coordinates.
    Z_blade = blade_surface[:, :, 0]       # axial coordinate (n_streamlines x n_points)
    R_blade = blade_surface[:, :, 1]       # radial coordinate
    Theta_blade = blade_surface[:, :, 2]   # twist angle (in degrees)
    
    # For each blade, rotate the entire blade surface by an offset so the blades are evenly spaced.
    for i in range(z):
        offset = i * (360.0 / z)
        # Add the offset to the twist angle.
        TotalTheta = Theta_blade + offset
        TotalTheta_rad = np.deg2rad(TotalTheta)
        X_blade = R_blade * np.cos(TotalTheta_rad)
        Y_blade = R_blade * np.sin(TotalTheta_rad)
        ax.plot_surface(X_blade, Y_blade, Z_blade, color='green', alpha=0.3,
                        rstride=1, cstride=1, linewidth=0)
    
    # Set labels and title.
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    
    # Set equal scaling for all axes.
    set_axes_equal(ax)
    plt.show()

def plot_impeller_views(imp):
    """
    Plots two different views of the impeller (given an Impeller object):
      1) Meridional (z-r plane)
      2) Plan (top view, x-y plane)

    The function uses imp.geometry["streamlines_3D"], which should be
    a list of streamlines. Each streamline is a list of (z, r, theta) tuples.

    The plan view also:
      - Draws circles for imp.d_2, imp.d_1, and imp.d_n.
      - Replicates the streamlines around the origin to match the number
        of blades (imp.z).
    """

    streamlines_3D = imp.geometry["streamlines_3D"]
    num_blades = imp.z

    # -------------------------
    # 1) Meridional View (Side)
    # -------------------------
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)

    # Plot each streamline in z-r
    for streamline in streamlines_3D:
        z_coords, r_coords, _ = zip(*streamline)  # we ignore theta here
        ax1.plot(z_coords, r_coords, color='b', lw=1)

    ax1.set_xlabel("Axial Coordinate (z)")
    ax1.set_ylabel("Radial Coordinate (r)")
    ax1.set_title("Meridional View (Side Projection)")
    ax1.grid(True)
    ax1.set_aspect('equal')

    # Show the first figure without blocking
    plt.show(block=False)
    plt.pause(0.1)  # ensures the figure is displayed properly

    # -------------------------
    # 2) Plan View (Top)
    # -------------------------
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)

    # Plot the streamlines in x-y plane, rotating them for each blade
    for i in range(num_blades):
        offset_deg = i * (360.0 / num_blades)  # even spacing
        for streamline in streamlines_3D:
            _, r_coords, theta_coords = zip(*streamline)
            # Convert to radians with offset
            theta_offset = np.deg2rad(np.array(theta_coords) + offset_deg)
            x_coords = np.array(r_coords) * np.cos(theta_offset)
            y_coords = np.array(r_coords) * np.sin(theta_offset)
            ax2.plot(x_coords, y_coords, color='g', lw=1)

    # Draw circles for outer diameter (d_2), suction eye (d_1), and inner diameter (d_n)
    # Using dashed lines for clarity
    diameters = {
        "Inner (d_n)": imp.d_n,
        "Suction Eye (d_1)": imp.d_1,
        "Outer (d_2)": imp.d_2,
    }
    colors = ["r", "orange", "blue"]
    legend_patches = []  # List to store legend handles

    for (label, diameter), color in zip(diameters.items(), colors):
        circle = Circle((0, 0), radius=diameter/2, fill=False, color=color, lw=1.5, alpha=0.8, linestyle='dashed')
        ax2.add_patch(circle)
        legend_patches.append(Patch(color=color, linestyle='dashed', label=label))  # Add to legend

    # Add the legend in a separate box
    ax2.legend(handles=legend_patches, loc="upper right", fontsize=9, frameon=True)

    ax2.set_xlabel("X Coordinate")
    ax2.set_ylabel("Y Coordinate")
    ax2.set_title("Plan View (Top Projection)")
    ax2.grid(True)
    ax2.set_aspect('equal')

    # Show the second figure (blocking=False if you want code to continue afterward)
    plt.show()
    #plt.pause(0.1)

    print("Both plots displayed, continuing execution...")
'''