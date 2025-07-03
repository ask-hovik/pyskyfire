import matplotlib.pyplot as plt

def plot_transport_properties(*cts,
                              properties=None,
                              labels=None,
                              title=None):
    """
    Plot specified transport-property maps from one or more CombustionTransport objects.

    Parameters
    ----------
    *cts : CombustionTransport
        One or more combustion‐transport objects with precomputed .x_domain and 
        maps: M_map, gamma_map, T_map, p_map, h_map, cp_map, k_map, mu_map, Pr_map, rho_map.
    properties : list of str, optional
        Which properties to plot. Possible keys:
            'M', 'gamma', 'T', 'p', 'h', 'cp', 'k', 'mu', 'Pr', 'rho'
        Default is all of the above in that order.
    labels : list of str, optional
        Legend labels corresponding to each cts. If omitted, no legend is drawn.
    title : str, optional
        Figure-level title (suptitle).
    """
    plt.style.use('ggplot')
    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Latin Modern Roman'],
        'font.size': 10, 
        'savefig.dpi': 600
    })
    # default ordering
    all_props = ['M','gamma','T','p','h','cp', 'cv', 'k','mu','Pr','rho']
    prop_info = {
        'M':     ('M', ''),
        'gamma': ('gamma', ''),
        'T':     ('T', 'K'),
        'p':     ('p', 'Pa'),
        'h':     ('h', 'J/kg'),
        'cp':    ('cp (mass)', 'J/(kg·K)'),
        'cv':    ('cv (mass)', 'J/(kg·K)'),
        'k':     ('k', 'W/(m·K)'),
        'mu':    ('mu', 'Pa·s'),
        'Pr':    ('Pr', ''),
        'rho':   ('rho', 'kg/m³'),
    }

    if properties is None:
        properties = all_props
    else:
        # validate
        for prop in properties:
            if prop not in all_props:
                raise ValueError(f"Unknown property '{prop}'. Valid keys: {all_props}")

    nplots = len(properties)
    fig, axs = plt.subplots(nplots, 1,
                             figsize=(8, 2.5*nplots),
                             sharex=True)
    # Ensure axs is iterable
    if nplots == 1:
        axs = [axs]

    for ax, prop in zip(axs, properties):
        key = f"{prop}_map"
        label_text, unit = prop_info[prop]
        for i, ct in enumerate(cts):
            x = ct.x_domain
            y = getattr(ct, key)
            lbl = labels[i] if labels and i < len(labels) else None
            ax.plot(x, y, label=lbl)
        ax.set_ylabel(f"{label_text}" + (f" ({unit})" if unit else ""))
        ax.grid(True)
        if labels:
            ax.legend(frameon=True, framealpha=0.9)

    axs[-1].set_xlabel("Axial position, x (m)")
    if title:
        fig.suptitle(title, y=1.02, fontsize=14)
    plt.tight_layout()
    return fig, axs
