import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Latin Modern Roman'],
    'font.size': 10,
    'mathtext.fontset': 'cm',                 # ← built-in Computer Modern
    'savefig.dpi': 600,
})

# ----------------------------------------------------------------------
# Core correlations
# ----------------------------------------------------------------------
RE_LAM_MAX = 2300
RE_TURB_MIN = 3500

def _f_laminar(Re):
    """Laminar (Hagen–Poiseuille) branch."""
    return 64.0 / Re

def _colebrook_white(Re, rel_rough):
    """Colebrook–White equation solved by Newton iteration."""
    def eqn(f):
        return 1.0/np.sqrt(f) + 2.0*np.log10(rel_rough/3.71 + 2.51/(Re*np.sqrt(f)))
    return fsolve(eqn, x0=0.02, xtol=1e-12, maxfev=200)[0]

def _f_turbulent(Re, rel_rough):
    """Turbulent branch: Colebrook–White (rough) or Putukhov/Blasius (smooth)."""
    if rel_rough == 0:
        return (0.79*np.log(Re) - 1.64)**-2        # smooth correlation
    return _colebrook_white(Re, rel_rough)

def darcy_friction(Re, rel_rough):
    """Piece-wise f(Re) with linear interpolation in the transitional band."""
    if Re < RE_LAM_MAX:
        return _f_laminar(Re)
    elif Re < RE_TURB_MIN:
        f_lam = _f_laminar(Re)
        f_tur = _f_turbulent(Re, rel_rough)
        return np.interp(Re, [RE_LAM_MAX, RE_TURB_MIN], [f_lam, f_tur])
    return _f_turbulent(Re, rel_rough)

# ----------------------------------------------------------------------
# Plotting utility
# ----------------------------------------------------------------------
def plot_darcy_vs_re(
    rel_rough_list=None,
    Re_min=7e2,           # <-- starts at 500
    Re_max=1e8,
    n_pts=400
):
    """
    Moody-style plot of Darcy friction factor vs Re for a list of ε/D values.

    Parameters
    ----------
    rel_rough_list : list[float], optional
        Relative roughness values ε/D to plot.  Defaults to
        [1e-6, 5e-6, ..., 5e-2].
    Re_min, Re_max : float
        Range of Reynolds numbers.
    n_pts : int
        Number of points per curve (log-spaced).
    """
    # -------- plotting style --------

    if rel_rough_list is None:
        rel_rough_list = [ 
                  3e-5,
            1e-4, 3e-4,
            1e-3, 3e-3,
            1e-2, 3e-2,
        ]

    Re_vals = np.logspace(np.log10(Re_min), np.log10(Re_max), n_pts)

    fig, ax = plt.subplots(figsize=(8, 5))

    for eps_D in rel_rough_list:
        f_vals = [darcy_friction(Re, eps_D) for Re in Re_vals]
        ax.loglog(Re_vals, f_vals, color='tab:blue', linewidth=1.1)

        # Position the label slightly up-and-left of the last point
        label = f"{eps_D:.6f}".rstrip('0').rstrip('.')
        x_txt = Re_vals[-1] * 0.88
        y_txt = f_vals[-1] * 1.05
        ax.text(
            x_txt, y_txt,
            label,
            ha='right', va='bottom',
            fontsize=8
        )
    from matplotlib.ticker import FixedLocator, FormatStrFormatter, NullLocator, NullFormatter

    # 1) Get current y‐limits
    ymin, ymax = ax.get_ylim()

    # 2) Compute the first and last 0.01‐aligned values
    start = np.ceil(ymin * 100) / 100
    stop  = np.floor(ymax * 100) / 100

    # 3) Create ticks at every 0.01
    yticks = np.arange(start, stop + 1e-8, 0.01)

    # 4) Apply them as the ONLY major ticks
    ax.yaxis.set_major_locator(FixedLocator(yticks))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    # 5) Suppress all minor ticks
    ax.yaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_formatter(NullFormatter())


    ax.set_xlabel(r"Reynolds number, $Re$")
    ax.set_ylabel(r"Friction factor, $f$")
    ax.grid(which='both')
    ax.grid(True, which='major', axis='y')

    ax.set_xlim([Re_min, Re_max])
    plt.tight_layout()
    return fig, ax


# ----------------------------------------------------------------------
# Example usage ---------------------------------------------------------
if __name__ == "__main__":
    plot_darcy_vs_re()
    plt.show()
