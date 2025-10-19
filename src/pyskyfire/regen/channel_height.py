from typing import Sequence, Callable, Union
import numpy as np
from pyskyfire.regen.thrust_chamber import Contour

def make_channel_height_fn(
    contour: Contour,
    region_fractions: Sequence[float],
    flat_heights: Sequence[float],
    pinch_factors: Sequence[float],
    transition_widths: Union[float, Sequence[float]],
    logistic_k: float = 10.0
) -> Callable[[float], float]:
    """Construct a smooth channel-height profile along a thrust chamber contour.

    Builds a continuous, piecewise-logistic blend of per-region channel heights
    across the chamber and nozzle. Each region is defined by a normalized axial
    coordinate fraction relative to the throat, and can optionally "pinch"
    channel height with radius.

    The returned function evaluates the local channel height at any axial
    coordinate ``x`` along the contour.

    Parameters
    ----------
    contour : Contour
        Thrust-chamber contour defining the wall coordinates ``x`` and ``r(x)``.
    region_fractions : Sequence[float]
        Normalized axial coordinates marking region boundaries:

        - ``-1.0`` → chamber inlet
        - ``0.0`` → throat
        - ``+1.0`` → nozzle exit
    flat_heights : Sequence[float]
        Base (unpinched) channel heights for each region [m].
    pinch_factors : Sequence[float]
        Fraction in ``[0, 1]`` describing how strongly each region’s height scales
        with local radius. ``0`` means constant height, ``1`` means fully proportional
        to radius.
    transition_widths : float or Sequence[float]
        Axial transition width(s) controlling how quickly heights blend between
        adjacent regions. If a single scalar is given, it is used for all boundaries.
    logistic_k : float, optional
        Steepness of the logistic blending function. Higher values yield sharper
        transitions. Default is ``10.0``.

    Returns
    -------
    Callable[[float], float]
        A callable function ``channel_height(x)`` that returns the local channel
        height [m] at axial coordinate ``x``.

    Raises
    ------
    ValueError
        If the number of transition widths does not match ``len(flat_heights) - 1``.

    Notes
    -----
    The normalized coordinate mapping is:

    - Negative region fractions → chamber side
    - Positive region fractions → nozzle side

    Transitions are blended smoothly using a logistic weighting function centered
    at each region boundary. Channel heights can pinch with radius according to
    ``pinch_factors``.

    Examples
    --------
    >>> channel_height_fn = psf.regen.make_channel_height_fn(
                contour=contour, 
                region_fractions=[-1.0, 0.25, 1.0], 
                flat_heights= [0.0032, 0.00134], 
                pinch_factors= [0.6, -5.0], 
                transition_widths=[0.1]
    >>> fn(0.12)  # Evaluate channel height at some axial x
    0.0034
    """
    # 1) find throat and endpoints
    x_start = contour.xs[0]
    x_end   = contour.xs[-1]
    # throat = x at minimum radius
    throat_idx = int(np.argmin(contour.rs))
    x_throat   = contour.xs[throat_idx]

    # 2) map each fraction to its absolute x-location
    bounds_x = []
    L_chamber = x_throat - x_start
    L_nozzle  = x_end   - x_throat
    for f in region_fractions:
        if f < 0:
            # negative: fraction of chamber length
            bounds_x.append(x_throat + f * L_chamber)
        else:
            # zero or positive: fraction of nozzle length
            bounds_x.append(x_throat + f * L_nozzle)

    # 3) radius scaling (unchanged)
    rmin, rmax = min(contour.rs), max(contour.rs)
    dr = (rmax - rmin) or 1e-12

    # 4) build per-zone height funcs (unchanged)
    H_funcs = []
    for h_flat, pinch in zip(flat_heights, pinch_factors):
        def H(x, h_flat=h_flat, pinch=pinch):
            r  = contour.r(x)
            sr = (r - rmin) / dr
            return (1 - pinch)*h_flat + pinch*(sr * h_flat)
        H_funcs.append(H)

    # 5) normalize transition_widths (unchanged)
    n_bounds = len(H_funcs) - 1
    if isinstance(transition_widths, (float, int)):
        tw_list = [float(transition_widths)] * n_bounds
    else:
        if len(transition_widths) != n_bounds:
            raise ValueError(f"Expected {n_bounds} widths, got {len(transition_widths)}")
        tw_list = list(transition_widths)

    # 6) final blended channel_height (unchanged)
    def channel_height(x: float) -> float:
        h_val = H_funcs[0](x)
        for xi, H_next, tw in zip(bounds_x[1:], H_funcs[1:], tw_list):
            k = logistic_k / tw
            w = 1.0 / (1.0 + np.exp(-k * (x - xi)))
            h_val = (1 - w)*h_val + w*H_next(x)
        return h_val

    return channel_height
