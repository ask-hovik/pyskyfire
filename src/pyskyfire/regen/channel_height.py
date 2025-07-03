from typing import Sequence, Callable, Union, Tuple
import numpy as np
from pyskyfire.regen import Contour

def make_channel_height_fn(
    contour: Contour,
    region_fractions: Sequence[float],
    flat_heights: Sequence[float],
    pinch_factors: Sequence[float],
    transition_widths: Union[float, Sequence[float]],
    logistic_k: float = 10.0
) -> Callable[[float], float]:
    """
    Like before, but now:
      - f = -1.0 maps to the chamber inlet,
      - f =  0.0 maps to the throat (min radius),
      - f = +1.0 maps to the nozzle exit.
    Values in-between interpolate linearly within the chamber or nozzle.
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
