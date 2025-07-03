# engine_network.py
# ---------------------------------------------------------------------
import numpy as np
from dataclasses import dataclass
from typing      import Dict, List, Tuple, Iterable, Any
import pyskyfire as psf
import matplotlib.pyplot as plt

# ────────────────────────────────────────────────────────────────────
#  Tiny helper –– you already have this in blocks.py, re-import here
# ────────────────────────────────────────────────────────────────────
@dataclass
class Station:
    p:    float
    T:    float
    mdot: float = float("nan")

# ────────────────────────────────────────────────────────────────────
class EngineNetwork:
    """
    Fixed-point network that stores *two* dictionaries:
        • stations : thermodynamic states        (Station objects)
        • signals  : scalar data (Δp, powers, targets, flags …)
    Every Block advertises four lists:
        station_inputs,  station_outputs,
        signal_inputs,   signal_outputs
    and its compute() returns (stations_out, signals_out).
    """

    # --------------------------------------------------------------
    def __init__(self,
                 stations: Dict[str, Station],
                 signals:  Dict[str, float],
                 blocks:   List["Block"]):

        self.stations = dict(stations)   # shallow copy so caller keeps their own
        self.signals  = dict(signals)
        self.blocks   = self._toposort(blocks)  # ← optionally reorder
        self.residuals: List[float] = []

    # --------------------------------------------------------------
    def _toposort(self, blocks: Iterable["Block"]) -> List["Block"]:
        """
        Very small dependency sorter.  If you’re sure you already passed
        the blocks in a safe order, just return list(blocks).
        """
        # SIMPLE IMPLEMENTATION – keeps given order:
        return list(blocks)

        # For a real sort, build a dependency graph using the four metadata
        # lists and run Kahn’s algorithm here.

    # --------------------------------------------------------------
    def _merge_and_residual(self,
                            store: Dict[str, Any],
                            delta: Dict[str, Any],
                            res:   float) -> float:
        """
        Merge *delta* into *store* and update the running maximum
        relative residual *res*.
        """
        for k, v_new in delta.items():
            v_old = store.get(k)
            if v_old is not None:
                # ---- Station residuals
                if isinstance(v_new, Station) and isinstance(v_old, Station):
                    res = max(
                        res,
                        abs(v_new.p    - v_old.p   ) / max(abs(v_old.p   ), 1e-10),
                        abs(v_new.T    - v_old.T   ) / max(abs(v_old.T   ), 1e-10),
                        0.0 if (np.isnan(v_old.mdot) or np.isnan(v_new.mdot))
                        else abs(v_new.mdot - v_old.mdot) / max(abs(v_old.mdot), 1e-10)
                    )
                # ---- Scalar residuals
                elif isinstance(v_new, (int, float)) and isinstance(v_old, (int, float)):
                    res = max(res, abs(v_new - v_old) / max(abs(v_old), 1e-8))
            # write value (new or overwrite)
            store[k] = v_new
        return res

    # --------------------------------------------------------------
    def update(self) -> float:
        """
        Executes one full sweep over all blocks.
        Returns the maximum relative change (residual) seen.
        """
        residual = 0.0
        for blk in self.blocks:
            ds, sg = blk.compute(self.stations, self.signals)
            residual = self._merge_and_residual(self.stations, ds, residual)
            residual = self._merge_and_residual(self.signals,  sg, residual)
        return residual

    # --------------------------------------------------------------
    def run_fixed_point(self, tol: float = 1e-6, max_iter: int = 100):
        """
        Iterate until the largest relative change among *all* advertised
        outputs is below *tol* or until *max_iter* is reached.
        """
        for i in range(1, max_iter + 1):
            
            res = self.update()

            print(f"[fixed-point] iter {i:3d} → residual = {res:.3e}")
            self.residuals.append(res)
            if res < tol:
                print(f"Converged in {i} iterations (residual = {res:.3e})")
                break
        else:
            raise RuntimeError(f"Did not converge after {max_iter} iterations; "
                           f"last residual = {res:.3e}")

        # -----------------------------
        #  POST-PROCESS SWEEP
        # -----------------------------
        self.block_results: dict[str, dict[str, Any]] = {}

        for blk in self.blocks:
            out = blk.post_process(self.stations, self.signals)   # always exists
            if out:                                               # skip empty dicts
                self.block_results[blk.name] = out