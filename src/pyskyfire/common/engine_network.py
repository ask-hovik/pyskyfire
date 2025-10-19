# engine_network.py
# ---------------------------------------------------------------------
import numpy as np
from dataclasses import dataclass
from typing      import Dict, List, Iterable, Any


@dataclass
class Station:
    """Thermodynamic state container used by the engine network.

    Attributes
    ----------
    p : float
        Static/stagnation pressure [Pa].
    T : float
        Static/stagnation temperature [K].
    mdot : float
        Mass flow rate [kg/s]. Default is ``NaN`` if unknown.
    """

    p:    float
    T:    float                     
    mdot: float = float("nan")      

# ────────────────────────────────────────────────────────────────────
class EngineNetwork:
    """Fixed-point engine network orchestrating stations and signals.

    The network manages two dictionaries:

    - **stations**: thermodynamic states as :class:`Station` objects.
    - **signals**: scalar data (e.g., pressure drops, powers, targets, flags).

    Blocks added to the network advertise four metadata lists
    (``station_inputs``, ``station_outputs``, ``signal_inputs``, ``signal_outputs``)
    and implement ``compute(stations, signals) -> (stations_out, signals_out)``.
    The network executes blocks in a chosen order and applies a fixed-point
    iteration until convergence.

    Parameters
    ----------
    stations : dict[str, Station]
        Initial station map.
    signals : dict[str, float]
        Initial scalar signal map.
    blocks : list[Block]
        Block instances participating in the network.

    Attributes
    ----------
    stations : dict[str, Station]
        Current station map keyed by name.
    signals : dict[str, float]
        Current scalar signal map keyed by name.
    blocks : list[Block]
        Processing blocks in execution order.
    residuals : list[float]
        History of per-iteration maximum relative residuals.
    block_results : dict[str, dict[str, Any]]
        Optional post-process results keyed by block name, populated after
        convergence.

    See Also
    --------
    FluidBlock
        Base class for blocks that transform fluid stations.
    SignalBlock
        Base class for blocks that operate on scalar signals only.

    Notes
    -----
    The default execution order preserves the given block list. A dependency-
    aware topological sort can replace it; see :meth:`_toposort`.
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
        """Compute an execution order for blocks.

        The default implementation preserves the given order. Replace with a
        dependency-aware topological sort (e.g., Kahn’s algorithm) using each
        block’s ``station_*`` and ``signal_*`` metadata for strict ordering.

        Parameters
        ----------
        blocks : Iterable[Block]
            Blocks to order.

        Returns
        -------
        list[Block]
            Execution order for the fixed-point sweep.

        Notes
        -----
        Sorting is currently not implemented; the user is responsible for
        providing a solvable order.
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
        """Merge updates into a store and track the max relative residual.

        For each key in ``delta``:

        - If values are :class:`Station`, compute componentwise relative changes
          for ``p``, ``T`` and (if both present) ``mdot``; update ``res`` with
          the maximum.
        - If values are scalars, compute a relative change and update ``res``.
        - Overwrite/write the value in ``store``.

        Parameters
        ----------
        store : dict[str, Any]
            Destination dictionary to update in-place.
        delta : dict[str, Any]
            New values produced by a block this pass.
        res : float
            Current running maximum relative residual.

        Returns
        -------
        float
            Updated maximum relative residual.

        Notes
        -----
        To avoid division by zero, denominators are clamped with small epsilons:
        ``1e-10`` for stations and ``1e-8`` for scalars. If either ``mdot`` is
        ``NaN``, the mass-flow term is skipped for that key.
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
        """Execute one full sweep over all blocks.

        Calls each block’s :meth:`compute` with the current ``stations`` /
        ``signals``, merges its outputs, and tracks the maximum relative change
        across all updated entries.

        Returns
        -------
        float
            Maximum relative residual observed this sweep.
        """
        residual = 0.0
        for blk in self.blocks:
            ds, sg = blk.compute(self.stations, self.signals)
            residual = self._merge_and_residual(self.stations, ds, residual)
            residual = self._merge_and_residual(self.signals,  sg, residual)
        return residual

    # --------------------------------------------------------------
    def run_fixed_point(self, tol: float = 1e-6, max_iter: int = 100):
        """Iterate fixed-point sweeps until convergence or iteration cap.

        Repeatedly calls :meth:`update` until the maximum relative residual falls
        below ``tol`` or until ``max_iter`` sweeps have been performed.

        Parameters
        ----------
        tol : float, optional
            Convergence threshold on max relative residual. Default is ``1e-6``.
        max_iter : int, optional
            Maximum number of sweeps. Default is ``100``.

        Raises
        ------
        RuntimeError
            If the network does not converge within ``max_iter``.

        Notes
        -----
        Appends per-iteration residuals to :attr:`residuals`. On convergence,
        performs a post-process sweep calling each block’s :meth:`post_process`
        and stores any non-empty results in :attr:`block_results`.
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