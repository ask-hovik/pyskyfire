import numpy as np

from pyskyfire.regen.design_inputs import (
    DesignContext,
    contour_providers,
    make_evaluator,
)

class ChannelPlacement:
    """Base class defining how coolant channels are wound around the contour.

    A placement owns the circumferential packing and the winding rate; where a
    section then sits radially is resolved from the cross section and the wall
    stack rather than from the placement. On its own this class describes
    plain axial channels, and subclasses (surface, internal) refine it.

    Parameters
    ----------
    n_channel_positions : int
        Number of circumferential channel locations (“leaves”).
    channel_angular_width : callable or None, optional
        Function returning angular width [rad]. May be `None` if uniform.
    occludes : bool, optional
        Whether this placement blocks part of the wall from hot gas view.

    Attributes
    ----------
    n_channel_positions : int
        Number of circumferential positions.
    channel_angular_width : callable or None
        Angular width function or constant.
    occludes : bool
        Whether the placement occludes hot surface area.

    See Also
    --------
    SurfacePlacement
    InternalPlacement
    """
    def __init__(
        self,
        n_channel_positions: int,
        channel_angular_width=None,
        occludes: bool = True,
        helix_angle=0.0,
    ):
        self.n_channel_positions = n_channel_positions   # << new home
        self.occludes = occludes

        self._channel_angular_width = channel_angular_width
        self._channel_angular_width_evaluator = (
            None
            if channel_angular_width is None
            else make_evaluator(
                channel_angular_width,
                slot="channel_angular_width",
                label=f"{type(self).__name__}.channel_angular_width",
            )
        )

        self._helix_angle = helix_angle
        self._helix_angle_evaluator = make_evaluator(
            helix_angle,
            slot="helix_angle",
            label=f"{type(self).__name__}.helix_angle",
        )

    def __setstate__(self, state):
        """Load placements serialized before the angular-width rename."""
        if (
            "channel_angular_width" not in state
            and "channel_width" in state
        ):
            state["channel_angular_width"] = state.pop("channel_width")
        state.pop("_circuit", None)
        self.__dict__.update(state)

    # ---------------- design inputs ----------------
    def _context(
        self,
        x,
        contour=None,
        circuit=None,
        *,
        include_section: bool,
    ) -> DesignContext:
        """Build the design-input context for one station.

        A placement may be shared by several circuits, so the circuit is passed
        in per call rather than stored. Without one, only tier A is available.
        """
        if circuit is not None:
            return circuit.design_context(x, include_section=include_section)
        return DesignContext(contour_providers(contour, float(x)))

    def helix_angle(self, x, contour=None, circuit=None) -> float:
        """Return the local helix angle [rad], measured from the meridian."""
        evaluator = self._helix_angle_evaluator
        if evaluator.constant_value is not None:
            return evaluator.constant_value
        # include_section reaches the channel height, which is built before any
        # helix angle is evaluated. The pitch quantities are still absent, so a
        # helix law cannot circularly depend on the pitch it sets.
        return evaluator(
            self._context(x, contour, circuit, include_section=True)
        )

    @property
    def has_channel_angular_width(self) -> bool:
        """Whether this placement defines its own angular-width law."""
        return self._channel_angular_width_evaluator is not None

    def channel_angular_width(self, x, contour=None, circuit=None) -> float:
        """Return the user-specified angular width [rad] at ``x``."""
        evaluator = self._channel_angular_width_evaluator
        if evaluator is None:
            raise AttributeError(
                f"{type(self).__name__} does not define channel_angular_width"
            )
        if evaluator.constant_value is not None:
            return evaluator.constant_value
        return evaluator(
            self._context(x, contour, circuit, include_section=True)
        )

    def channel_count(self) -> int:
        return self.n_channel_positions
    
    def dtheta_dx(self, x: float, contour, circuit=None) -> float:
        """Return dθ/dx [rad/m] for the representative lane centerline.
        Default: 0 (axial/vertical channels).
        """
        return 0.0

class SurfacePlacement(ChannelPlacement):
    """Placement model for surface-mounted cooling channels.

    Channels are positioned just outside the main wall stack,
    offset by total thickness and corrected for local contour angle.

    Parameters
    ----------
    n_channel_positions : int
        Number of channels around the circumference.

    Attributes
    ----------
    n_channel_positions : int
        Number of circumferential leaves.
    n_channels_per_leaf : int
        Fixed to 1 for surface channels.
    occludes : bool
        Always True — these channels block hot-side area.

    See Also
    --------
    ChannelPlacement
    InternalPlacement
    """
    def __init__(self, n_channel_positions: int, helix_angle=0.0):
        self.n_channels_per_leaf = 1
        super().__init__(
            n_channel_positions,
            channel_angular_width=None,
            occludes=True,
            helix_angle=helix_angle,
        )

    def dtheta_dx(self, x, contour, circuit=None):
        gamma = self.helix_angle(x, contour, circuit)
        r = float(contour.r(x))
        drdx = contour.dr_dx(x)
        # guard r≈0 near axis (should not happen for chamber wall, but robust)
        r = max(1e-12, r)
        return np.tan(gamma) * np.sqrt(1.0 + drdx**2) / r

class InternalPlacement(ChannelPlacement):
    """Placement model for in-wall or in-chamber heat-exchanger channels.

    Allows multiple stacked channels per angular leaf and optional
    user-defined angular-width laws.

    Parameters
    ----------
    n_channel_positions : int
        Number of circumferential leaves.
    n_channels_per_leaf : int
        Channels stacked radially per leaf.
    channel_angular_width : callable
        Function returning angular spacing between channel rows [rad].
    occludes : bool, optional
        Whether this placement occludes the hot-side wall. Default False.

    Attributes
    ----------
    n_channel_positions : int
        Number of circumferential leaves.
    n_channels_per_leaf : int
        Channels per leaf.
    channel_angular_width : callable
        Angular spacing function.
    occludes : bool
        Occlusion flag.

    See Also
    --------
    SurfacePlacement
    ChannelPlacement
    """

    def __init__(self,
                    n_channel_positions: int,          # leaves
                    n_channels_per_leaf: int,          # radial stack in each leaf
                    *,
                    channel_angular_width,                  # row-to-row θ-spacing
                    occludes: bool=False,
                    helix_angle=0.0):
            
        super().__init__(n_channel_positions,
                        channel_angular_width=channel_angular_width,
                        occludes=occludes,
                        helix_angle=helix_angle)

        self.n_channels_per_leaf = n_channels_per_leaf

    def dtheta_dx(self, x, contour, circuit=None):
        gamma = self.helix_angle(x, contour, circuit)
        return np.tan(gamma) * np.sqrt(1 + float(contour.dr_dx(x))**2) / float(contour.r(x))
