from abc import ABC, abstractmethod
import numpy as np

class ChannelPlacement(ABC):
    """Abstract base class defining how coolant channels are positioned.

    Subclasses implement different placement strategies (surface, internal)
    that compute the radial coordinate of the channel centerline.

    Parameters
    ----------
    n_channel_positions : int
        Number of circumferential channel locations (“leaves”).
    channel_width : callable or None, optional
        Function returning angular width [rad]. May be `None` if uniform.
    occludes : bool, optional
        Whether this placement blocks part of the wall from hot gas view.

    Attributes
    ----------
    n_channel_positions : int
        Number of circumferential positions.
    channel_width : callable or None
        Angular width function or constant.
    occludes : bool
        Whether the placement occludes hot surface area.

    See Also
    --------
    SurfacePlacement
    InternalPlacement
    """
    def __init__(self, n_channel_positions: int, channel_width=None, occludes: bool = True):
        self.n_channel_positions = n_channel_positions   # << new home
        self.channel_width = channel_width
        self.occludes = occludes

    @abstractmethod
    def compute_centerline_radius(self,
                                  x: float,
                                  contour) -> float:
        """Return the centerline radius r(x) for this placement."""
        raise NotImplementedError

    def channel_count(self) -> int:
        return self.n_channel_positions
    
    def dtheta_dx(self, x: float, contour) -> float:
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
    def __init__(self, n_channel_positions: int, helix_angle = lambda *args: 0.0):
        self.n_channels_per_leaf = 1
        super().__init__(n_channel_positions, channel_width = None, occludes=True)
        self.helix_angle = helix_angle  # float or callable(x)->float

    def compute_centerline_radius(self, x, contour):
        r_hot   = contour.r(x)
        #alpha       = contour.normal_angle(x)
        #t_total = wall_group.total_thickness(x)
        #TODO: This is approximately correct. Will want to reimplement the effect of the channel wall thickness in the future. 
        return r_hot #+ t_total/np.cos(alpha)
    
    def dtheta_dx(self, x, contour):
        gamma = self.helix_angle(x) if callable(self.helix_angle) else float(self.helix_angle)
        r = self.compute_centerline_radius(x, contour)
        drdx = contour.dr_dx(x)
        # guard r≈0 near axis (should not happen for chamber wall, but robust)
        r = max(1e-12, r)
        return np.tan(gamma) * np.sqrt(1.0 + drdx**2) / r

class InternalPlacement(ChannelPlacement):
    """Placement model for in-wall or in-chamber heat-exchanger channels.

    Allows multiple stacked channels per angular leaf and optional
    user-defined width laws.

    Parameters
    ----------
    n_channel_positions : int
        Number of circumferential leaves.
    n_channels_per_leaf : int
        Channels stacked radially per leaf.
    channel_width : callable
        Function returning angular spacing between channel rows [rad].
    occludes : bool, optional
        Whether this placement occludes the hot-side wall. Default False.

    Attributes
    ----------
    n_channel_positions : int
        Number of circumferential leaves.
    n_channels_per_leaf : int
        Channels per leaf.
    channel_width : callable
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
                    channel_width,                          # row-to-row θ-spacing
                    occludes: bool=False):
            
        super().__init__(n_channel_positions,
                        channel_width=channel_width,
                        occludes=occludes)
        
        self.n_channels_per_leaf = n_channels_per_leaf

    def compute_centerline_radius(self, x, contour):
        return 1 # TODO: I need a better way to implement this

    # override: total = leaves ×  channels / leaf
    """def channel_count(self) -> int:
        return self.n_channel_positions# * self.n_channels_per_leaf"""
