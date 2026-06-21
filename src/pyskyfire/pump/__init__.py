# __init__.py

# Explicit imports for accessibility
#from .impeller import Impeller
from .impeller_new import Impeller
#from .plot import plot_impeller_views, plot_impeller_3D

# Expose submodules for direct access if needed
from . import impeller, utils, constants

# Define what gets imported when using "import centrifugal_pump as pump"
__all__ = ["Impeller", "plot_impeller_views", "impeller", "utils", "constants"]
