from importlib.util import find_spec


_VIZ_DEPENDENCIES = ("plotly", "pyvista")
_missing_dependencies = [
    dependency
    for dependency in _VIZ_DEPENDENCIES
    if find_spec(dependency) is None
]

if _missing_dependencies:
    missing = ", ".join(_missing_dependencies)
    raise ImportError(
        "Pyskyfire's visualisation tools require optional dependencies "
        f"that are not installed: {missing}. Install them with "
        "`pip install 'pyskyfire[viz]'`."
    )

from .plot_regen import *
from .plot_film_cooling import *
from .report import *
from .core import *
from .lookup_tables import *
from .plot_thermo import *
from .plot_common import *
from .engine_viz import *
from .impeller_viz import *
from .turbine_viz import *
from .network_viz import *
from .sankey import *
from .scene import *
from .components import *
from .heat_exchanger import *
from .heat_exchanger_component import *
from .injector_viz import *
from .plot_ducts import *
