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
from .plot_skycea import *
from .plot_common import *
from .engine_viz import *
from .impeller_viz import *
from .network_viz import *
from .sankey import *
