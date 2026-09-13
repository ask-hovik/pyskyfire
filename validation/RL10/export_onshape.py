"""RL10 entry point for the shared example Onshape exporter."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
_SPEC = spec_from_file_location(
    "pyskyfire_example_onshape_export",
    _SCRIPT_DIR.parents[1] / "examples" / "advanced" / "export_onshape.py",
)
_EXPORTER = module_from_spec(_SPEC)
_SPEC.loader.exec_module(_EXPORTER)
build_onshape_data = _EXPORTER.build_onshape_data
export_onshape_json = _EXPORTER.export_onshape_json

if __name__ == "__main__":
    _EXPORTER.main(
        script_dir=_SCRIPT_DIR,
        results_filename="regen_results.pkl",
        output_filename="RL10-cooling-channels.json",
    )
