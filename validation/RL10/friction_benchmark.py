"""Compare coolant friction backends against the Binder et al. RL10 data.

Runs the full RL10 regen case once per friction model and overlays the
resulting coolant static-pressure curves on the digitised reference data. The
two models share everything except the Darcy friction factor and its curvature
multiplier, so the gap between the curves isolates the friction correlation's
contribution to the predicted pressure drop.

Run directly to produce ``standalone/friction-benchmark.html``::

    python validation/RL10/friction_benchmark.py
"""

import json
import os

import numpy as np
import plotly.graph_objects as go

import regen_sim

script_dir = os.path.dirname(os.path.abspath(__file__))
reference_dir = os.path.join(script_dir, "reference_data")
standalone_dir = os.path.join(script_dir, "standalone")

# The RL10A-3-3 design-report curve is a different engine build and is excluded
# from the comparison plots for the same reason post_process.py excludes it.
EXCLUDED_REFERENCE_MODELS = ("PWA_FR_1769",)

MODEL_STYLE = {
    "native": dict(color="#2563eb", dash="solid"),
    "fluids": dict(color="#d97706", dash="dash"),
}

MODEL_LABEL = {
    "native": "pyskyfire native (RL10 / Ito)",
    "fluids": "fluids (Clamond + helical-coil)",
}


def load_reference_pressure():
    """Return ``[(label, x, p_static), ...]`` from the digitised reference."""
    path = os.path.join(reference_dir, "reference_coolant_static_pressure.json")
    with open(path, "r", encoding="utf-8") as handle:
        data = json.load(handle)

    series = []
    for key, entry in data.items():
        if key == "units" or key in EXCLUDED_REFERENCE_MODELS:
            continue
        series.append(
            (
                entry.get("name", key),
                np.asarray(entry["x"], dtype=float),
                np.asarray(entry["p_static"], dtype=float),
            )
        )
    return series


def circuit_pressure_drop(cooling_data):
    """Return total stagnation pressure drop across a circuit [Pa]."""
    return float(cooling_data.p_stagnation[0] - cooling_data.p_stagnation[-1])


def build_figure(runs):
    """Overlay every friction model's circuit-B static pressure on the data."""
    fig = go.Figure()

    for label, x, p in load_reference_pressure():
        fig.add_trace(
            go.Scatter(
                x=x,
                y=p / 1e5,
                mode="markers",
                name=label,
                marker=dict(size=8, symbol="circle-open", line=dict(width=2)),
            )
        )

    for model, cooling_data in runs.items():
        # Circuit B is the full-length pass the reference curve covers.
        result = cooling_data[1]
        fig.add_trace(
            go.Scatter(
                x=np.asarray(result.x_coolant),
                y=np.asarray(result.p_static) / 1e5,
                mode="lines",
                name=MODEL_LABEL.get(model, model),
                line=dict(width=2.5, **MODEL_STYLE.get(model, {})),
            )
        )

    fig.update_layout(
        template="plotly_white",
        title=(
            "RL10 coolant static pressure: friction backends vs. Binder et al."
        ),
        xaxis=dict(title="Axial position (m)"),
        yaxis=dict(title="Static pressure (bar)"),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0),
        margin=dict(l=70, r=30, t=90, b=60),
    )
    return fig


def report_pressure_drops(runs):
    """Print the per-circuit and total pressure drop for each model."""
    print("\n" + "=" * 66)
    print("Stagnation pressure drop by friction model")
    print("=" * 66)
    header = f"{'model':<10}{'circuit A':>14}{'circuit B':>14}{'total':>14}"
    print(header)
    print("-" * 66)

    totals = {}
    for model, cooling_data in runs.items():
        dp_a = circuit_pressure_drop(cooling_data[0])
        dp_b = circuit_pressure_drop(cooling_data[1])
        totals[model] = dp_a + dp_b
        print(
            f"{model:<10}{dp_a / 1e5:>12.3f} bar{dp_b / 1e5:>12.3f} bar"
            f"{(dp_a + dp_b) / 1e5:>12.3f} bar"
        )

    if "native" in totals and "fluids" in totals:
        native, fluid = totals["native"], totals["fluids"]
        print("-" * 66)
        print(
            f"fluids / native = {fluid / native:.4f}  "
            f"({(fluid - native) / 1e5:+.3f} bar)"
        )
    print("=" * 66 + "\n")
    return totals


def main():
    chamber, runs = run_benchmark()
    report_pressure_drops(runs)

    fig = build_figure(runs)
    os.makedirs(standalone_dir, exist_ok=True)
    out_path = os.path.join(standalone_dir, "friction-benchmark.html")
    fig.write_html(out_path, include_plotlyjs="cdn")
    print(f"wrote {out_path}")
    return chamber, runs, fig


def run_benchmark():
    """Run the RL10 case with every friction backend and save the results."""
    chamber, runs = regen_sim.run_friction_benchmark(use_reference_nodes=True)
    regen_sim.save_friction_benchmark(chamber, runs)
    return chamber, runs


if __name__ == "__main__":
    main()
