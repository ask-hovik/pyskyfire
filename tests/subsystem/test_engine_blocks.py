"""Mass, energy, and signal contracts for engine-cycle blocks."""

from types import SimpleNamespace

import CoolProp.CoolProp as CP
import pytest

from pyskyfire.common import (
    MassFlowMergerBlock,
    MassFlowSplitterBlock,
    PipeBlock,
    RegenBlock,
    SimpleDuctBlock,
    Station,
    TransmissionBlock,
)
from pyskyfire.common import blocks
from pyskyfire.ducts import Pipe, PipeStation, Straight


def test_adiabatic_duct_preserves_enthalpy_and_mass_flow() -> None:
    inlet = Station(1.0e6, 350.0, 2.0)
    block = SimpleDuctBlock("duct", "in", "out", 0.8, "Water")

    stations, signals = block.compute({"in": inlet}, {})
    outlet = stations["out"]
    h_in = CP.PropsSI("H", "P", inlet.p, "T", inlet.T, "Water")
    h_out = CP.PropsSI("H", "P", outlet.p, "T", outlet.T, "Water")

    assert outlet.p == pytest.approx(0.8 * inlet.p)
    assert outlet.mdot == pytest.approx(inlet.mdot)
    assert h_out == pytest.approx(h_in, rel=1e-10)
    assert signals["dp_duct"] == pytest.approx(0.2 * inlet.p)


def test_pipe_block_adapts_the_public_pipe_solver(monkeypatch) -> None:
    pipe = Pipe(
        [PipeStation((0, 0, 0), 0.02), PipeStation((1, 0, 0), 0.02)],
        [Straight()],
    )
    transport = SimpleNamespace(fluid="Water")
    results = [
        SimpleNamespace(p=[3.0e6, 2.8e6], T=[320.0, 320.0], total_pressure_loss=2.0e5),
        SimpleNamespace(p=[3.0e6, 2.75e6], T=[320.0, 320.0], total_pressure_loss=2.5e5),
    ]
    calls = []

    def fake_solve(self, **kwargs):
        calls.append(kwargs)
        return results[len(calls) - 1]

    monkeypatch.setattr(Pipe, "solve", fake_solve)
    block = PipeBlock("feed", "in", "out", pipe, transport, friction_model="native")
    inlet = Station(3.0e6, 320.0, 1.2)

    stations, signals = block.compute({"in": inlet}, {})

    assert block.station_inputs == ["in"]
    assert block.station_outputs == ["out"]
    assert block.signal_inputs == []
    assert block.signal_outputs == ["dp_feed"]
    assert stations == {"out": Station(2.8e6, 320.0, 1.2)}
    assert signals == {"dp_feed": 2.0e5}
    assert calls[0]["inlet_state"].p == inlet.p
    assert calls[0]["inlet_state"].T == inlet.T
    assert calls[0]["mass_flow"] == inlet.mdot
    assert calls[0]["transport"] is transport
    assert calls[0]["friction_model"] == "native"

    assert block.post_process({"in": inlet}, {}) is results[1]


def test_splitter_and_merger_conserve_mass_and_enthalpy() -> None:
    inlet = Station(1.0e6, 350.0, 5.0)
    splitter = MassFlowSplitterBlock(
        "split", "in", ["a", "b"], "Water", fractions=[0.3, 0.7]
    )
    branches, _ = splitter.compute({"in": inlet}, {})

    assert sum(station.mdot for station in branches.values()) == pytest.approx(inlet.mdot)

    branches["a"] = Station(1.0e6, 320.0, branches["a"].mdot)
    branches["b"] = Station(0.9e6, 380.0, branches["b"].mdot)
    merger = MassFlowMergerBlock("merge", ["a", "b"], "out", "Water")
    merged, _ = merger.compute(branches, {})
    outlet = merged["out"]

    h_expected = sum(
        station.mdot
        * CP.PropsSI("Hmass", "T", station.T - 1e-3, "P", station.p, "Water")
        for station in branches.values()
    ) / inlet.mdot
    h_out = CP.PropsSI("Hmass", "T", outlet.T, "P", outlet.p, "Water")
    assert outlet.mdot == pytest.approx(inlet.mdot)
    assert outlet.p == pytest.approx(0.9e6)
    assert h_out == pytest.approx(h_expected, rel=1e-10)


def test_dynamic_split_and_transmission_validate_runtime_signals() -> None:
    splitter = MassFlowSplitterBlock(
        "split", "in", ["a", "b"], "Water", frac_keys=["fa", "fb"]
    )
    with pytest.raises(ValueError, match="do not sum"):
        splitter.compute({"in": Station(1e6, 300.0, 1.0)}, {"fa": 0.4, "fb": 0.5})

    transmission = TransmissionBlock("shaft", ["pump_a", "pump_b"], ["required"])
    _, signals = transmission.compute({}, {"pump_a": 10.0, "pump_b": 15.0})
    assert signals == {"required": 25.0}


def test_regen_block_forwards_inlet_state_and_returns_solver_outlet(monkeypatch) -> None:
    captured = {}
    result = SimpleNamespace(
        T_stagnation=[410.0],
        p_stagnation=[2.7e6],
    )

    def fake_analysis(chamber, **kwargs):
        captured.update(kwargs)
        return result

    monkeypatch.setattr(blocks, "coupled_steady_heating_analysis", fake_analysis)
    nodes = {
        "wall": [0.0, 0.4, 1.0],
        "heat_flux": [0.0, 0.4, 1.0],
        "coolant": [0.0, 0.4, 1.0],
    }
    block = RegenBlock(
        "jacket",
        "in",
        "out",
        1,
        object(),
        "Water",
        nodes=nodes,
        heat_curvature_correction=False,
    )
    inlet = Station(3.0e6, 300.0, 2.5)

    stations, signals = block.compute({"in": inlet}, {})

    boundary = captured["boundary_conditions"]
    assert (boundary.T_coolant_in, boundary.p_coolant_in, boundary.mdot_coolant) == (
        inlet.T,
        inlet.p,
        inlet.mdot,
    )
    assert captured["circuit_index"] == 1
    assert captured["nodes"] == nodes
    assert captured["heat_curvature_correction"] is False
    assert stations["out"] == Station(2.7e6, 410.0, 2.5)
    assert signals == {"dp_jacket": pytest.approx(0.3e6)}
