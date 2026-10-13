"""Mass, energy, and signal contracts for engine-cycle blocks."""

from types import SimpleNamespace

import CoolProp.CoolProp as CP
import pytest

from pyskyfire.network import (
    MassFlowMergerBlock,
    InjectorBlock,
    MassFlowSplitterBlock,
    OrificeBlock,
    PipeBlock,
    RegenBlock,
    SimpleDuctBlock,
    Station,
    TransmissionBlock,
    VoluteBlock,
)
from pyskyfire.network import blocks
from pyskyfire.injector import ElementSet, ImpingingElement, Injector, Jet
from pyskyfire.ducts import (
    ConstantAngularMomentum,
    EllipticalSlottedSection,
    ExitCone,
    Orifice,
    OrificePlate,
    Pipe,
    PathStation,
    Straight,
    Volute,
    VoluteInlet,
)


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
        [PathStation((0, 0, 0)), PathStation((1, 0, 0))],
        [Straight()],
    diameter=0.02)
    transport = SimpleNamespace(fluid="Water")
    results = [
        SimpleNamespace(p_stagnation=[3.0e6, 2.8e6], T_stagnation=[320.0, 320.0], total_pressure_loss=2.0e5),
        SimpleNamespace(p_stagnation=[3.0e6, 2.75e6], T_stagnation=[320.0, 320.0], total_pressure_loss=2.5e5),
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
    assert calls[0]["inlet_swirl"] == 0.0

    assert block.post_process({"in": inlet}, {}) is results[1]


def test_pipe_block_carries_a_fixed_inlet_swirl_into_the_pipe() -> None:
    pipe = Pipe([PathStation((0, 0, 0)), PathStation((1, 0, 0))], [Straight()],
                diameter=0.04, core_diameter=0.02)
    transport = SimpleNamespace(fluid="Water", get_rho=lambda T, p: 1000.0,
                                get_mu=lambda T, p: 1.0e-3, get_cp=lambda T, p: 4180.0)
    inlet = Station(3.0e6, 320.0, 2.0)
    plain = PipeBlock("feed", "in", "out", pipe, transport)
    swirling = PipeBlock("feed", "in", "out", pipe, transport, inlet_swirl=5.0)

    result = swirling.post_process({"in": inlet}, {})
    _, signals = swirling.compute({"in": inlet}, {})

    assert result.swirl_velocity[0] == pytest.approx(5.0)
    assert signals["dp_feed"] == pytest.approx(result.total_pressure_loss)
    assert signals["dp_feed"] > plain.compute({"in": inlet}, {})[1]["dp_feed"]


def test_orifice_block_emits_the_orifice_loss_at_constant_enthalpy() -> None:
    orifice = Orifice.circular(2e-3, cd=0.7)
    block = OrificeBlock("inj", "in", "out", orifice, "Water")
    inlet = Station(50e5, 300.0, 0.05)

    stations, signals = block.compute({"in": inlet}, {})
    outlet = stations["out"]
    result = block.post_process({"in": inlet}, {})

    assert block.station_inputs == ["in"]
    assert block.station_outputs == ["out"]
    assert block.signal_outputs == ["dp_inj"]
    assert block.medium == "Water"
    assert signals["dp_inj"] == pytest.approx(result.pressure_loss)
    assert outlet.p == pytest.approx(inlet.p - signals["dp_inj"])
    assert outlet.mdot == inlet.mdot
    h_in = CP.PropsSI("H", "P", inlet.p, "T", inlet.T, "Water")
    h_out = CP.PropsSI("H", "P", outlet.p, "T", outlet.T, "Water")
    assert h_out == pytest.approx(h_in)


def test_orifice_block_books_a_plates_permanent_loss() -> None:
    plate = OrificePlate(0.1, 0.05, cd=0.6)
    block = OrificeBlock("meter", "in", "out", plate, "Water")

    stations, signals = block.compute({"in": Station(10e5, 300.0, 20.0)}, {})
    result = block.post_process({"in": Station(10e5, 300.0, 20.0)}, {})

    assert signals["dp_meter"] == pytest.approx(result.pressure_loss)
    assert signals["dp_meter"] < result.differential_pressure
    assert stations["out"].p == pytest.approx(10e5 - signals["dp_meter"])


def test_orifice_block_rejects_viscosity_for_a_bare_orifice() -> None:
    with pytest.raises(ValueError, match="viscosity"):
        OrificeBlock("inj", "in", "out", Orifice.circular(2e-3, cd=0.7), "Water",
                     viscosity=1e-3)


def test_injector_block_solves_one_feed_to_the_face() -> None:
    element = ImpingingElement([
        Jet("fuel", Orifice.circular(1e-3, cd=0.7), 0.5),
        Jet("ox", Orifice.circular(1.2e-3, cd=0.7), 0.5, 3.14159),
    ], impingement_height=5e-3)
    injector = Injector([ElementSet(element, [0.01, 0.02], [0.0, 1.0], name="main")])
    block = InjectorBlock("inj_ox", "ox_plenum", "ox_face", injector, "ox", "Water")
    plenum = Station(20e5, 300.0, 0.02)

    stations, signals = block.compute({"ox_plenum": plenum}, {})
    result = block.post_process({"ox_plenum": plenum}, {})

    assert block.signal_outputs == ["dp_inj_ox"]
    assert signals["dp_inj_ox"] == pytest.approx(result.pressure_drop)
    assert stations["ox_face"].p == pytest.approx(plenum.p - signals["dp_inj_ox"])
    assert stations["ox_face"].mdot == plenum.mdot
    assert result.set_mass_flows == {"main": pytest.approx(0.02)}


def test_injector_block_rejects_a_feed_the_injector_lacks() -> None:
    element = ImpingingElement([Jet("fuel", Orifice.circular(1e-3, cd=0.7), 0.5)],
                               impingement_height=5e-3)
    injector = Injector([ElementSet(element, [0.01], [0.0], name="main")])
    with pytest.raises(ValueError, match="feed"):
        InjectorBlock("inj", "in", "out", injector, "film", "Water")


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


def test_volute_block_adapts_the_public_volute_solver(monkeypatch) -> None:
    volute = Volute(
        EllipticalSlottedSection(inlet_width=0.012, fillet_radius=0.003),
        VoluteInlet(radius=0.060, C_m=6.0, C_u=45.0),
        ConstantAngularMomentum(),
        rotation="ccw",
        n_stations=13,
    )
    transport = SimpleNamespace(fluid="Water")
    results = [
        SimpleNamespace(
            states=(SimpleNamespace(p=3.0e6), SimpleNamespace(p=2.8e6)),
            total_pressure_loss=2.0e5,
        ),
        SimpleNamespace(
            states=(SimpleNamespace(p=3.0e6), SimpleNamespace(p=2.9e6)),
            total_pressure_loss=1.0e5,
        ),
    ]
    calls = []

    def fake_solve(self, **kwargs):
        calls.append(kwargs)
        return results[len(calls) - 1]

    monkeypatch.setattr(Volute, "solve", fake_solve)
    block = VoluteBlock(
        "discharge", "in", "out", volute, transport,
        exit_diameter=0.05, roughness=1.0e-5, friction_model="native",
    )
    inlet = Station(3.0e6, 320.0, 1.2)

    stations, signals = block.compute({"in": inlet}, {})

    assert block.station_inputs == ["in"]
    assert block.station_outputs == ["out"]
    assert block.signal_inputs == []
    assert block.signal_outputs == ["dp_discharge"]
    assert stations == {"out": Station(2.8e6, 320.0, 1.2)}
    assert signals == {"dp_discharge": 2.0e5}
    assert calls[0]["inlet_state"].p == inlet.p
    assert calls[0]["inlet_state"].T == inlet.T
    assert calls[0]["mass_flow"] == inlet.mdot
    assert calls[0]["transport"] is transport
    assert calls[0]["exit_diameter"] == 0.05
    assert calls[0]["roughness"] == 1.0e-5
    assert calls[0]["friction_model"] == "native"

    assert block.post_process({"in": inlet}, {}) is results[1]


def test_volute_block_defaults_exit_diameter_from_an_attached_exit_cone(monkeypatch) -> None:
    cone = ExitCone(final_diameter=0.04, first_angle="tangential")
    volute = Volute(
        EllipticalSlottedSection(inlet_width=0.012, fillet_radius=0.003),
        VoluteInlet(radius=0.060, C_m=6.0, C_u=45.0),
        ConstantAngularMomentum(),
        rotation="ccw",
        exit_cone=cone,
        n_stations=13,
    )
    transport = SimpleNamespace(fluid="Water")
    calls = []

    def fake_solve(self, **kwargs):
        calls.append(kwargs)
        return SimpleNamespace(
            states=(SimpleNamespace(p=3.0e6), SimpleNamespace(p=2.8e6)),
            total_pressure_loss=2.0e5,
        )

    monkeypatch.setattr(Volute, "solve", fake_solve)
    block = VoluteBlock("discharge", "in", "out", volute, transport)
    block.compute({"in": Station(3.0e6, 320.0, 1.2)}, {})
    assert calls[0]["exit_diameter"] is None


def test_volute_block_requires_exit_diameter_without_an_exit_cone() -> None:
    volute = Volute(
        EllipticalSlottedSection(inlet_width=0.012, fillet_radius=0.003),
        VoluteInlet(radius=0.060, C_m=6.0, C_u=45.0),
        ConstantAngularMomentum(),
        rotation="ccw",
        n_stations=13,
    )
    with pytest.raises(ValueError, match="exit_diameter is required"):
        VoluteBlock("discharge", "in", "out", volute, SimpleNamespace(fluid="Water"))
