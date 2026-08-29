"""Fixed-point orchestration contracts for engine networks."""

from types import SimpleNamespace

import pytest

from pyskyfire.common import EngineNetwork, Station


class _ConstantBlock:
    name = "constant"

    def compute(self, stations, signals):
        return {"out": Station(2.0e6, 350.0, 3.0)}, {"target": 4.0}

    def post_process(self, stations, signals):
        return {"final_target": signals["target"]}


class _OscillatingBlock:
    name = "oscillator"

    def compute(self, stations, signals):
        return {}, {"value": 1.0 - signals["value"]}

    def post_process(self, stations, signals):
        return {}


def test_fixed_point_converges_and_collects_post_process_results(capsys) -> None:
    network = EngineNetwork(
        stations={"out": Station(1.0e6, 300.0, 1.0)},
        signals={"target": 1.0},
        blocks=[_ConstantBlock()],
    )

    network.run_fixed_point(tol=1e-12, max_iter=3)

    assert len(network.residuals) == 2
    assert network.residuals[-1] == pytest.approx(0.0)
    assert network.stations["out"] == Station(2.0e6, 350.0, 3.0)
    assert network.block_results == {"constant": {"final_target": 4.0}}
    assert "Converged in 2 iterations" in capsys.readouterr().out


def test_fixed_point_reports_nonconvergence() -> None:
    network = EngineNetwork(
        stations={},
        signals={"value": 0.0},
        blocks=[_OscillatingBlock()],
    )

    with pytest.raises(RuntimeError, match="Did not converge after 3 iterations"):
        network.run_fixed_point(tol=1e-12, max_iter=3)

    assert len(network.residuals) == 3


def test_station_residual_tracks_pressure_temperature_and_mass_flow() -> None:
    network = EngineNetwork({}, {}, [])
    old = Station(100.0, 200.0, 4.0)
    store = {"station": old}

    residual = network._merge_and_residual(
        store,
        {"station": Station(110.0, 180.0, 5.0)},
        0.0,
    )

    assert residual == pytest.approx(0.25)
