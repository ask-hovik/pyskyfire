"""Shared scene placement, mesh reuse, and surface-specific exchanger fields."""

from types import SimpleNamespace
import json
import numpy as np
import pytest

from pyskyfire import common, ducts, regen, heat_exchanger as hx
from pyskyfire.viz import HeatExchangerComponent, Scene, SceneNode, Transform
from pyskyfire.viz.backends.threejs import scene_payload


def case():
    path = ducts.Path.from_nodes(
        [
            ducts.PathNode((0.2, 0, 0)),
            ducts.PathNode((0.6, 0, 0), bend=ducts.Continuous(0.1)),
            ducts.PathNode((0.6, 0.4, 0)),
        ]
    )
    phi = 0.7
    placements = [
        ducts.PathPlacement((0, 0, 0)),
        ducts.PathPlacement(
            (0, 0.2, 0.3),
            rotation=ducts.rotation_about_axis((1, 0, 0), phi),
            initial_floor_direction=(0, -np.cos(phi), -np.sin(phi)),
        ),
        ducts.PathPlacement((0, 0, 0.4), initial_floor_direction=(0, 0, -1)),
    ]
    c = hx.HeatExchangerCircuit(
        name="stream",
        path=path,
        placements=placements,
        cross_section=regen.CrossSectionSquared(
            0.001, floor="line", roof="line", exposure={"floor": 0.5, "roof": 1}
        ),
        wall=regen.Wall(common.Material("metal", k=common.ConstantModel(30.0)), 0.001),
        transport=SimpleNamespace(fluid="fixture"),
        dimensions=dict(height=0.006, floor_width=0.01, roof_width=0.012),
    )
    n = len(path.stations)
    states = tuple(
        SimpleNamespace(T=300 + 10 * i, p=1e6 - i * 100, h=1000 * i) for i in range(n)
    )
    r = hx.CircuitResult(
        c, states, 0.3, np.arange(n, dtype=float), *(np.ones(n - 1) for _ in range(6))
    )
    walls = tuple(
        hx.WallResult(
            "stream", i, "floor", None, None, None, 400 + i, 600 + i, 10.0, 0.001
        )
        for i in range(n - 1)
    )
    result = hx.HeatExchangerResult({"stream": r}, walls, 2, 0.0, 0.0, 0.0)
    return hx.HeatExchanger([c]), c, result


def test_instances_draw_the_wetted_passage_and_only_share_equal_local_roll():
    exchanger, c, result = case()
    scene = Scene(color_by="T_coolant_wall")
    root = scene.add(HeatExchangerComponent(exchanger, result, samples=4))
    channels = root.children[0].children
    assert channels[0].mesh is channels[1].mesh
    assert channels[0].mesh is not channels[2].mesh
    for i, node in enumerate(channels):
        np.testing.assert_allclose(
            node.world_transform.apply(node.mesh.points),
            c.section_points(i, exterior=False, samples=4).reshape(-1, 3),
            atol=1e-12,
        )
        np.testing.assert_allclose(
            node.anchor("inlet").position, c.placed_paths[i].centerline[0], atol=1e-12
        )
        np.testing.assert_allclose(
            node.anchor("outlet").position, c.placed_paths[i].centerline[-1], atol=1e-12
        )
    np.testing.assert_allclose(
        root.anchor("circuit_0.inlet").position, c.placed_paths[0].centerline[0]
    )
    assert len(scene_payload(scene)["meshes"]) == 2
    json.dumps(scene_payload(scene), allow_nan=False)


def test_wall_temperatures_follow_segment_surface_and_partial_exposure():
    exchanger, c, result = case()
    root = HeatExchangerComponent(exchanger, result, samples=4).build()
    node = root.children[0].children[0]
    np.testing.assert_allclose(
        node.fields["T_static"].values.reshape(-1, 16)[:, 0],
        result.circuits["stream"].T,
    )
    assert "T_hot_wall" not in node.fields
    wall = node.fields["T_coolant_wall"]
    assert wall.association == "cell" and wall.units == "K"
    values = wall.values.reshape(-1, 4, 4, 2)
    np.testing.assert_allclose(values[:, 0, 1, 0], 400 + np.arange(len(c.s) - 1))
    assert np.isnan(values[:, 0, [0, 3], :]).all()
    assert np.isnan(values[:, 1:, :, :]).all()


def test_installation_and_parent_transform_are_applied_once():
    exchanger, c, result = case()
    installation = hx.ChamberInstallation.__new__(hx.ChamberInstallation)
    installation.exchanger = exchanger
    installation.gas = SimpleNamespace(
        rotation=ducts.rotation_about_axis((0, 0, 1), 0.4),
        position=np.array([0.7, 0.8, 0.9]),
    )
    scene = Scene()
    parent = scene.add(
        SimpleNamespace(build=lambda: SceneNode("parent")),
        position=(2, 3, 4),
        rotation=(0, 0.3, 0),
    )
    root = scene.add(HeatExchangerComponent(installation), parent=parent)
    node = root.children[0].children[1]
    actual = node.world_transform.apply(node.mesh.points)
    expected = (
        c.section_points(1, exterior=False, samples=8).reshape(-1, 3) @ installation.gas.rotation.T
        + installation.gas.position
    )
    np.testing.assert_allclose(
        actual, parent.world_transform.apply(expected), atol=1e-12
    )
    assert not node.fields
    # Scene.add explicit placement replaces the installed root placement.
    raw = Scene().add(HeatExchangerComponent(installation), position=(10, 0, 0))
    np.testing.assert_allclose(raw.transform.matrix, Transform((10, 0, 0)).matrix)


def test_result_identity_and_resolution_are_validated():
    exchanger, c, result = case()
    different, _, _ = case()
    with pytest.raises(ValueError, match="different circuit"):
        HeatExchangerComponent(different, result).build()
    with pytest.raises(ValueError, match="samples"):
        HeatExchangerComponent(exchanger, samples=1).build()
