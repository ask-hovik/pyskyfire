"""Scene and face-view contracts for injectors."""

import math

import numpy as np
import pytest

pytest.importorskip("plotly")
pytest.importorskip("pyvista")

from pyskyfire.ducts import Orifice
from pyskyfire.injector import ElementSet, FeedCondition, ImpingingElement, Injector, Jet
from pyskyfire.thermo import Incompressible
from pyskyfire.viz import FEED_COLORS, InjectorComponent, Scene, plot_injector_face


def injector():
    element = ImpingingElement([
        Jet("fuel", Orifice.circular(1e-3, cd=0.7, length=4e-3), math.radians(30)),
        Jet("ox", Orifice.circular(1.2e-3, cd=0.7), math.radians(30), math.pi),
    ], impingement_height=5e-3)
    return Injector([ElementSet(element, [0.01, 0.02, 0.02], [0.0, 0.0, math.pi], name="main")])


def test_holes_are_grouped_by_set_and_feed_and_coloured_by_feed():
    root = InjectorComponent(injector(), face_radius=0.03, n_theta=8).build()

    assert [child.name for child in root.children] == ["face", "main"]
    main = root.children[1]
    # Only wetted surface by default: no jets, no arrows.
    assert [child.name for child in main.children] == ["fuel bores", "ox bores"]
    assert main.children[0].color == FEED_COLORS[0]
    assert main.children[1].color == FEED_COLORS[1]
    # Three holes per feed, one open tube of 2 × 8 points each.
    assert len(main.children[0].mesh.points) == 3*16
    np.testing.assert_allclose(root.anchor("face").position, [0, 0, 0])


def test_bores_end_on_the_face_plate():
    inj = injector()
    bores = InjectorComponent(inj, n_theta=8).build().children[0].children[0].mesh.points
    holes = inj.holes()

    # Nothing wetted pokes into the chamber.
    assert bores[:, 0].max() == pytest.approx(0.0, abs=1e-15)
    for k in range(3):
        inlet, face = bores[16*k:16*k + 8], bores[16*k + 8:16*k + 16]
        # The downstream ring lies in the face plane, centred on the exit.
        np.testing.assert_allclose(face[:, 0], 0.0, atol=1e-15)
        np.testing.assert_allclose(face.mean(axis=0), holes.exit[2*k], atol=1e-12)
        # The upstream end is square to the axis, the orifice length back.
        np.testing.assert_allclose(inlet.mean(axis=0),
                                   holes.exit[2*k] - 4e-3*holes.direction[2*k], atol=1e-12)
        np.testing.assert_allclose((inlet - inlet.mean(axis=0)) @ holes.direction[2*k], 0.0,
                                   atol=1e-15)


def test_impingement_arrows_are_optional_and_point_at_the_target():
    inj = injector()
    main = InjectorComponent(inj, show_impingement=True, n_theta=8).build().children[0]
    names = [child.name for child in main.children]
    assert names == ["fuel bores", "fuel impingement", "ox bores", "ox impingement"]

    arrows = main.children[1]
    assert arrows.color == "dimgray" and not arrows.fields
    tips = arrows.mesh.points
    # Every arrow ends at the impingement point and starts at the face.
    assert tips[:, 0].max() == pytest.approx(5e-3)
    assert tips[:, 0].min() >= -1e-4
    np.testing.assert_allclose(tips[np.argmax(tips[:, 0])], inj.holes().target[0], atol=1e-12)


def test_a_result_becomes_hole_flow_fields():
    inj = injector()
    water = Incompressible(1000.0)
    result = inj.solve({"fuel": FeedCondition(water, 300.0, pressure=12e5),
                        "ox": FeedCondition(water, 300.0, pressure=13e5)}, p_face=10e5)
    fuel = InjectorComponent(inj, result=result, n_theta=8).build().children[0].children[0]

    per_hole = result.feeds["fuel"].mass_flow/3
    assert fuel.fields["hole_mass_flow"].values == pytest.approx(np.full(3*16, per_hole))
    assert fuel.fields["jet_velocity"].values[0] == pytest.approx(result.feeds["fuel"].velocity)
    scene = Scene()
    scene.add(InjectorComponent(inj, result=result), name="injector")
    assert "hole_mass_flow" in scene.to_html()


def test_more_feeds_than_colours_need_explicit_colours():
    element = ImpingingElement([Jet(f"f{i}", Orifice.circular(1e-3, cd=0.7), 0.3, i)
                                for i in range(len(FEED_COLORS) + 1)], impingement_height=3e-3)
    inj = Injector([ElementSet(element, [0.01], [0.0], name="many")])
    with pytest.raises(ValueError, match="colors"):
        InjectorComponent(inj).build()
    colors = {feed: "black" for feed in inj.feeds}
    assert InjectorComponent(inj, colors=colors).build().children[0].children[0].color == "black"


def test_face_view_has_one_legend_entry_per_feed_and_hover_per_hole():
    fig = plot_injector_face(injector(), face_radius=0.03)
    names = [trace.name for trace in fig.data if trace.showlegend is not False]
    assert names == ["Chamber wall", "fuel", "ox", "Impingement points"]
    hover = [trace for trace in fig.data if trace.hovertemplate]
    assert sum(len(trace.x) for trace in hover) == 6
    assert "main, element 1" in hover[0].text[0]


def test_each_feed_records_the_fluid_it_ran_on():
    inj = injector()
    result = inj.solve({"fuel": FeedCondition("Ethanol", 300.0, pressure=12e5),
                        "ox": FeedCondition(Incompressible(1000.0), 300.0, pressure=13e5)},
                       p_face=10e5)
    assert result.feeds["fuel"].fluid == "Ethanol"
    # A constant density names no fluid.
    assert result.feeds["ox"].fluid is None
    fuel, ox = InjectorComponent(inj, result=result, n_theta=8).build().children[0].children
    assert fuel.metadata["fluid"] == "Ethanol"
    assert "fluid" not in ox.metadata
