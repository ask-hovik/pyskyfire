"""Scene contracts: placement, result mapping, reusable seeds, and both outputs."""
import base64
import json
from importlib.resources import files
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest

from pyskyfire import ducts
from pyskyfire.regen.channel_layout import ChannelLayout, ResolvedChannel
from pyskyfire.regen.cross_section import CrossSectionSquared, CrossSectionRounded
from pyskyfire.viz import (
    Anchor, FieldData, ImpellerComponent, MeshData, PipeComponent, Scene,
    SceneNode, ThrustChamberComponent, Transform, VoluteComponent, duct_fields,
    volute_fields,
)
from pyskyfire.viz.backends.threejs import scene_payload


def chamber_case(section=None):
    x = np.linspace(0., 0.04, 9)
    path = np.column_stack((x, 0.1+0.2*x, 2*x))
    circuit = SimpleNamespace(
        name="Fuel", direction=-1, reference_path=path, x_domain=x,
        cross_section=section or CrossSectionSquared((0.001) / 2),
        channel_heights=np.full_like(x, 0.012),
        channel_local_sector=np.full_like(x, 0.1),
        channel_wall_thickness=lambda v: np.full_like(v, 0.001),
        hot_side_layer_thickness=lambda v: np.zeros_like(v),
    )
    chamber = SimpleNamespace(cooling_circuits=[circuit])
    seed = ResolvedChannel("fuel", 0, 0, path, 6, 2*np.pi/6)
    layout = ChannelLayout("analytical", (seed,))
    xc = np.array([0.04, 0.02, 0.0])
    xw = np.linspace(0, 0.04, 5)
    result = SimpleNamespace(
        circuit_index=0, circuit_name="Fuel", x=xw, x_coolant=xc,
        p_static=1e6+1e7*xc, T_static=300+100*xc,
        x_wall=xw, T=np.column_stack((300+xw, 400+xw, 900+xw)),
        x_heat_flux=np.array([0., 0.04]), qpp_hot=np.array([1e6, 2e6]),
    )
    return chamber, layout, result


def pipe():
    return SimpleNamespace(centerline=np.array([[0., 0, 0], [1., 0, 0], [1., 1., 0]]), radius=0.01)


@pytest.mark.parametrize("section", [CrossSectionSquared((0.001) / 2), CrossSectionRounded()])
def test_each_seed_is_built_once_and_rotations_match_legacy(monkeypatch, section):
    from pyskyfire.viz import components
    from pyskyfire.viz.engine_viz import build_sections_for_centerlines
    chamber, layout, _ = chamber_case(section)
    build = Mock(wraps=components.build_channel_construction)
    monkeypatch.setattr(components, "build_channel_construction", build)
    scene = Scene()
    root = scene.add(ThrustChamberComponent(chamber, layout=layout, stride=3))
    channels = root.children[0].children
    assert build.call_count == 1
    assert len(channels) == 6
    assert all(c.mesh is channels[0].mesh for c in channels)
    assert len(channels[0].mesh.points) == 4*16  # includes final partial stride
    expanded = layout.expanded_channels()
    reference = build_sections_for_centerlines([c.reference_path for c in expanded], [0]*6, chamber, stride=3)
    for channel, sections in zip(channels, reference, strict=True):
        np.testing.assert_allclose(channel.world_transform.apply(channel.mesh.points), np.vstack(sections), atol=1e-12)
    assert len(scene_payload(scene)["meshes"]) == 1


def test_fields_use_native_grids_flow_direction_and_wall_column_order():
    chamber, layout, result = chamber_case()
    scene = Scene(color_by="p_static")
    root = scene.add(ThrustChamberComponent(chamber, result, layout=layout, stride=3))
    channel = root.children[0].children[0]
    x = channel.mesh.points[:, 0]
    np.testing.assert_allclose(channel.fields["p_static"].values, 1e6+1e7*x)
    np.testing.assert_allclose(channel.fields["T_hot_wall"].values, 900+x)
    np.testing.assert_allclose(channel.fields["T_coolant_wall"].values, 400+x)
    np.testing.assert_allclose(channel.fields["qpp_hot"].values, 1e6+1e6*x/0.04)
    assert channel.fields["p_static"].units == "Pa"
    inlet, outlet = root.anchor("circuit_0.inlet"), root.anchor("circuit_0.outlet")
    assert inlet.position[0] == pytest.approx(0.04)
    assert outlet.position[0] == pytest.approx(0.)
    assert inlet.transform.matrix[0, 0] > 0
    assert outlet.transform.matrix[0, 0] < 0


def test_missing_data_is_not_extrapolated_or_interpolated_across_nan():
    chamber, layout, result = chamber_case()
    result.x_coolant = np.array([0.01, 0.02, 0.03])
    result.p_static = np.array([1., np.nan, 3.])
    scene = Scene()
    root = scene.add(ThrustChamberComponent(chamber, result, layout=layout, stride=1))
    values = root.children[0].children[0].fields["p_static"].values[::16]
    assert np.isnan(values[[0, 1, 3, 4, 5, 7, 8]]).all()
    np.testing.assert_allclose(values[[2, 6]], [1, 3])
    json.dumps(scene_payload(scene), allow_nan=False)


@pytest.mark.parametrize("form", [lambda r: r, lambda r: [r], lambda r: {"arbitrary": r}])
def test_result_input_forms_use_circuit_identity(form):
    chamber, layout, result = chamber_case()
    root = ThrustChamberComponent(chamber, form(result), layout=layout).build()
    assert "p_static" in root.children[0].children[0].fields


@pytest.mark.parametrize("problem", ["duplicate", "name", "index", "grid"])
def test_invalid_result_associations_are_rejected(problem):
    chamber, layout, result = chamber_case()
    data = result
    if problem == "duplicate":
        data = [result, result]
    elif problem == "name":
        result.circuit_name = "Wrong"
    elif problem == "index":
        result.circuit_index = 1
    else:
        result.x_coolant = [0, 0, 0.04]
    scene = Scene()
    with pytest.raises(ValueError):
        scene.add(ThrustChamberComponent(chamber, data, layout=layout))
    assert not scene.nodes


def test_transform_hierarchy_and_visual_anchor_mating():
    scene = Scene()
    assembly = scene.add(SimpleNamespace(build=lambda: SceneNode("assembly")), position=(2, 3, 4), rotation=(0, 0, np.pi/2))
    first = scene.add(PipeComponent(pipe()), parent=assembly, name="first", position=(1, 0, 0))
    np.testing.assert_allclose(first.world_transform.apply([0, 0, 0]), [2, 4, 4], atol=1e-12)
    second = scene.add(PipeComponent(pipe()), name="second", parent=assembly,
                       anchor="inlet", attach_to=first.anchor("outlet"))
    np.testing.assert_allclose(second.anchor("inlet").position, first.anchor("outlet").position, atol=1e-12)
    np.testing.assert_allclose(second.anchor("inlet").transform.matrix[:3, 0],
                               -first.anchor("outlet").transform.matrix[:3, 0], atol=1e-12)
    assert second.id == "assembly/second"
    # No relationship to fluid connectivity and no placement requirement.
    scene.add(PipeComponent(pipe()), name="unplaced")
    scene.add(PipeComponent(pipe()), name="also_unplaced")


@pytest.mark.parametrize("kwargs", [
    {"transform": Transform(), "position": (0, 0, 0)},
    {"anchor": "inlet"},
    {"anchor": "inlet", "attach_to": Anchor(), "rotation": (0, 0, 0)},
])
def test_overconstrained_or_incomplete_placements_fail_atomically(kwargs):
    scene = Scene()
    with pytest.raises(ValueError):
        scene.add(PipeComponent(pipe()), **kwargs)
    assert not scene.nodes


def test_duplicate_names_invalid_meshes_fields_and_transforms():
    scene = Scene()
    scene.add(PipeComponent(pipe()), name="feed")
    with pytest.raises(ValueError, match="Duplicate"):
        scene.add(PipeComponent(pipe()), name="feed")
    with pytest.raises(ValueError, match="rigid"):
        Transform(matrix=np.diag([2, 1, 1, 1]))
    with pytest.raises(ValueError):
        MeshData(np.zeros((3, 3)), np.array([[0, 1, 3]]))
    root = PipeComponent(pipe()).build()
    root.fields["bad"] = FieldData([1], "Pa", "Pressure")
    with pytest.raises(ValueError, match="length"):
        Scene().add(SimpleNamespace(build=lambda: root))


def test_pipe_uses_only_supplied_route_and_station_fields():
    model = pipe()
    original = model.centerline.copy()
    root = PipeComponent(model, {"p_static": FieldData([3e6, 2e6, 1e6], "Pa", "Static pressure")}).build()
    rings = root.mesh.points.reshape(3, 16, 3)
    np.testing.assert_allclose(rings.mean(axis=1), original, atol=1e-12)
    np.testing.assert_allclose(np.linalg.norm(rings-original[:, None, :], axis=2), 0.01, atol=1e-12)
    np.testing.assert_array_equal(model.centerline, original)
    np.testing.assert_array_equal(root.fields["p_static"].values, np.repeat([3e6, 2e6, 1e6], 16))


def bent_pipe(radius=0.05, diameter=0.01):
    """A solver pipe whose single corner is one lumped bend between tangents."""
    stations, connections = ducts.generate_path_stations([
        ducts.PathNode((0, 0, 0)),
        ducts.PathNode((1, 0, 0), bend=ducts.RoundedBend(radius=radius)),
        ducts.PathNode((1, 1, 0)),
    ])
    return ducts.Pipe(stations, connections, diameter=diameter)


class ConstantTransport:
    def get_rho(self, T, p):
        return 1000.0

    def get_mu(self, T, p):
        return 1.0e-3

    def get_cp(self, T, p):
        return 4180.0


def solved_bent_pipe():
    from pyskyfire.thermo import PHASE_UNKNOWN, FluidState
    pipe = bent_pipe()
    state = FluidState(T=300., p=2e6, h=np.nan, quality=np.nan, rho=np.nan,
                       phase=PHASE_UNKNOWN, T_sat=np.nan)
    return pipe, pipe.solve(inlet_state=state, mass_flow=0.1, transport=ConstantTransport())


def test_solver_pipe_is_drawn_through_its_own_refined_bend_geometry():
    pipe = bent_pipe()
    bend = pipe.connections[1]
    root = PipeComponent(pipe, n_theta=8, segments=6).build()
    rings = root.mesh.points.reshape(-1, 8, 3)
    centers = rings.mean(axis=1)

    # The corner is swept as an arc, not as the solver centerline's chord.
    assert len(centers) == len(pipe.display_path(segments=6).points) > len(pipe.stations)
    arc = centers[2:8]
    np.testing.assert_allclose(
        np.linalg.norm(arc-np.asarray(bend.center), axis=1), bend.radius, atol=1e-9)
    np.testing.assert_allclose(centers[0], pipe.centerline[0], atol=1e-12)
    np.testing.assert_allclose(centers[-1], pipe.centerline[-1], atol=1e-12)
    np.testing.assert_allclose(
        np.linalg.norm(rings-centers[:, None, :], axis=2), 0.005, atol=1e-12)
    assert root.anchors["inlet"].position.tolist() == [0, 0, 0]


def test_pipe_result_becomes_station_point_fields_and_connection_cell_fields():
    pipe, result = solved_bent_pipe()
    root = PipeComponent(pipe, result, n_theta=8, segments=4).build()
    path = pipe.display_path(segments=4)

    assert set(duct_fields(result)) == set(root.fields)
    for name in ("p_stagnation", "p_static", "T_stagnation", "T_static", "density", "velocity"):
        assert root.fields[name].association == "point"
        assert len(root.fields[name].values) == len(root.mesh.points)
    for name in ("pressure_loss", "reynolds", "darcy_factor", "loss_coefficient"):
        assert root.fields[name].association == "cell"
        assert len(root.fields[name].values) == len(root.mesh.triangles)

    # Station values land on their own stations; the bend interpolates between.
    pressure = root.fields["p_static"].values.reshape(-1, 8)[:, 0]
    np.testing.assert_allclose(pressure[[0, 1, 5, 6]], result.p_static, rtol=1e-12)
    assert np.all(np.diff(pressure) <= 0)
    # A lumped bend defines no Darcy factor, and the NaN survives the mapping.
    darcy = root.fields["darcy_factor"].values.reshape(-1, 2*8)[:, 0]
    np.testing.assert_array_equal(np.isnan(darcy), np.isnan(result.darcy_factor[path.connection]))
    assert np.isnan(result.darcy_factor[1])


def test_pipe_fields_reject_wrong_association_lengths():
    pipe, result = solved_bent_pipe()
    with pytest.raises(ValueError, match="one value per centerline station"):
        PipeComponent(pipe, {"p": FieldData([1., 2.], "Pa", "P")}).build()
    with pytest.raises(ValueError, match="one value per pipe connection"):
        PipeComponent(pipe, {"p": FieldData([1., 2.], "Pa", "P", "cell")}).build()
    with pytest.raises(TypeError, match="DuctResult"):
        duct_fields(result.stations)


def test_solver_pipe_shares_a_color_field_with_a_chamber():
    chamber, layout, chamber_result = chamber_case()
    pipe, result = solved_bent_pipe()
    scene = Scene(color_by="p_static")
    scene.add(ThrustChamberComponent(chamber, chamber_result, layout=layout), name="chamber")
    duct = scene.add(PipeComponent(pipe, result), name="feed", position=(0, 0, 1))
    limits = scene_payload(scene)["style"]["limits"]
    assert limits[0] <= float(result.p_static.min()) and limits[1] >= float(result.p_static.max())
    assert duct.metadata["kind"] == "pipe"


def test_impeller_can_share_scene_with_channels_and_pipe():
    pytest.importorskip("pyturbo")
    from pyskyfire.pump import Impeller
    imp = Impeller(Q=0.01, H=100, n=10000, shaft_diameter=0.02,
                   inlet_hub_diameter=0.03, blade_count=6, beta2_deg=25,
                   blade_thickness=0.0015, hydraulic_efficiency=0.85)
    imp.build_geometry(span_points=5, chord_points=33)
    scene = Scene()
    pump = scene.add(ImpellerComponent(imp, n_theta=12), position=(1, 2, 3))
    assert {"hub", "shroud", "blade_0"} <= {c.name for c in pump.children}
    np.testing.assert_allclose(pump.anchor("shaft").position, [1+imp.geometry.hub_curve[0, 0], 2, 3])
    chamber, layout, result = chamber_case()
    scene.add(ThrustChamberComponent(chamber, result, layout=layout))
    scene.add(PipeComponent(pipe()))
    payload = scene_payload(scene)
    assert len(payload["nodes"]) == 3
    assert len(payload["meshes"]) == 5  # hub, shroud, shared blade, chamber, pipe


def test_html_contains_fields_transform_hierarchy_and_safe_text(tmp_path):
    chamber, layout, result = chamber_case()
    scene = Scene(title="</script><evil>", color_by="p_static", clim=(0, 2e6))
    scene.add(ThrustChamberComponent(chamber, result, layout=layout), position=(1, 2, 3))
    output = scene.save_html(tmp_path/"nested"/"engine.html")
    document = output.read_text()
    payload = json.loads(document.split('<script id="psf-scene" type="application/json">')[1].split('</script>')[0])
    assert payload["style"]["limits"] == [0, 2e6]
    assert payload["nodes"][0]["transform"][12:15] == [1, 2, 3]
    assert "p_static" in payload["nodes"][0]["children"][0]["children"][0]["fields"]
    assert "</script><evil>" not in document
    assert base64.b64decode(scene.data_url.split(",", 1)[1]).decode() == document


def test_desktop_keeps_scalar_arrays_ids_and_world_coordinates():
    scene = Scene(color_by="p_static")
    node = scene.add(PipeComponent(pipe(), {"p_static": FieldData([1e6, 2e6, 3e6], "Pa", "Static pressure")}), position=(2, 0, 0))
    pl = scene.to_plotter(off_screen=True)
    try:
        actor = pl.actors[node.id]
        mesh = actor.mapper.dataset
        np.testing.assert_allclose(mesh.points, node.world_transform.apply(node.mesh.points))
        np.testing.assert_allclose(mesh.point_data["p_static"], node.fields["p_static"].values)
        assert mesh.field_data["scene_node_id"][0] == node.id
        assert actor.mapper.scalar_range == (1e6, 3e6)
        image = pl.screenshot(return_img=True)
        assert image.shape[2] == 3 and image.min() < 200
    finally:
        pl.close()
        pl.deep_clean()


def test_shared_field_range_missing_fields_and_cell_values():
    scene = Scene(color_by="p_static")
    first = scene.add(PipeComponent(pipe(), {"p_static": FieldData([1., 2., 3.], "Pa", "Pressure")}), name="first")
    second = scene.add(PipeComponent(pipe()), name="second")
    second.fields["p_static"] = FieldData(np.full(len(second.mesh.triangles), 8), "Pa", "Pressure", "cell")
    assert scene_payload(scene)["style"]["limits"] == [1, 8]
    second.visible = False
    assert scene_payload(scene)["style"]["limits"] == [1, 3]
    first.visible = False
    with pytest.raises(ValueError, match="No visible"):
        scene_payload(scene)


def test_legacy_wrapper_warns_and_uses_scene(monkeypatch):
    from pyskyfire.viz.engine_viz import make_engine_3d
    chamber, layout, _ = chamber_case()
    plotter = Mock()
    monkeypatch.setattr(Scene, "to_plotter", lambda *a, **kw: plotter)
    with pytest.warns(DeprecationWarning):
        viewer = make_engine_3d(chamber, layout=layout, show=True)
    plotter.show.assert_called_once_with(auto_close=False)
    assert "psf-scene" in base64.b64decode(viewer.data_url.split(",", 1)[1]).decode()
    viewer.close()


def test_axes_triad_is_opt_in_for_both_backends():
    assert scene_payload(Scene())["show_axes"] is False
    scene = Scene(show_axes=True)
    scene.add(PipeComponent(pipe()))
    assert scene_payload(scene)["show_axes"] is True
    assert "data.show_axes" in files("pyskyfire.viz").joinpath("assets/scene_viewer.js").read_text()
    pl = scene.to_plotter(off_screen=True)
    try:
        assert pl.renderer.axes_widget.GetEnabled()
    finally:
        pl.close()
        pl.deep_clean()


def test_payload_offers_every_field_and_the_placement_arguments():
    scene = Scene(color_by="p_static", clim=(0, 5e6))
    fields = {"p_static": FieldData([1e6, 2e6, 3e6], "Pa", "Static pressure"),
              "T_static": FieldData([120., 200., 300.], "K", "Static temperature")}
    scene.add(PipeComponent(pipe(), fields), name="duct", position=(1, 2, 3), rotation=(0, 0.5, 0))
    payload = scene_payload(scene)
    # clim pins the selected field only; the rest range over their own values.
    assert payload["styles"]["p_static"]["limits"] == [0, 5e6]
    assert payload["styles"]["T_static"]["limits"] == [120, 300]
    assert payload["styles"]["T_static"]["units"] == "K"
    placement = payload["nodes"][0]["placement"]
    np.testing.assert_allclose(placement["position"], [1, 2, 3])
    np.testing.assert_allclose(placement["rotation"], [0, 0.5, 0], atol=1e-12)


def test_placement_arguments_reproduce_an_anchored_transform():
    scene = Scene()
    first = scene.add(PipeComponent(pipe()), name="first")
    second = scene.add(PipeComponent(pipe()), name="second", anchor="inlet",
                       attach_to=first.anchor("outlet"))
    placement = scene_payload(scene)["nodes"][1]["placement"]
    rebuilt = Transform(placement["position"], placement["rotation"])
    np.testing.assert_allclose(rebuilt.matrix, second.transform.matrix, atol=1e-12)


def test_a_field_with_inconsistent_units_is_not_offered():
    scene = Scene()
    scene.add(PipeComponent(pipe(), {"q": FieldData([1., 2., 3.], "Pa", "Quantity")}), name="first")
    scene.add(PipeComponent(pipe(), {"q": FieldData([1., 2., 3.], "K", "Quantity")}), name="second")
    assert scene_payload(scene)["styles"] == {}
    scene.color_by = "q"
    with pytest.raises(ValueError, match="consistent SI units"):
        scene_payload(scene)


def volute_case(tilt=0.0):
    section = ducts.EllipticalSlottedSection(0.012, 0.003, throat_length=0.002)
    inlet = ducts.VoluteInlet(radius=0.060, C_m=6.0, C_u=45.0)
    return ducts.Volute(section, inlet, ducts.ConstantAngularMomentum(),
                        rotation="ccw", tilt=tilt, n_stations=25)


def _mesh_area(mesh):
    corners = mesh.points[mesh.triangles]
    return 0.5*np.linalg.norm(np.cross(corners[:, 1]-corners[:, 0],
                                       corners[:, 2]-corners[:, 0]), axis=1).sum()


def test_volute_component_separates_its_openings_from_its_walls():
    """Walls, slot mouth and flange come back as three nodes, not one mesh."""
    volute = volute_case()
    node = VoluteComponent(volute, n_points=8).build()
    assert node.metadata["kind"] == "volute"
    assert set(node.anchors) == {"inlet", "outlet"}
    assert node.mesh is None
    assert [child.name for child in node.children] == ["passage", "inlet", "outlet"]

    passage, inlet, outlet = node.children
    # The openings carry no fields, which is what keeps them their own colour
    # whichever field the scene is coloured by.
    assert not inlet.fields and not outlet.fields
    assert inlet.color == outlet.color == "red"
    assert passage.color != "red"

    # The slot mouth is the strip of inlet-radius wall the flow enters through:
    # the slot width times the swept circumference, which the mesh resolves as
    # an inscribed polygon on the volute's own station count.
    segments = len(volute.geometry.theta) - 1
    expected = (volute.section.inlet_width
                * 2*segments*volute.inlet.radius*np.sin(np.pi/segments))
    assert _mesh_area(inlet.mesh) == pytest.approx(expected, rel=1e-9)


def test_volute_opening_colour_is_configurable():
    node = VoluteComponent(volute_case(), n_points=8, opening_color="#00ff00").build()
    assert [child.color for child in node.children[1:]] == ["#00ff00", "#00ff00"]


def test_volute_fields_are_repeated_around_each_section_ring():
    volute = volute_case()

    class Transport:
        def get_rho(self, T, p):
            return 1000.0

        def get_mu(self, T, p):
            return 1.0e-3

    state = ducts.pipe._state(300.0, 8.0e6, 1000.0)
    result = volute.solve(inlet_state=state, mass_flow=27.14,
                          transport=Transport(), exit_diameter=0.05)
    fields = volute_fields(result)
    assert fields["area"].values == pytest.approx(volute.geometry.area)
    assert fields["collected"].values[-1] == pytest.approx(1.0)
    # Section geometry that says nothing about the flow is off the dropdown.
    assert "passage_radius" not in fields and "centroid_radius" not in fields
    assert {"p_static", "p_stagnation", "velocity"} <= set(fields)

    passage = VoluteComponent(volute, result, n_points=8).build().children[0]
    for field in passage.fields.values():
        assert len(field.values) == len(passage.mesh.points)
    scene = Scene()
    scene.add(VoluteComponent(volute, result, n_points=8), name="volute")
    assert scene_payload(scene)["nodes"]


def test_volute_component_rejects_mismatched_field_lengths():
    volute = volute_case()
    bad = {"junk": FieldData(np.zeros(3), "1", "Junk")}
    with pytest.raises(ValueError, match="one value per station"):
        VoluteComponent(volute, bad, n_points=8).build()


class NamedTransport(ConstantTransport):
    def __init__(self, fluid):
        self.fluid = fluid


def test_tags_label_added_nodes_and_travel_to_the_viewer():
    scene = Scene()
    feed = scene.add(PipeComponent(pipe()), name="feed", tags="fuel path")
    ret = scene.add(PipeComponent(pipe()), name="return", tags=["fuel path", "hot"],
                    position=(0, 0, 1))
    ret.children.append(SceneNode("liner", tags={"hot"}))
    assert feed.tags == {"fuel path"}
    payload = scene_payload(scene)
    assert [node["tags"] for node in payload["nodes"]] == [["fuel path"], ["fuel path", "hot"]]
    with pytest.raises(ValueError, match="tags"):
        scene.add(PipeComponent(pipe()), name="bad", tags=["", "x"])
    assert [node.name for node in scene.nodes] == ["feed", "return"]


def test_solved_components_record_the_fluid_they_carry():
    from pyskyfire.thermo import PHASE_UNKNOWN, FluidState, Fluid
    state = FluidState(T=300., p=2e6, h=np.nan, quality=np.nan, rho=np.nan,
                       phase=PHASE_UNKNOWN, T_sat=np.nan)
    duct = bent_pipe()
    # thermo.Fluid spells a pure fluid as a one-part HEOS mixture; an alias
    # such as H2 names the same fluid, so both land in one group.
    hydrogen = NamedTransport(Fluid("pure", ["Hydrogen"], [1.0]).coolprop_string())
    solved = duct.solve(inlet_state=state, mass_flow=0.1, transport=hydrogen)
    assert solved.fluid == "HEOS::Hydrogen[1.000]"
    assert PipeComponent(duct, solved).build().metadata["fluid"] == "Hydrogen"
    alias = duct.solve(inlet_state=state, mass_flow=0.1, transport=NamedTransport("H2"))
    assert PipeComponent(duct, alias).build().metadata["fluid"] == "Hydrogen"
    # Constant properties, or fields given by hand, name no fluid.
    _, unnamed = solved_bent_pipe()
    assert unnamed.fluid is None
    assert "fluid" not in PipeComponent(duct, unnamed).build().metadata
    assert "fluid" not in PipeComponent(duct, duct_fields(solved)).build().metadata

    chamber, layout, result = chamber_case()
    chamber.cooling_circuits[0].coolant_transport = NamedTransport("Methane")
    root = ThrustChamberComponent(chamber, result, layout=layout).build()
    assert root.children[0].metadata["fluid"] == "Methane"
    assert "fluid" not in ThrustChamberComponent(chamber, layout=layout).build().children[0].metadata


def test_mixtures_are_labelled_by_mass():
    from pyskyfire.thermo import Fluid
    from pyskyfire.viz.components import _fluid_label
    # thermo.Fluid hands CoolProp mole fractions; the label reads by mass again.
    ethalox = Fluid("fuel", ["Ethanol", "Water"], [0.8, 0.2]).coolprop_string()
    assert ethalox == "HEOS::Ethanol[0.610]&Water[0.390]"
    assert _fluid_label(ethalox) == "Ethanol 80 wt% & Water 20 wt%"
    # Unknown components keep the spelling they came with.
    assert _fluid_label("HEOS::Foo[0.5]&Water[0.5]") == "Foo[0.5]&Water[0.5]"


def test_parts_take_the_default_colour_unless_given_one():
    from pyskyfire.viz.scene import DEFAULT_COLOR
    chamber, layout, _ = chamber_case()
    channel = ThrustChamberComponent(chamber, layout=layout).build().children[0].children[0]
    assert channel.color == DEFAULT_COLOR
    painted = ThrustChamberComponent(chamber, layout=layout, circuit_colors=["red"]).build()
    assert painted.children[0].children[0].color == "red"


def test_repeated_channels_share_their_field_values_in_the_page():
    chamber, layout, result = chamber_case()
    scene = Scene(color_by="p_static")
    scene.add(ThrustChamberComponent(chamber, result, layout=layout), name="chamber")
    payload = scene_payload(scene)
    channels = payload["nodes"][0]["children"][0]["children"]
    assert len(channels) == 6
    # Six channels from one seed carry one copy of each field's values.
    refs = {channel["fields"]["p_static"]["values"] for channel in channels}
    assert len(refs) == 1 and len(payload["values"]) == len(channels[0]["fields"])
    values = payload["values"][refs.pop()]
    assert len(values) == len(scene.nodes[0].children[0].children[0].mesh.points)
