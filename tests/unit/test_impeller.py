"""Sizing, pyTurbo convention mapping, and physical blade/scene contracts."""
import subprocess
import sys

import numpy as np
import pytest

from pyskyfire.pump import Impeller, MeridionalContours, SplitterSettings
from pyskyfire.pump.geometry import cylindrical_to_cartesian


def design(**changes):
    args = dict(Q=0.01, H=100., n=10000., shaft_diameter=0.02,
                inlet_hub_diameter=0.03, blade_count=6, beta2_deg=25.,
                blade_thickness=0.0015)
    args.update(changes)
    return Impeller(**args)


def test_eye_continuity_uses_fixed_hub_and_internal_flow():
    imp = design(volumetric_efficiency=0.95)
    d = imp.dimensions
    area = np.pi/4*(d.d1**2-d.inlet_hub_diameter**2)
    assert 0.20*d.u1*area == pytest.approx(0.01/0.95)
    assert imp.geometry is None
    assert d.shaft_diameter == 0.02
    assert d.inlet_hub_diameter == 0.03


def test_inlet_triangles_incidence_prewhirl_and_override():
    imp = design(inlet_prewhirl=lambda span_fraction: 2*span_fraction,
                 incidence_deg=2.)
    for triangle in imp.velocity_triangles.inlet:
        assert triangle.u == pytest.approx(2*np.pi*imp.n/60*triangle.radius)
        assert triangle.w == pytest.approx(np.hypot(triangle.cm, triangle.u-triangle.cu))
        flow = np.degrees(np.arctan2(triangle.cm, triangle.u-triangle.cu))
        assert triangle.beta_metal_deg == pytest.approx(flow+2)
    override = design(beta1_deg=lambda span_fraction: 20+10*span_fraction)
    assert [s.beta1_deg for s in override.blade_definition.sections] == [20,25,30]


def test_meridional_passage_stays_open_and_preserves_end_dimensions():
    imp = design()
    d = imp.dimensions
    mer = imp.blade_definition.meridional
    np.testing.assert_allclose(mer.hub[0], [0, d.inlet_hub_diameter/2])
    np.testing.assert_allclose(mer.shroud[0], [0, d.d1/2])
    assert mer.hub[-1,1] == pytest.approx(d.d2/2)
    assert mer.shroud[-1,1] == pytest.approx(d.d2/2)
    assert mer.hub[-1,0]-mer.shroud[-1,0] == pytest.approx(d.b2)
    cuts = mer.shroud-mer.hub
    assert np.all(np.linalg.norm(cuts, axis=1) > 0)
    # Tangents at inlet are axial, at outlet radial.
    for curve in (mer.hub, mer.shroud):
        assert abs(np.diff(curve[:2],axis=0)[0,1]/np.diff(curve[:2],axis=0)[0,0]) < 0.02
        assert abs(np.diff(curve[-2:],axis=0)[0,0]/np.diff(curve[-2:],axis=0)[0,1]) < 0.02
    custom = design(meridional=mer)
    np.testing.assert_array_equal(custom.blade_definition.meridional.hub, mer.hub)
    bad = mer.shroud.copy()
    bad[-1] = mer.hub[-1]
    with pytest.raises(ValueError, match="separated"):
        MeridionalContours(mer.hub, bad)


def test_supplied_efficiencies_and_pressure_rise_constructor():
    imp = Impeller.from_pressure_rise(mdot=10, rho=1000, dp=980665, n=10000,
        shaft_diameter=.02, inlet_hub_diameter=.03, blade_count=6, beta2_deg=25,
        blade_thickness=.0015, volumetric_efficiency=.97,
        hydraulic_efficiency=.85, mechanical_efficiency=.98)
    assert imp.Q == pytest.approx(.01)
    assert imp.H == pytest.approx(100)
    assert imp.overall_efficiency == pytest.approx(.97*.85*.98)
    assert imp.shaft_power == pytest.approx(980665*.01/(.97*.85*.98))
    assert design().shaft_power is None
    assert design().velocity_triangles.hydraulic_head is None


@pytest.mark.parametrize("changes, message", [
    ({"Q":np.nan}, "Q"), ({"n":0}, "n"),
    ({"blade_count":6.5}, "blade_count"),
    ({"inlet_hub_diameter":.01}, "shaft_diameter"),
    ({"d1":.5}, "requires"), ({"d1":.02}, "requires"),
    ({"volumetric_efficiency":1.1}, "volumetric_efficiency"),
    ({"hydraulic_efficiency":np.nan}, "hydraulic_efficiency"),
    ({"blade_thickness":0}, "blade_thickness"),
    ({"blade_thickness":lambda unknown: 1.}, "unknown argument"),
    ({"blade_thickness":lambda: [1,2]}, "finite scalar"),
    ({"beta1_deg":0}, "metal angles"),
    ({"wrap_deg":-30}, "wrap_deg"),
    ({"tip_clearance":1}, "tip_clearance"),
])
def test_invalid_design_inputs_are_rejected(changes, message):
    with pytest.raises(ValueError, match=message):
        design(**changes)


def test_sizing_does_not_import_optional_backend_or_plotting():
    script = """
import sys
from pyskyfire.pump import Impeller
i = Impeller(Q=.01,H=100,n=10000,shaft_diameter=.02,inlet_hub_diameter=.03,
             blade_count=6,beta2_deg=25,blade_thickness=.0015)
assert i.geometry is None
assert not {'pyturbo', 'matplotlib.pyplot', 'plotly', 'pyvista'} & sys.modules.keys()
"""
    subprocess.run([sys.executable,"-c",script], check=True, capture_output=True, text=True)


def test_missing_backend_has_install_hint(monkeypatch):
    from pyskyfire.pump import _pyturbo
    def absent(name):
        raise _pyturbo.PackageNotFoundError(name)
    monkeypatch.setattr(_pyturbo, "version", absent)
    with pytest.raises(ImportError, match=r"pyskyfire\[impeller\]"):
        design().build_geometry()


@pytest.fixture(scope="module")
def built():
    pytest.importorskip("pyturbo")
    imp = design(wrap_deg=160., beta1_deg=25., beta2_deg=30.)
    imp.build_geometry(span_points=9, chord_points=97)
    return imp


def test_actual_pyturbo_wrap_sweep_and_thickness(built):
    g = built.geometry
    blade = g.main_blade
    assert blade.suction.shape == (9,97,3)
    for side in (blade.suction, blade.pressure):
        np.testing.assert_allclose(np.degrees(side[:,-1,2]-side[:,0,2]), -160, atol=1e-8)
        assert not side.flags.writeable
    assert g.thickness_max_relative_error < .03
    # Independent Cartesian measurement, including interpolated span sections.
    ss = cylindrical_to_cartesian(blade.suction)
    ps = cylindrical_to_cartesian(blade.pressure)
    gap = np.linalg.norm(ss-ps, axis=-1)
    np.testing.assert_allclose(gap[:,20:75], .0015, rtol=.08)
    np.testing.assert_allclose(ss[:,[0,-1]], ps[:,[0,-1]], atol=1e-10)
    # Recover metal angles from the physical mean-surface endpoint tangents.
    # Finite sampling and the fitted rounded edges introduce a few degrees of
    # error, but a wrong upstream angle convention must fail this check.
    mean = (blade.suction+blade.pressure)/2
    dm = np.linalg.norm(np.diff(mean[...,:2], axis=1), axis=-1)
    r_dtheta = .5*(mean[:,:-1,1]+mean[:,1:,1])*np.diff(mean[...,2], axis=1)
    beta = np.degrees(np.arctan2(dm, -r_dtheta))
    np.testing.assert_allclose(beta[:,0], 25., atol=3.)
    np.testing.assert_allclose(beta[:,-1], 30., atol=3.)


def test_callable_thickness_clearance_splitters_and_repeat_build():
    pytest.importorskip("pyturbo")
    from pyturbo.aero import Centrif
    original_patterns = list(Centrif.patterns)
    imp = design(blade_thickness=lambda chord_fraction,span_fraction: .0015*(1-.3*span_fraction)*(1+.2*chord_fraction),
                 splitter=SplitterSettings(start=.5,pitch_fraction=.4), tip_clearance=.0002)
    a = imp.build_geometry(span_points=5,chord_points=65)
    b = imp.build_geometry(span_points=5,chord_points=65)
    np.testing.assert_array_equal(a.main_blade.suction,b.main_blade.suction)
    assert Centrif.patterns == original_patterns
    assert a.splitter is not None
    assert a.thickness_max_relative_error < .03
    # Generated tip points lie on the trimmed passage; original shroud is kept.
    mer = imp.blade_definition.meridional
    tip = mer.with_clearance(.0002)
    np.testing.assert_allclose(np.linalg.norm(tip.shroud-mer.shroud,axis=1),.0002)
    from scipy.spatial import cKDTree
    samples = tip.sample(np.linspace(0,1,20000),1.)
    distances = cKDTree(samples).query(a.main_blade.suction[-1,:,:2])[0]
    assert distances.max() < 1e-5


def test_splitter_cannot_be_placed_inside_a_main_blade():
    pytest.importorskip("pyturbo")
    imp = design(splitter=SplitterSettings(pitch_fraction=.001))
    with pytest.raises(ValueError, match="Splitter overlaps"):
        imp.build_geometry(span_points=3, chord_points=65)
    assert imp.geometry is None


def test_scene_reuses_closed_blades_and_export_preserves_surfaces(built,tmp_path):
    from pyskyfire.viz import Scene, ImpellerComponent
    from pyskyfire.viz.backends.threejs import scene_payload
    from pyskyfire.viz.impeller_viz import make_blade_meshes, export_impeller_surfaces
    scene = Scene()
    node = scene.add(ImpellerComponent(built), position=(1,2,3))
    blades = [child for child in node.children if child.name.startswith("blade_")]
    assert len(blades) == 6
    assert all(b.mesh is blades[0].mesh for b in blades)
    assert len(scene_payload(scene)["meshes"]) == 3
    mesh = make_blade_meshes(built)[0]
    assert mesh.is_manifold and mesh.volume > 0
    faces = blades[0].mesh.triangles
    # Each shared edge has two opposite directed incidences.
    edges = np.concatenate((faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]))
    _,counts = np.unique(np.sort(edges,axis=1),axis=0,return_counts=True)
    assert np.all(counts == 2)
    assert len(np.unique(edges,axis=0)) == len(edges)
    archive = built.geometry.export_npz(tmp_path/"geometry.npz")
    with np.load(archive,allow_pickle=False) as saved:
        np.testing.assert_array_equal(saved["suction"],built.geometry.main_blade.suction)
    assert export_impeller_surfaces(built,tmp_path/"impeller.stl").stat().st_size > 1000
    assert scene.save_html(tmp_path/"scene.html").stat().st_size > 1000


def test_visualizer_requires_explicit_build():
    from pyskyfire.viz import ImpellerComponent
    with pytest.raises(ValueError,match="build_geometry"):
        ImpellerComponent(design()).build()


def test_sections_stack_on_a_radial_trailing_edge():
    pytest.importorskip("pyturbo")
    imp = design(splitter=SplitterSettings())
    geometry = imp.build_geometry(span_points=5, chord_points=65)
    wraps = [s.wrap_deg for s in imp.blade_definition.sections]
    # The angle law wraps the hub section much further than the tip section.
    assert wraps[0]-wraps[-1] > 30
    for blade in (geometry.main_blade, geometry.splitter):
        for side in (blade.suction, blade.pressure):
            # Trailing edge is one axial line: constant theta from hub to tip.
            np.testing.assert_allclose(side[:,-1,2], 0., atol=1e-9)
    # The leading edge carries the spanwise wrap difference as sweep, and each
    # section still turns through the wrap resolved for its own span.
    for side in (geometry.main_blade.suction, geometry.main_blade.pressure):
        np.testing.assert_allclose(np.degrees(side[[0,2,-1],0,2]), wraps, atol=1e-6)


def test_slanted_leading_edge_repairs_the_passage_and_evens_the_wrap():
    from pyskyfire.pump import BladeSettings
    plain = design(beta2_deg=35.)
    slanted = design(beta2_deg=35., blade_settings=BladeSettings(start=(.34, .05)))
    assert slanted.blade_definition.settings.start == (.34, .05)
    assert BladeSettings().start == (.05, .05)
    straight, slant = (imp.blade_definition.meridional for imp in (plain, slanted))
    # The hub keeps its stations; only the pairing of shroud stations moves.
    np.testing.assert_array_equal(slant.hub, straight.hub)
    from scipy.interpolate import PchipInterpolator
    from scipy.spatial import cKDTree
    dense = PchipInterpolator(straight.stations, straight.shroud)(np.linspace(0, 1, 20000))
    assert cKDTree(dense).query(slant.shroud)[0].max() < 1e-6
    # The leading edge now runs from a hub point well downstream to the eye.
    hub_le, tip_le = (slant.sample(.34, s) for s in (0., 1.))
    np.testing.assert_allclose(tip_le, straight.sample(.05, 1.), atol=5e-6)
    assert hub_le[0] > 5*tip_le[0] and hub_le[1] > straight.sample(.05, 0.)[1]
    # Hub and tip sections then wrap by comparable amounts.
    spread = [s[0].wrap_deg-s[-1].wrap_deg
              for s in (imp.blade_definition.sections for imp in (plain, slanted))]
    assert spread[0] > 70 and abs(spread[1]) < 20
    # A blade starting in the inlet plane still resolves its angles there.
    assert design(blade_settings=BladeSettings(start=0.)).velocity_triangles.inlet[0].cm > 0
    # The inlet triangle is taken on the leading edge, not at the inlet plane.
    for triangle in slanted.velocity_triangles.inlet:
        assert triangle.radius == pytest.approx(slant.sample(.34, triangle.span_fraction)[1])
        assert triangle.u == pytest.approx(2*np.pi*slanted.n/60*triangle.radius)


def test_slanted_leading_edge_reaches_pyturbo():
    pytest.importorskip("pyturbo")
    from pyskyfire.pump import BladeSettings
    imp = design(beta2_deg=35., blade_settings=BladeSettings(start=(.34, .05)))
    geometry = imp.build_geometry(span_points=5, chord_points=65)
    edge = geometry.main_blade.suction[:,0,:2]
    hub_le, tip_le = (imp.blade_definition.meridional.sample(.34, s) for s in (0., 1.))
    np.testing.assert_allclose(edge[[0,-1]], [hub_le, tip_le], atol=1e-5)
    # Sections still stack on the trailing edge; the leading edge sweep is now
    # small because hub and tip wrap by comparable amounts.
    np.testing.assert_allclose(geometry.main_blade.suction[:,-1,2], 0., atol=1e-9)
    assert abs(np.degrees(geometry.main_blade.suction[0,0,2]-geometry.main_blade.suction[-1,0,2])) < 20


def test_leading_edge_area_does_not_depend_on_the_pairing():
    from pyskyfire.pump import BladeSettings
    # Continuity through the leading edge is a property of the passage, so
    # slanting it must not move the flow area at a fixed station.
    cm = [design(blade_settings=BladeSettings(start=start)).velocity_triangles.inlet[0].cm
          for start in (.05, (.05, .02), (.25, .05), (.34, .05))]
    assert max(cm)/min(cm) < 1.05


@pytest.mark.parametrize("start, message", [
    ((.05, .34), "no earlier on the hub"),
    ((.05, .05, .05), "one station or a"),
    ((1.2, .05), "Blade extent"),
])
def test_invalid_blade_start_is_rejected(start, message):
    from pyskyfire.pump import BladeSettings
    with pytest.raises(ValueError, match=message):
        BladeSettings(start=start)


def test_leading_edge_slant_must_fit_the_passage():
    from pyskyfire.pump import BladeSettings
    with pytest.raises(ValueError, match="does not fit this passage"):
        design(blade_settings=BladeSettings(start=(.95, .01)))
