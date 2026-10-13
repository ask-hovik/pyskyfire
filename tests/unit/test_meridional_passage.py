"""The shared meridional passage: section properties, wall motion and clearances."""

import numpy as np
import pytest

from pyskyfire.common import Clearance, MeridionalPassage, WallMotion
from pyskyfire.pump import MeridionalContours


def annulus(r_hub=.05, r_shroud=.08, length=.1, n=5, **kwargs):
    x = np.linspace(0, length, n)
    return MeridionalPassage(np.column_stack((x, np.full(n, r_hub))),
                             np.column_stack((x, np.full(n, r_shroud))), **kwargs)


def test_axial_annulus_area_and_hydraulic_diameter():
    p = annulus()
    np.testing.assert_allclose(p.area(), np.pi*(.08**2 - .05**2))
    np.testing.assert_allclose(p.hydraulic_diameter(), .06)
    assert p.area(.37) == pytest.approx(np.pi*(.08**2 - .05**2))


def test_radial_passage_area_is_a_cylinder():
    r = np.linspace(.1, .2, 6)
    p = MeridionalPassage(np.column_stack((np.full(6, .01), r)), np.column_stack((np.zeros(6), r)))
    np.testing.assert_allclose(p.area(), 2*np.pi*r*.01)
    np.testing.assert_allclose(p.hydraulic_diameter(), .02)


def test_two_point_walls_are_straight_and_sample_by_span():
    p = MeridionalPassage([[0, .05], [.1, .04]], [[0, .08], [.1, .09]])
    np.testing.assert_allclose(p.sample(.5, .5), [.05, .065])
    np.testing.assert_allclose(p.stations, [0, 1])
    np.testing.assert_allclose(p.wall_tangents(0.)[0], np.array([.1, -.01])/np.hypot(.1, .01))


def test_a_hub_on_the_axis_is_a_plain_pipe():
    p = annulus(r_hub=0.)
    np.testing.assert_allclose(p.area(), np.pi*.08**2)


def test_wall_motion_intervals():
    motion = WallMotion(hub=[(.6, .9), (.2, .5)])
    assert motion.hub == ((.2, .5), (.6, .9))
    np.testing.assert_array_equal(motion.rotating("hub", [.1, .2, .55, .7]), [False, True, False, True])
    assert not motion.rotating("shroud", .5)
    np.testing.assert_array_equal(annulus(motion=motion).rotating("hub"),
                                  [False, True, True, True, False])
    assert WallMotion.whole("hub") == WallMotion(hub=[(0, 1)])


@pytest.mark.parametrize("kwargs, message", [
    (dict(hub=[(.2, .5), (.4, .6)]), "overlap"),
    (dict(hub=[(.5, .5)]), "0 <= start < end <= 1"),
    (dict(shroud=[(-.1, .5)]), "0 <= start < end <= 1"),
])
def test_wall_motion_rejects_bad_intervals(kwargs, message):
    with pytest.raises(ValueError, match=message):
        WallMotion(**kwargs)


def test_clearances_are_recorded_and_checked_against_the_passage():
    tip = Clearance("R1", "shroud", .4, .8, .0005)
    p = annulus(clearances=[tip, Clearance("S1", "hub", .1, .3, .0004)])
    assert p.clearance("R1") == (tip,)
    with pytest.raises(ValueError, match="not smaller"):
        annulus(clearances=[Clearance("R1", "shroud", .4, .8, .03)])
    with pytest.raises(ValueError, match="wall must be"):
        Clearance("R1", "casing", .4, .8, .001)
    with pytest.raises(ValueError, match="positive finite gap"):
        Clearance("R1", "shroud", .4, .8, 0.)


def test_with_clearance_moves_one_wall_and_keeps_the_record():
    motion = WallMotion.whole("hub")
    tip = Clearance("R1", "shroud", 0., 1., .001)
    p = annulus(motion=motion, clearances=[tip])
    blade = p.with_clearance(.001)
    np.testing.assert_allclose(blade.shroud[:, 1], .079)
    np.testing.assert_array_equal(blade.hub, p.hub)
    assert blade.motion == motion and blade.clearances == (tip,)
    np.testing.assert_allclose(p.with_clearance(.002, wall="hub").hub[:, 1], .052)
    with pytest.raises(ValueError, match="smaller than the passage height"):
        p.with_clearance(.05)


def test_surface_revolves_the_wall_about_x():
    xyz = annulus().surface("shroud", n_theta=12)
    assert xyz.shape == (5, 12, 3)
    np.testing.assert_allclose(np.hypot(xyz[..., 1], xyz[..., 2]), .08)


@pytest.mark.parametrize("hub, shroud, message", [
    ([[0, .05]], [[0, .08]], "N >= 2"),
    ([[0, .05], [.1, .05]], [[0, .08], [.1, .08], [.2, .08]], "paired"),
    ([[0, .05], [.1, .05]], [[0, .04], [.1, .08]], "below hub"),
    ([[0, .05], [0, .05]], [[0, .08], [.1, .08]], "distinct"),
    ([[0, -.01], [.1, .05]], [[0, .08], [.1, .08]], "negative"),
    ([[0, .05], [.1, .08]], [[0, .08], [.1, .08]], "separated"),
])
def test_passage_rejects_bad_walls(hub, shroud, message):
    with pytest.raises(ValueError, match=message):
        MeridionalPassage(hub, shroud)


def test_impeller_contours_are_passages_and_keep_their_record_when_repaired():
    t = np.linspace(0, np.pi/2, 20)
    hub = np.column_stack((.03*np.sin(t), .015 + .06*(1 - np.cos(t))))
    shroud = np.column_stack((.02*np.sin(t), .03 + .045*(1 - np.cos(t))))
    motion = WallMotion.whole("hub")
    mer = MeridionalContours(hub, shroud, motion=motion)
    assert isinstance(mer, MeridionalPassage)
    slanted = mer.with_leading_edge(.3, .1)
    assert type(slanted) is MeridionalContours and slanted.motion == motion
    with pytest.raises(ValueError, match="nondecreasing"):
        MeridionalContours(hub[::-1], shroud[::-1])


def test_a_passage_pickled_before_wall_motion_existed_stands_still():
    p = annulus()
    state = dict(p.__dict__)
    del state["motion"], state["clearances"]
    old = object.__new__(MeridionalPassage)
    old.__dict__.update(state)
    assert old.motion == WallMotion() and old.clearances == ()
    assert not old.rotating("hub").any()
