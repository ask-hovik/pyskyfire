"""Resolved section geometry is built once per geometry, not once per query.

A single-station profile carries every field the N-station bundle does, so the
only way to get one is to build the bundle and interpolate it. That is correct
but expensive: assembling the bundle walks every station of every wall and
evaluates each design-input function a station at a time. The solver reaches a
station through ``section_profiles_at`` on every residual it forms, so without
the caches here the RL10 validation case spent roughly 80% of its runtime
rebuilding geometry that had not changed.

These tests pin the caching down from both sides: that a repeated query is
free, and that a geometry change is still seen.
"""

import numpy as np
import pytest

import pyskyfire as psf
from pyskyfire.regen.cross_section import SectionProfiles


class _Cylinder:
    """Constant-radius contour, so the section is predictable everywhere."""

    def __init__(self, radius=0.1, length=1.0):
        self.radius = float(radius)
        self.xs = np.array([0.0, length])
        self.x_t = 0.0
        self.r_t = float(radius)
        self.r_c = float(radius)
        self.r_e = float(radius)

    def r(self, x):
        return np.full_like(np.asarray(x, dtype=float), self.radius)

    def dr_dx(self, x):
        return np.zeros_like(np.asarray(x, dtype=float))

    def normal_angle(self, x):
        return np.zeros_like(np.asarray(x, dtype=float))


class _Material:
    @staticmethod
    def get_k(T):
        return 350.0


def _chamber(*, cross_section=None, wall_thickness=0.001):
    contour = _Cylinder()
    if cross_section is None:
        cross_section = psf.regen.CrossSectionSquared(sidewall_half_thickness=(5.0e-4) / 2)
    circuit = psf.regen.CoolingCircuit(
        name="test",
        contour=contour,
        cross_section=cross_section,
        span=[0.0, 1.0],
        placement=psf.regen.SurfacePlacement(n_channel_positions=60),
        channel_height=0.004,
        walls=[psf.regen.Wall(_Material(), wall_thickness)],
        coolant_transport=None,
        roughness=0.0,
    )
    chamber = psf.regen.ThrustChamber(
        contour=contour,
        cooling_circuits=[circuit],
        combustion_transport=None,
        n_nodes=20,
        compute_gas=False,
    )
    return chamber, circuit


def _count_calls(monkeypatch, owner, name):
    """Replace ``owner.name`` with a counting wrapper and return the counter."""
    original = getattr(owner, name)
    calls = []

    def counted(self, *args, **kwargs):
        calls.append(None)
        return original(self, *args, **kwargs)

    monkeypatch.setattr(owner, name, counted)
    return calls


# --- the bundle -------------------------------------------------------------


def test_the_profile_bundle_is_built_once_per_path(monkeypatch) -> None:
    _, circuit = _chamber()
    builds = _count_calls(monkeypatch, psf.regen.CoolingCircuit, "_build_prof")
    # Assembling the chamber has already populated the cache, so start clean.
    circuit._invalidate_section_cache()

    path = circuit.cross_section_base_point_path
    first = circuit._prof(path)
    for _ in range(5):
        assert circuit._prof(path) is first
    assert len(builds) == 1


def test_a_height_change_rebuilds_the_bundle() -> None:
    _, circuit = _chamber()
    path = circuit.cross_section_base_point_path
    before = circuit._prof(path)

    circuit.set_channel_height(np.full_like(circuit.channel_heights, 0.006))
    after = circuit._prof(path)

    assert after is not before
    np.testing.assert_allclose(after.h, 0.006)


def test_a_new_path_does_not_reuse_the_old_bundle() -> None:
    _, circuit = _chamber()
    bundle = circuit._prof(circuit.cross_section_base_point_path)
    other = circuit._prof(circuit.reference_path)
    assert other is not bundle


# --- single stations --------------------------------------------------------


def test_a_repeated_station_query_is_served_from_the_cache(monkeypatch) -> None:
    _, circuit = _chamber()
    builds = _count_calls(
        monkeypatch, psf.regen.CoolingCircuit, "_build_section_profiles_at"
    )

    first = circuit.section_profiles_at(0.5)
    for _ in range(20):
        assert circuit.section_profiles_at(0.5) is first
    assert len(builds) == 1

    circuit.section_profiles_at(0.6)
    assert len(builds) == 2


def test_cached_stations_do_not_outlive_their_geometry() -> None:
    _, circuit = _chamber()
    before = circuit.section_profiles_at(0.5)
    np.testing.assert_allclose(before.h, 0.004)

    circuit.set_channel_height(np.full_like(circuit.channel_heights, 0.006))
    after = circuit.section_profiles_at(0.5)

    assert after is not before
    np.testing.assert_allclose(after.h, 0.006)


def test_a_station_matches_an_uncached_rebuild() -> None:
    """The cache must not be the only thing that makes the fast path agree."""
    _, circuit = _chamber()
    cached = circuit.section_profiles_at(0.37)
    fresh = circuit._build_section_profiles_at(0.37)
    for name in ("h", "theta", "t_wall", "cross_section_base_point_path",
                 "projection", "azimuthal_pitch", "winding_rate"):
        expected, actual = getattr(fresh, name), getattr(cached, name)
        if expected is None:
            assert actual is None
        else:
            np.testing.assert_array_equal(actual, expected)


# --- derived geometry on a profile -----------------------------------------


def test_a_section_integrates_its_boundary_once_per_profile(monkeypatch) -> None:
    import pyskyfire.regen.cross_section as cs

    _, circuit = _chamber()
    prof = circuit.section_profiles_at(0.5)
    section = circuit.cross_section

    integrations = []
    original = cs.arc_boundary_metrics

    def counted(*args, **kwargs):
        integrations.append(None)
        return original(*args, **kwargs)

    monkeypatch.setattr(cs, "arc_boundary_metrics", counted)

    area = section.A_coolant(prof)
    for _ in range(10):
        np.testing.assert_array_equal(section.A_coolant(prof), area)
        section.P_coolant(prof)
        section.Dh_coolant(prof)
    assert len(integrations) == 1


def test_derived_geometry_is_not_shared_between_profiles() -> None:
    _, circuit = _chamber()
    section = circuit.cross_section
    thin = circuit.section_profiles_at(0.5)
    thin_area = section.A_coolant(thin)

    circuit.set_channel_height(np.full_like(circuit.channel_heights, 0.008))
    thick_area = section.A_coolant(circuit.section_profiles_at(0.5))

    assert thick_area[0] > thin_area[0]


def test_a_replaced_profile_starts_with_an_empty_memo() -> None:
    from dataclasses import replace

    _, circuit = _chamber()
    prof = circuit.section_profiles_at(0.5)
    circuit.cross_section.A_coolant(prof)
    assert prof._derived

    assert not replace(prof, h=prof.h * 1.5)._derived


def test_returned_arrays_are_not_views_into_the_memo() -> None:
    """Callers have always been handed arrays of their own; keep it that way."""
    _, circuit = _chamber()
    section = circuit.cross_section
    prof = circuit.section_profiles_at(0.5)

    area = section.A_coolant(prof)
    perimeter = section.P_coolant(prof)
    width = section.sidewall_half_thickness(prof)
    area[:] = -1.0
    perimeter[:] = -1.0
    width[:] = -1.0

    assert section.A_coolant(prof)[0] > 0.0
    assert section.P_coolant(prof)[0] > 0.0
    assert section.sidewall_half_thickness(prof)[0] > 0.0


# --- serialisation ----------------------------------------------------------


def test_pickling_a_circuit_leaves_the_caches_behind() -> None:
    # cloudpickle, because that is what `Results.save` uses and the design-input
    # evaluators are closures that plain pickle cannot take.
    import cloudpickle as pickle

    _, circuit = _chamber()
    circuit.section_profiles_at(0.5)
    assert circuit._section_cache is not None

    restored = pickle.loads(pickle.dumps(circuit))
    assert restored._section_cache is None
    np.testing.assert_allclose(
        restored.section_profiles_at(0.5).h,
        circuit.section_profiles_at(0.5).h,
    )


def test_a_pickled_profile_rebuilds_its_memo() -> None:
    import cloudpickle as pickle

    _, circuit = _chamber()
    prof = circuit.section_profiles_at(0.5)
    circuit.cross_section.A_coolant(prof)

    restored = pickle.loads(pickle.dumps(prof))
    assert restored._derived == {}
    np.testing.assert_allclose(restored.h, prof.h)


# --- the solver path --------------------------------------------------------


def test_the_rounded_section_reaches_the_same_resistance_either_way() -> None:
    """``R_coolant_per_len`` is the solver's entry point into all of this."""
    _, circuit = _chamber(
        cross_section=psf.regen.CrossSectionRounded(height_input="base"),
    )
    cached = circuit.R_coolant_per_len(0.5, h_c=1000.0, T_wall_rep=400.0)

    circuit._invalidate_section_cache()
    assert circuit.R_coolant_per_len(0.5, h_c=1000.0, T_wall_rep=400.0) == cached


def test_repeated_resistance_queries_do_not_rebuild_geometry(monkeypatch) -> None:
    _, circuit = _chamber()
    circuit.R_coolant_per_len(0.5, h_c=1000.0, T_wall_rep=400.0)

    builds = _count_calls(monkeypatch, psf.regen.CoolingCircuit, "_build_prof")
    stations = _count_calls(
        monkeypatch, psf.regen.CoolingCircuit, "_build_section_profiles_at"
    )
    # A wall node holds x fixed while its temperatures iterate, which is the
    # shape of the loop that made this expensive.
    for temperature in np.linspace(300.0, 900.0, 40):
        circuit.R_coolant_per_len(0.5, h_c=1000.0, T_wall_rep=float(temperature))
    assert not builds
    assert not stations
