"""Tests for numerically stable cooling-channel geometry helpers."""

import warnings

import numpy as np

from pyskyfire.regen.channel_height import make_channel_height_fn


class _Contour:
    xs = np.array([-1.0, 0.0, 1.0])
    rs = np.array([2.0, 1.0, 2.0])

    @staticmethod
    def r(x):
        return 1.0 + np.abs(x)


def test_narrow_channel_height_transition_does_not_overflow() -> None:
    channel_height = make_channel_height_fn(
        contour=_Contour(),
        region_fractions=[-1.0, 0.0, 1.0],
        flat_heights=[0.0025, 0.004],
        pinch_factors=[-1.0, -4.0],
        transition_widths=[1.0e-6],
    )

    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        heights = [channel_height(-1.0), channel_height(1.0)]

    assert np.all(np.isfinite(heights))
