"""Impinging elements: jets drilled so that they meet at one point.

An impinging element is a set of jets aimed at a common impingement point a
height above the face. Like and unlike doublets, triplets and pentads are the
same class with different jets: a like doublet has two jets on one feed, an
unlike doublet one on each, an O-F-O triplet a fuel jet on the axis between two
oxidizer jets. Each jet's inclination is its own, so an unlike doublet can
lean its jets unequally to point the resultant momentum down the axis.

A jet's hole leaves the face at ``impingement_height·tan(inclination)`` from
the element centre, at its azimuth, and points at the impingement point.

The mixing quantities are the element's resultant momentum, its angle from the
chamber axis, and the element mixture ratio. An element with exactly one fuel
and one oxidizer jet also reports Rupe's ratio, ``ρ_o v_o² d_o / ρ_f v_f² d_f``,
which Rupe found gives the most uniform mixing near one, and the
momentum-flux ratio ``ρ_o v_o² / ρ_f v_f²``.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Hashable

import numpy as np

from pyskyfire.ducts.orifice import Orifice

from .element import Element, Hole, JetFlow

__all__ = ["Jet", "ImpingingElement", "ImpingingDiagnostics"]


def _diameter(orifice: Orifice) -> float:
    if orifice.hydraulic_diameter is not None:
        return orifice.hydraulic_diameter
    return math.sqrt(4.0 * orifice.area / math.pi)


@dataclass(frozen=True)
class Jet:
    """One jet of an impinging element.

    Parameters
    ----------
    feed : hashable
        What supplies the hole: feed hardware such as a
        :class:`~pyskyfire.ducts.Manifold`, or a name.
    orifice : Orifice
        The hole.
    inclination : float
        Angle between the jet and the chamber axis [rad], in [0, π/2). Zero
        is a jet down the element axis.
    azimuth : float, optional
        Where the hole sits around the element centre, from the element's
        orientation [rad].
    """

    feed: Hashable
    orifice: Orifice
    inclination: float
    azimuth: float = 0.0

    def __post_init__(self):
        hash(self.feed)
        if not isinstance(self.orifice, Orifice):
            raise TypeError("orifice must be a pyskyfire.ducts.Orifice")
        inclination = float(self.inclination)
        if not (math.isfinite(inclination) and 0.0 <= inclination < 0.5 * math.pi):
            raise ValueError("inclination must lie in [0, pi/2)")
        azimuth = float(self.azimuth)
        if not math.isfinite(azimuth):
            raise ValueError("azimuth must be finite")
        object.__setattr__(self, "inclination", inclination)
        object.__setattr__(self, "azimuth", azimuth)


@dataclass(frozen=True)
class ImpingingDiagnostics:
    """Mixing quantities of one impinging element.

    ``momentum`` is the magnitude of the summed jet momentum [N] and
    ``resultant_angle`` its angle from the chamber axis [rad].
    ``mixture_ratio`` is the element's oxidizer over fuel mass flow, NaN
    without fuel or without any oxidizer feed named. ``rupe`` and
    ``momentum_flux_ratio`` are oxidizer over fuel, and NaN unless the element
    has exactly one fuel and one oxidizer jet.
    """

    momentum: float
    resultant_angle: float
    mixture_ratio: float
    rupe: float
    momentum_flux_ratio: float


class ImpingingElement(Element):
    """Jets drilled to meet at a point ``impingement_height`` above the face.

    Parameters
    ----------
    jets : sequence of Jet
        The element's jets.
    impingement_height : float
        Distance from the face to the impingement point [m].
    """

    def __init__(self, jets, impingement_height):
        jets = tuple(jets)
        if not jets:
            raise ValueError("an impinging element needs at least one jet")
        if not all(isinstance(jet, Jet) for jet in jets):
            raise TypeError("jets must contain only Jet instances")
        height = float(impingement_height)
        if not (math.isfinite(height) and height > 0.0):
            raise ValueError("impingement_height must be finite and positive")
        self.jets = jets
        self.impingement_height = height

    def __repr__(self):
        return (f"ImpingingElement({len(self.jets)} jets, "
                f"impingement_height={self.impingement_height:g})")

    def holes(self) -> tuple[Hole, ...]:
        target = np.array([self.impingement_height, 0.0, 0.0])
        holes = []
        for jet in self.jets:
            offset = self.impingement_height * math.tan(jet.inclination)
            exit_ = np.array([0.0, offset * math.cos(jet.azimuth), offset * math.sin(jet.azimuth)])
            direction = (target - exit_) / np.linalg.norm(target - exit_)
            holes.append(Hole(jet.feed, jet.orifice, exit_, direction, target.copy()))
        return tuple(holes)

    def free_jet_length(self, jet: Jet) -> float:
        """Distance a jet travels from its hole to the impingement point [m]."""
        return self.impingement_height / math.cos(jet.inclination)

    def diagnostics(self, flows, oxidizer) -> ImpingingDiagnostics:
        flows = tuple(flows)
        if len(flows) != len(self.jets):
            raise ValueError("flows must give one JetFlow per jet")
        momentum = np.zeros(3)
        for hole, flow in zip(self.holes(), flows):
            momentum += flow.mass_flow * flow.velocity * hole.direction
        magnitude = float(np.linalg.norm(momentum))
        angle = math.atan2(float(np.linalg.norm(momentum[1:])), float(momentum[0]))

        ox = [(jet, flow) for jet, flow in zip(self.jets, flows) if jet.feed in oxidizer]
        fu = [(jet, flow) for jet, flow in zip(self.jets, flows) if jet.feed not in oxidizer]
        mdot_ox = sum(flow.mass_flow for _, flow in ox)
        mdot_fu = sum(flow.mass_flow for _, flow in fu)
        mixture_ratio = mdot_ox / mdot_fu if oxidizer and mdot_fu > 0.0 else math.nan

        rupe = flux_ratio = math.nan
        if len(ox) == 1 and len(fu) == 1:
            (jet_o, o), (jet_f, f) = ox[0], fu[0]
            flux_ratio = (o.density * o.velocity**2) / (f.density * f.velocity**2)
            rupe = flux_ratio * _diameter(jet_o.orifice) / _diameter(jet_f.orifice)
        return ImpingingDiagnostics(magnitude, angle, mixture_ratio, rupe, flux_ratio)
