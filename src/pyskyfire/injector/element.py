"""What every injector element type offers the injector it sits in.

An element is one repeated unit on the faceplate: an impinging doublet, a
coaxial swirl element, a pintle. Each type has its own class, its own
parameters and its own modelling. The injector asks each of them only two
things: where its holes are and which feed supplies each one
(:meth:`Element.holes`), and what its flows mean for mixing
(:meth:`Element.diagnostics`).

A *feed* is whatever supplies a hole: a :class:`~pyskyfire.ducts.Manifold`,
another piece of feed hardware, or a plain name where the hardware is not
modelled yet. Feeds are only ever compared and used as dictionary keys, so any
hashable object works. Hardware objects compare by identity, so two physically
separate manifolds stay separate even if they are built alike.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Hashable

import numpy as np

from pyskyfire.ducts.orifice import Orifice

__all__ = ["Element", "Hole", "JetFlow"]


@dataclass(frozen=True)
class Hole:
    """One hole of an element, placed in the element's own frame.

    The frame has its origin at the element centre on the face, +x along the
    chamber axis into the chamber, and +y along the element's orientation.
    ``exit`` is where the hole meets the face [m], ``direction`` a unit
    vector along the discharge, and ``target`` where the jet is drawn to and
    checked against the wall: the impingement point for an impinging element.
    """

    feed: Hashable
    orifice: Orifice
    exit: np.ndarray
    direction: np.ndarray
    target: np.ndarray


@dataclass(frozen=True)
class JetFlow:
    """The flow out of one hole: mass flow [kg/s], velocity [m/s], density [kg/m³]."""

    mass_flow: float
    velocity: float
    density: float


class Element(ABC):
    """Base class for injector element types."""

    @abstractmethod
    def holes(self) -> tuple[Hole, ...]:
        """Every hole of one element, in the element frame."""

    @abstractmethod
    def diagnostics(self, flows: tuple[JetFlow, ...], oxidizer: frozenset):
        """Mixing quantities of one element.

        ``flows`` is aligned with :meth:`holes`. ``oxidizer`` is the set of
        feeds carrying oxidizer; every other feed carries fuel.
        """

    @property
    def feeds(self) -> tuple:
        """The feeds supplying this element, in the order its holes list them."""
        return tuple(dict.fromkeys(hole.feed for hole in self.holes()))

    def orifices(self, feed) -> list[Orifice]:
        """Every orifice fed from ``feed``: one entry per hole, repeats included."""
        return [hole.orifice for hole in self.holes() if hole.feed == feed]
