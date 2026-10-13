"""Injectors: element sets on a face, and the flow from each feed through them.

Each element type is its own class. :class:`ImpingingElement` is the first:
jets aimed at a common point, which covers like and unlike doublets, triplets
and pentads. Holes are :class:`~pyskyfire.ducts.Orifice` objects, so sizing,
discharge coefficients, choking and cavitation come from :mod:`pyskyfire.ducts`.
This subpackage depends on :mod:`pyskyfire.ducts` and :mod:`pyskyfire.thermo`
only; the engine network adapts it through
:class:`~pyskyfire.network.InjectorBlock`.
"""

from .element import Element, Hole, JetFlow
from .impinging import ImpingingDiagnostics, ImpingingElement, Jet
from .injector import (ElementSet, FeedCondition, FeedResult, Injector, InjectorChecks,
                       InjectorHoles, InjectorResult, SetResult)

__all__ = ["Element", "Hole", "JetFlow", "Jet", "ImpingingElement", "ImpingingDiagnostics",
           "ElementSet", "Injector", "InjectorHoles", "FeedCondition", "FeedResult",
           "SetResult", "InjectorResult", "InjectorChecks"]
