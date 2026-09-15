"""Placement and volume accounting inside a prescribed thrust-chamber flow."""

from dataclasses import dataclass
import numpy as np
from .circuit import ChamberGas

__all__ = ["ChamberInstallation", "VolumeChange"]


@dataclass(frozen=True)
class VolumeChange:
    displaced: float
    envelope_total: float
    outside: float
    chamber_displaced: float
    """Volume removed upstream of the throat [m³]."""


class ChamberInstallation:
    """An exchanger in chamber coordinates, without gas-solution feedback.

    Volume is the sum of nonoverlapping swept outer envelopes, including
    passages. Users must avoid overlapping copies. Partially intersecting
    sections use deterministic triangle/axial quadrature, not a CAD Boolean.
    ``volume_change.chamber_displaced`` counts only x <= contour.x_t;
    ``displaced`` includes nozzle volume too.
    """

    def __init__(
        self,
        chamber,
        exchanger,
        *,
        position=(0.0, 0.0, 0.0),
        rotation=None,
        coefficient=None,
    ):
        self.chamber, self.exchanger = chamber, exchanger
        self.gas = ChamberGas(
            chamber, position=position, rotation=rotation, coefficient=coefficient
        )
        self.volume_change = self._volume_change()

    def solve(self, inlets, **kwargs):
        if "boundaries" in kwargs:
            raise ValueError(
                "a chamber installation supplies its own prescribed gas boundary"
            )
        if self.volume_change.outside > max(
            1e-12, 1e-7 * self.volume_change.envelope_total
        ):
            raise ValueError(
                "thermal solve requires fully contained channel envelopes; clipped installations support volume accounting only"
            )
        result = self.exchanger.solve(
            inlets,
            boundaries={name: self.gas for name in self.exchanger.circuits},
            **kwargs,
        )
        result.displaced_volume = self.volume_change.displaced
        return result

    def _volume_change(self):
        contour = self.chamber.contour
        xmin, xmax = float(np.min(contour.xs)), float(np.max(contour.xs))
        throat = float(contour.x_t)

        def inside(points, xmax_limit=xmax):
            x = points[..., 0]
            radii = contour.r(np.clip(x, xmin, xmax))
            return (
                (x >= xmin)
                & (x <= xmax_limit)
                & (np.linalg.norm(points[..., 1:], axis=-1) <= radii)
            )

        # Uniform-area triangle quadrature; refined only for clipped cells.
        bary = []
        n = 12
        for i in range(n):
            for j in range(n - i):
                bary.append(((i + 1 / 3) / n, (j + 1 / 3) / n))
                if i + j < n - 1:
                    bary.append(((i + 2 / 3) / n, (j + 2 / 3) / n))
        bary = np.asarray(bary)
        gx, gw = np.polynomial.legendre.leggauss(6)
        fractions = (gx + 1) / 2
        weights = gw / 2
        displaced = chamber_displaced = 0.0
        for c in self.exchanger.circuits.values():
            for placement in range(c.n_channels):
                sections = c.section_points(placement, samples=16)
                sections = sections @ self.gas.rotation.T + self.gas.position
                centers = (
                    c.placed_paths[placement].centerline @ self.gas.rotation.T
                    + self.gas.position
                )
                for i, length in enumerate(c.path.lengths):
                    cell_volume = 0.5 * (c.outer_area[i] + c.outer_area[i + 1]) * length
                    vertices = sections[i : i + 2]
                    # Exact full-envelope integral for cells safely contained.
                    if (
                        inside(vertices).all()
                        and inside((sections[i] + sections[i + 1]) / 2).all()
                    ):
                        displaced += cell_volume
                        if np.max(vertices[..., 0]) <= throat:
                            chamber_displaced += cell_volume
                            continue
                        if np.min(vertices[..., 0]) >= throat:
                            continue
                    else:
                        cell_volume = None
                    for fraction, weight in zip(fractions, weights):
                        points = (1 - fraction) * sections[i] + fraction * sections[
                            i + 1
                        ]
                        center = (1 - fraction) * centers[i] + fraction * centers[i + 1]
                        a, b = points - center, np.roll(points, -1, axis=0) - center
                        triangle_area = np.linalg.norm(np.cross(a, b), axis=1) / 2
                        samples = (
                            center
                            + a[:, None, :] * bary[None, :, 0, None]
                            + b[:, None, :] * bary[None, :, 1, None]
                        )
                        area = (1 - fraction) * c.outer_area[
                            i
                        ] + fraction * c.outer_area[i + 1]
                        patch_weight = triangle_area / triangle_area.sum()
                        if cell_volume is None:
                            displaced += (
                                length
                                * weight
                                * area
                                * np.dot(patch_weight, inside(samples).mean(axis=1))
                            )
                        chamber_displaced += (
                            length
                            * weight
                            * area
                            * np.dot(patch_weight, inside(samples, throat).mean(axis=1))
                        )
        total = self.exchanger.displaced_volume
        return VolumeChange(
            float(displaced),
            float(total),
            float(max(0, total - displaced)),
            float(chamber_displaced),
        )
