"""Inspect local exchanger sections without constructing a thrust chamber.

Run with: python examples/cross_sections/cross_sections.py
Dimensions are internal endpoint widths; all quantities use SI units.
"""

import numpy as np

from pyskyfire.regen import CrossSectionRounded, CrossSectionSquared, SectionProfiles


def build_sections():
    profiles = SectionProfiles.from_dimensions(
        floor_width=0.006, roof_width=0.006, height=0.010,
        wall_thickness=0.0005,
    )
    all_surfaces = dict(floor=1.0, roof=1.0, side_left=1.0, side_right=1.0)
    rounded = CrossSectionRounded(exposure=all_surfaces)
    squared = CrossSectionSquared(
        sidewall_half_thickness=0.0004, roof_thickness=0.0008,
        floor="line", roof="line", exposure=all_surfaces,
    )
    curved_profiles = SectionProfiles.from_dimensions(
        floor_width=0.006, roof_width=0.007, height=0.010,
        wall_thickness=0.0005, radius=0.040,
    )
    curved = CrossSectionSquared(
        sidewall_half_thickness=0.0004, floor="arc", roof="line",
        exposure={"floor": 0.75, "roof": 1.0},
    )
    sector_profiles = SectionProfiles.from_dimensions(
        radius=0.100, angle=0.100, height=0.015, wall_thickness=0.0005,
    )
    chamber = CrossSectionRounded(exposure_angles={"floor": 112*np.pi/180})
    return [("Rounded parallel sides", rounded, profiles),
            ("Rectangle", squared, profiles),
            ("Arc floor, straight roof", curved, curved_profiles),
            ("Rounded sector", chamber, sector_profiles)]


def main():
    for name, section, profiles in build_sections():
        lower, upper = section.widths(profiles)
        print(name)
        print(f"  floor/roof widths: {lower[0]*1000:.3f} / {upper[0]*1000:.3f} mm")
        print(f"  passage area: {section.A_coolant(profiles)[0]*1e6:.3f} mm²")
        print(f"  hydraulic diameter: {section.Dh_coolant(profiles)[0]*1000:.3f} mm")
        print(f"  active inner perimeter: {section.P_coolant(profiles)[0]*1000:.3f} mm")
        print(f"  area inside exterior boundary: {section.A_exterior(profiles)[0]*1e6:.3f} mm²")


if __name__ == "__main__":
    main()
