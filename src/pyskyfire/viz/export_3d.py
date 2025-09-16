import numpy as np
import gmsh

def make_engine_gmsh(thrust_chamber, filename="engine", display_channels=None):
    gmsh.initialize()
    gmsh.model.add(filename)

    thrust_chamber.build_channel_centerlines(mode="plot")
    for circuit in thrust_chamber.cooling_circuit_group.circuits:
        circuit.compute_single_centerline()

    # Build one tube (seed) from wires
    for circuit in range(len(thrust_chamber.cooling_circuit_group.circuits)):
        wire_list = thrust_chamber.cooling_circuit_group.circuits[circuit].wires
        tube = [dt for dt in gmsh.model.occ.addThruSections(wire_list, makeSolid=True) if dt[0] == 3][0]

        # Full pattern count (spacing basis)
        c = thrust_chamber.cooling_circuit_group.circuits[circuit]
        n_chan = c.placement.n_channel_positions * c.placement.n_channels_per_leaf

        # How many to actually display?
        # None or >= n_chan -> show all; else show first display_channels positions.
        if display_channels is None or display_channels >= n_chan:
            to_show = n_chan
        else:
            to_show = max(1, int(display_channels))  # at least the seed

        parts = [tube]  # seed is position k = 0
        # Copy only the next (to_show - 1) instances, at the same angular step as if all n_chan existed
        for k in range(1, to_show):
            cp = gmsh.model.occ.copy([tube])
            gmsh.model.occ.rotate(cp, 0, 0, 0, 0, 0, 1, 2 * np.pi * k / n_chan)
            parts.extend(cp)

    gmsh.model.occ.synchronize()
    gmsh.write(f"{filename}")
    gmsh.finalize()
    print(f"Wrote {filename} with {to_show}/{n_chan} channels placed.")
