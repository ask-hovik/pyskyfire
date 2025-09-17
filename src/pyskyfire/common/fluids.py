import CoolProp.CoolProp as CP

class Fluid:
    def __init__(self, type, propellants, fractions, basis="mass", precision=3):
        """
        Args:
            type (str): "fuel" | "oxidizer" | "coolant" etc.
            propellants (list[str]): Component names (CoolProp canonical names)
            fractions (list[float]): Fractions (mass or mole, depending on basis)
            basis (str): "mass" or "mole"
            precision (int): number of decimal digits when exporting mole fractions
        """
        self.type = type
        self.propellants = propellants
        self.fractions = fractions
        self.basis = basis
        self.precision = precision

    def molar_masses(self):
        """Return molar masses [kg/mol] for each propellant."""
        return [CP.PropsSI("M", p) for p in self.propellants]

    def as_mole_fractions(self):
        """Convert stored fractions to mole fractions."""
        if self.basis == "mole":
            mole = self.fractions
        elif self.basis == "mass":
            M = self.molar_masses()
            mole_basis = [f / M_i for f, M_i in zip(self.fractions, M)]
            total = sum(mole_basis)
            mole = [x / total for x in mole_basis]
        else:
            raise ValueError("basis must be 'mass' or 'mole'")

        # Normalize + round
        s = sum(mole)
        mole = [x / s for x in mole]
        return [round(x, self.precision) for x in mole]

    def coolprop_string(self):
        """Return a CoolProp HEOS string like 'HEOS::Ethanol[0.800]&Water[0.200]'."""
        mole_fracs = self.as_mole_fractions()
        parts = [f"{p}[{mf:.{self.precision}f}]" for p, mf in zip(self.propellants, mole_fracs)]
        return "HEOS::" + "&".join(parts)
