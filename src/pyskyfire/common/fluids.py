import CoolProp.CoolProp as CP

class Fluid:
    """Mixture helper for fluids.

    Stores a set of component names and composition.

    Parameters
    ----------
    type : str
        Role of the mixture (e.g., ``"fuel"``, ``"oxidizer"``, ``"coolant"``).
    propellants : list[str]
        Component names (CoolProp canonical names, e.g., ``"Ethanol"``, ``"Water"``).
    fractions : list[float]
        Composition fractions corresponding to ``propellants``. Interpreted
        as mass or mole fractions depending on ``basis``.
    basis : str, optional
        ``"mass"`` or ``"mole"``. Default is ``"mass"``.
    precision : int, optional
        Number of decimal digits when exporting mole fractions (e.g., in
        :meth:`coolprop_string`). Default is ``3``.

    Attributes
    ----------
    type : str
        Role marker for the mixture.
    propellants : list[str]
        CoolProp component names.
    fractions : list[float]
        Stored composition in the given ``basis``.
    basis : str
        Basis of the stored composition (``"mass"`` or ``"mole"``).
    precision : int
        Decimal places for exported mole fractions.

    See Also
    --------
    CoolProp
        External property library used for molar masses and HEOS string format.

    Notes
    -----
    Fractions are re-normalized internally when converting to mole basis to
    guard against small rounding errors.
    """
    def __init__(self, type, propellants, fractions, basis="mass", precision=3):
        self.type = type
        self.propellants = propellants
        self.fractions = fractions
        self.basis = basis
        self.precision = precision

    def molar_masses(self):
        """Return molar masses for each component.

        Returns
        -------
        list[float]
            Molar masses in ``kg/mol`` for each name in :attr:`propellants`.

        Notes
        -----
        Values are retrieved via ``CoolProp.PropsSI("M", name)``.
        """
        return [CP.PropsSI("M", p) for p in self.propellants]

    def as_mole_fractions(self):
        """Convert the stored composition to mole fractions.

        Returns
        -------
        list[float]
            Mole fractions for each component in :attr:`propellants`, rounded to
            :attr:`precision` decimals and re-normalized to sum to 1.0.

        Raises
        ------
        ValueError
            If :attr:`basis` is not ``"mass"`` or ``"mole"``.

        Notes
        -----
        For a mass-basis input, the conversion uses
        :math:`x_i = \\frac{w_i/M_i}{\\sum_j w_j/M_j}` with molar masses
        :math:`M_i` from :meth:`molar_masses`.
        """
        
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
        """Format a CoolProp HEOS mixture string.

        Returns
        -------
        str
            A string such as ``'HEOS::Ethanol[0.800]&Water[0.200]'``, where the
            bracketed values are **mole fractions** with :attr:`precision` digits.

        Notes
        -----
        The composition is always exported on a **mole** basis, using
        :meth:`as_mole_fractions` for conversion when needed.
        """
        mole_fracs = self.as_mole_fractions()
        parts = [f"{p}[{mf:.{self.precision}f}]" for p, mf in zip(self.propellants, mole_fracs)]
        return "HEOS::" + "&".join(parts)
