from pyskyfire.common.fluids import Fluid

class TransportProperties:
    """Container for transport-property models (constants or callables).

    This class stores transport properties that can be either constant scalars
    or callables of temperature and pressure in the order ``(T [K], p [Pa])``.
    It is suitable for **walls**, **coolants**, or generic working fluids. When
    ``gamma_coolant`` is supplied, the instance is treated as a *compressible
    coolant*.

    Parameters
    ----------
    Pr : float or callable
        Prandtl number [-] or a function ``Pr(T, p) -> float``.
    mu : float or callable
        Dynamic (absolute) viscosity [Pa·s] or a function ``mu(T, p) -> float``.
    k : float or callable
        Thermal conductivity [W/(m·K)] or a function ``k(T, p) -> float``.
    cp : float or callable, optional
        Isobaric specific heat [J/(kg·K)] or ``cp(T, p) -> float``. Required when
        this object represents a coolant.
    rho : float or callable, optional
        Density [kg/m³] or ``rho(T, p) -> float``. Required when this object
        represents a coolant.
    gamma_coolant : float or callable, optional
        Ratio of specific heats ``γ = cp/cv`` [-] or ``gamma(T, p) -> float``.
        If provided, the object is assumed to represent a **compressible coolant**.

    Attributes
    ----------
    compressible_coolant : bool
        Whether the instance represents a compressible coolant (``gamma_coolant``
        was provided).
    _Pr, _mu, _k, _cp, _rho, _gamma_coolant : Any
        Backing values/callables as provided.

    See Also
    --------
    :class:`~pyskyfire.common.fluids.Fluid`
        Fluid identifier utilities used elsewhere in the package.

    Notes
    -----
    - All callables must accept ``(T [K], p [Pa])`` in that order.
    - For coolants, both ``cp`` and ``rho`` should be provided to allow
    consistent thermo-fluid calculations.
    """

    def __init__(self, Pr, mu, k, cp = None, rho = None, gamma_coolant = None):

        self.type = type
        self._Pr = Pr
        self._mu = mu
        self._k = k
        self._rho = rho
        self._cp = cp
        self._gamma_coolant = gamma_coolant

        if gamma_coolant is None:
            self.compressible_coolant = False
        else:
            self.compressible_coolant = True

    def Pr(self, T, p):
        """Prandtl number.

        Parameters
        ----------
        T : float
            Temperature [K].
        p : float
            Pressure [Pa].

        Returns
        -------
        float
            Prandtl number [-].
        """

        if callable(self._Pr):
            return self._Pr(T, p)
        
        else:
            return self._Pr

    def mu(self, T, p):
        """Dynamic (absolute) viscosity.

        Parameters
        ----------
        T : float
            Temperature [K].
        p : float
            Pressure [Pa].

        Returns
        -------
        float
            Viscosity [Pa·s].
        """

        if callable(self._mu):
            return self._mu(T, p)
        
        else:
            return self._mu

    def k(self, T, p):
        """Thermal conductivity.

        Parameters
        ----------
        T : float
            Temperature [K].
        p : float
            Pressure [Pa].

        Returns
        -------
        float
            Thermal conductivity [W/(m·K)].
        """

        if callable(self._k):
            return self._k(T, p)
        
        else:
            return self._k

    def rho(self, T, p):
        """Density.

        Parameters
        ----------
        T : float
            Temperature [K].
        p : float
            Pressure [Pa].

        Returns
        -------
        float
            Density [kg/m³].

        Raises
        ------
        ValueError
            If ``rho`` was not provided at construction.
        """

        if self._rho is None:
            raise ValueError("TransportProperties object does not have its density 'rho' defined. If you tried to use this TransportProperties object for a coolant, you need to specify the 'rho' input.")

        if callable(self._rho):
            return self._rho(T, p)
        
        else:
            return self._rho

    def cp(self, T, p):
        """Isobaric specific heat capacity.

        Parameters
        ----------
        T : float
            Temperature [K].
        p : float
            Pressure [Pa].

        Returns
        -------
        float
            ``c_p`` [J/(kg·K)].

        Raises
        ------
        ValueError
            If ``cp`` was not provided at construction.
        """

        if self._cp is None:
            raise ValueError("TransportProperties object does not have its isobaric specific heat capacity 'cp' defined. If you tried to use this TransportProperties object for a coolant, you need to specify the 'cp' input.")

        if callable(self._cp):
            return self._cp(T, p)
        
        else:
            return self._cp

    def gamma_coolant(self, T, p):
        """Ratio of specific heats for a compressible coolant.

        Parameters
        ----------
        T : float
            Temperature [K].
        p : float
            Pressure [Pa].

        Returns
        -------
        float
            ``γ = c_p / c_v`` [-].

        Raises
        ------
        ValueError
            If ``gamma_coolant`` was not provided at construction.
        """


        if self._gamma_coolant is None:
            raise ValueError("TransportProperties object does not have its compressibgle coolant gamma 'gamma_coolant' defined.")

        if callable(self._gamma_coolant):
            return self._gamma_coolant(T, p)
        
        else:
            return self._gamma_coolant

import CoolProp.CoolProp as CP       

class CoolantTransport:
    """CoolProp-backed transport/thermo properties for a single fluid.

    Thin wrapper around **CoolProp** that exposes convenient getters using a
    :class:`~pyskyfire.common.fluids.Fluid` to obtain the CoolProp fluid string.

    Parameters
    ----------
    fluid : pyskyfire.common.fluids.Fluid
        Fluid descriptor; its ``.coolprop_string()`` is passed to CoolProp.

    Attributes
    ----------
    fluid : str
        CoolProp fluid identifier string.

    See Also
    --------
    CoolProp.CoolProp.PropsSI
        Underlying property evaluator.
    :class:`~pyskyfire.common.fluids.Fluid`
        Source of the CoolProp-compatible fluid string.

    Notes
    -----
    - All getters accept ``(T [K], p [Pa])`` in that order and return SI units.
    """

    def __init__(self, fluid: Fluid):
        self.fluid = fluid.coolprop_string()
        
    def get_Pr(self, T, p):
        return CP.PropsSI("PRANDTL", "T", T, "P", p, self.fluid)

    def get_mu(self, T, p):
        return CP.PropsSI("VISCOSITY", "T", T, "P", p, self.fluid)

    def get_k(self, T, p):
        return CP.PropsSI("CONDUCTIVITY", "T", T, "P", p, self.fluid)

    def get_cp(self, T, p):
        return CP.PropsSI("CPMASS", "T", T, "P", p, self.fluid)

    def get_rho(self, T, p):
        return CP.PropsSI("DMASS", "T", T, "P", p, self.fluid)
    
    def get_cv(self, T, p):
        return CP.PropsSI("CVMASS", "T", T, "P", p, self.fluid)
    
    def get_gamma(self, T, p):
        return CP.PropsSI("ISENTROPIC_EXPONENT", "T", T, "P", p, self.fluid)
