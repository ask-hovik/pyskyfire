class TransportProperties:
    def __init__(self, Pr, mu, k, cp = None, rho = None, gamma_coolant = None):
        """
        Container for specifying your transport properties. Each input can either be a function of temperature (K) and pressure (Pa) in that order, e.g. mu(T, p). Otherwise they can be constant floats.

        Args:
            Pr (float or callable): Prandtl number.
            mu (float or callable): Absolute viscosity (Pa s).
            k (float or callable): Thermal conductivity (W/m/K).
            cp (float or callable, optional): Isobaric specific heat capacity (J/kg/K) - only required for coolants.
            rho (float or callable, optional): Density (kg/m^3) - only required for coolants.
            gamma_coolant (float or callable, optional): Ratio of specific heats (cp/cv) for a compressible coolant. If this is submitted, it is assumed that this object represents a compressible coolant.
        
        Attributes:
            compressible_coolant (bool): Whether or not this TransportProperties object represents a compressible coolant.
        """

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

        Args:
            T (float): Temperature (K)
            p (float): Pressure (Pa)

        Returns:
            float: Prandtl number
        """
        if callable(self._Pr):
            return self._Pr(T, p)
        
        else:
            return self._Pr

    def mu(self, T, p):
        """Absolute viscosity (Pa s)

        Args:
            T (float): Temperature (K)
            p (float): Pressure (Pa)

        Returns:
            float: Absolute viscosity (Pa s)
        """
        if callable(self._mu):
            return self._mu(T, p)
        
        else:
            return self._mu

    def k(self, T, p):
        """Thermal conductivity (W/m/K)

        Args:
            T (float): Temperature (K)
            p (float): Pressure (Pa)

        Returns:
            float: Thermal conductivity (W/m/K)
        """
        if callable(self._k):
            return self._k(T, p)
        
        else:
            return self._k

    def rho(self, T, p):
        """Density (kg/m^3)
        Args:
            T (float): Temperature (K)
            p (float): Pressure (Pa)
        Returns:
            float: Density (kg/m^3)
        """
        if self._rho is None:
            raise ValueError("TransportProperties object does not have its density 'rho' defined. If you tried to use this TransportProperties object for a coolant, you need to specify the 'rho' input.")

        if callable(self._rho):
            return self._rho(T, p)
        
        else:
            return self._rho

    def cp(self, T, p):
        """Isobaric specific heat capacity (J/kg/K)

        Args:
            T (float): Temperature (K)
            p (float): Pressure (Pa)

        Returns:
            float: Isobaric specific heat capacity (J/kg/K)
        """

        if self._cp is None:
            raise ValueError("TransportProperties object does not have its isobaric specific heat capacity 'cp' defined. If you tried to use this TransportProperties object for a coolant, you need to specify the 'cp' input.")

        if callable(self._cp):
            return self._cp(T, p)
        
        else:
            return self._cp

    def gamma_coolant(self, T, p):
        """Ratio of specific heat capacities for a compressible coolant.

        Args:
            T (float): Temperature (K)
            p (float): Pressure (Pa)

        Returns:
            float: Ratio of specific heat capacities (cp/cv).
        """

        if self._gamma_coolant is None:
            raise ValueError("TransportProperties object does not have its compressibgle coolant gamma 'gamma_coolant' defined.")

        if callable(self._gamma_coolant):
            return self._gamma_coolant(T, p)
        
        else:
            return self._gamma_coolant

import CoolProp.CoolProp as CP       

class CoolantTransport:
    def __init__(self, propellant):
        self.propellant = propellant
        
    def get_Pr(self, T, p):
        return CP.PropsSI("PRANDTL", "T", T, "P", p, self.propellant)

    def get_mu(self, T, p):
        return CP.PropsSI("VISCOSITY", "T", T, "P", p, self.propellant)

    def get_k(self, T, p):
        return CP.PropsSI("CONDUCTIVITY", "T", T, "P", p, self.propellant)

    def get_cp(self, T, p):
        return CP.PropsSI("CPMASS", "T", T, "P", p, self.propellant)

    def get_rho(self, T, p):
        return CP.PropsSI("DMASS", "T", T, "P", p, self.propellant)
    
    def get_cv(self, T, p):
        return CP.PropsSI("CVMASS", "T", T, "P", p, self.propellant)
    
    def get_gamma(self, T, p):
        return CP.PropsSI("ISENTROPIC_EXPONENT", "T", T, "P", p, self.propellant)
