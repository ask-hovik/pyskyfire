
class Material:
    """Class used to specify a material and its properties. For calculating temperatures, only 'k' must be defined. For stresses, you also need E, alpha, and poisson.

    Args:
        k (float): Thermal conductivity (W/m/K)
        
    Keyword Args:
        E (float): Young's modulus (Pa)
        alpha (float): Thermal expansion coefficient (strain/K)
        poisson (float): Poisson's ratio
    """
    def __init__(self, **kwargs): # TODO: make the different properties dependent on temperature like they are in real life.
        self.k = kwargs["k"]               
        self.E = kwargs["E"]
        self.alpha = kwargs["alpha"]
        self.poisson = kwargs["poisson"]
        self.rho = kwargs["rho"]

# Solids
CopperC106 = Material(
    E       = 117e9,
    poisson = 0.34,
    alpha   = 16.9e-6,
    k       = 391.2,
    rho     = 8960
)

GRCop42 = Material(
    E       = 117e9,
    poisson = 0.34,
    alpha   = 16.9e-6,
    k       = 312,
    rho     = 8960
)

StainlessSteel304 = Material(
    E       = 193e9,
    poisson = 0.29,
    alpha   = 16e-6,
    k       = 14.0,
    rho     = 8000
)

Graphite = Material(
    E       = float('NaN'),
    poisson = float('NaN'),
    alpha   = float('NaN'),
    k       = 63.81001,
    rho     = float('NaN')
)

Chromium = Material(
    E       = 140e9,
    poisson = 0.33,
    alpha   = 16e-6,
    k       = 96.069,
    rho     = 7190
)

ZirconiumOxide = Material(
    E       = float('NaN'),
    poisson = float('NaN'),
    alpha   = float('NaN'),
    k       = 2.2,
    rho     = 5600
)


# Fluids
#Water = TransportProperties(Pr = 6.159, mu = 0.89307e-3, k = 0.60627, cp = 4181.38, rho =  997.085)         # Water at 298 K and 1 bar [1]
#Ethanol = TransportProperties(Pr = 16.152, mu = 1.0855e-3, k = 0.163526, cp = 2433.31, rho = 785.26)        # Ethanol at 298 K and 1 bar [1]
#CO2 = TransportProperties(mu = 3.74e-5, k =  0.0737, Pr = 0.72)                                             # Representative values for CO2 gas



# ========================================
# Some old experimental classes under here
# ========================================

class CanteraTransport:
    """
    Cantera‑based transport / thermo wrapper for a *liquid‑vapor* fluid
    or an ideal‑solution mixture of such fluids.

    Parameters
    ----------
    yaml_files : str | list[str]
        Path(s) to YAML files that contain the species / phase you need.
        For the built‑ins just pass 'liquidvapor.yaml'.
    phase      : str
        The phase name to load (e.g. 'hydrogen', 'oxygen', or your own mixture).
    X          : dict, optional
        Mixture composition (mole fractions).  If omitted the phase
        must contain only one species.

    Notes
    -----
    *One* ct.Solution is created and cached; every call to `mu(T,P)` etc.
    just does a fast `TPX` reset (no re‑parsing, no equilibrate).
    """

    def __init__(self, yaml_files='liquidvapor.yaml', phase='hydrogen', X=None):
        self._gas = ct.Solution(yaml_files, phase=phase)
        # default composition = pure component
        self._X0  = X if X is not None else {self._gas.species_names[0]: 1.0}

        # quick cache for property tuples keyed on (rounded T, P)
        self._cache = {}

    # ---------- private helper ---------------------------------------------
    def _set_state(self, T, P):
        """Reset state & return key used for the small cache."""
        # round to 0.2 K / 100 Pa so we don't explode the cache
        key = (round(T, 1), round(P, 2))
        if key not in self._cache:
            self._gas.TPX = T, P, self._X0
            props = (
                self._gas.viscosity,                # µ  [Pa·s]
                self._gas.thermal_conductivity,     # k  [W/m/K]
                self._gas.cp_mass,                  # cp [J/kg/K]
                self._gas.density,                  # ρ  [kg/m³]
                self._gas.cp_mass / self._gas.cv_mass,  # γ
            )
            self._cache[key] = props
        return self._cache[key]

    # ---------- public API (same names as your old class) -------------------
    def mu(self, T, P):          # viscosity
        return self._set_state(T, P)[0]

    def k(self, T, P):           # thermal conductivity
        return self._set_state(T, P)[1]

    def cp(self, T, P):          # specific heat (constant P)
        return self._set_state(T, P)[2]

    def rho(self, T, P):         # density
        return self._set_state(T, P)[3]

    def gamma_coolant(self, T, P):  # cp/cv for compressible coolant
        return self._set_state(T, P)[4]

    def Pr(self, T, P):          # Prandtl number (computed, not stored)
        mu_, k_, cp_ = self.mu(T,P), self.k(T,P), self.cp(T,P)
        return cp_ * mu_ / k_

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

class NucleateBoiling:
    def __init__(self, vapour_transport, liquid_transport, sigma, h_fg, C_sf):
        """Class for representing the information needed to model nucleate boiling. Not currently used.

        Args:
            vapour_transport (TransportProperties): The transport properties of the vapour phase.
            liquid_transport (TransportProperties): The transport properties of the liquid phase.
            sigma (callable): Surface tension of the liquid-vapour interface (N/m), as a function of temperature (K) and pressure (Pa) in the form sigma(T,p).
            h_fg (callable): Enthalpy between vapour and liquid phases, as a function of pressure (Pa). h_fg = h_g - h_f. (J/kg/K)
            C_sf (float): Surface-fluid coefficient. Will be different for different material + fluid combinations. Some examples are available in References [4] and [6] given in bamboo.circuit.py.
        """

        raise ValueError("NucleateBoiling class is not yet implemented")
    
import cantera as ct

class CanteraFluidTransport:
    """
    Wraps a single ct.Solution('liquidvapor.yaml') so that you can
    query mixture properties of any combination of its species.
    """
    _sol = None

    def __init__(self, composition, basis='mass', mech='liquidvapor.yaml'):
        """
        composition : dict
            Either mass fractions (basis='mass') or mole fractions (basis='mole'),
            e.g. {'H2O':0.5, 'C2H5OH':0.5} or {'H2O':1, 'C2H5OH':1}.
        basis : 'mass' or 'mole'
        mech : path to YAML (default the built‑in liquidvapor.yaml)
        """
        # Lazy‑load the Solution once for all instances
        if CanteraFluidTransport._sol is None:
            CanteraFluidTransport._sol = ct.Solution(mech)
        self._gas = CanteraFluidTransport._sol

        # Normalize and store the desired composition vector
        species = list(self._gas.species_names)
        comp = {sp: 0.0 for sp in species}
        if basis.lower().startswith('m'):
            # interpret as mass fractions
            total = sum(composition.values())
            for sp, w in composition.items():
                comp[sp] = w/total
            self._basis = 'Y'
        else:
            # interpret as mole fractions
            total = sum(composition.values())
            for sp, x in composition.items():
                comp[sp] = x/total
            self._basis = 'X'
        self._comp = comp

    def _set_state(self, T, p):
        """Internal helper: set T, P, and X or Y on the Solution."""
        if self._basis == 'Y':
            self._gas.TPY = T, p, list(self._comp.values())
        else:
            self._gas.TPX = T, p, list(self._comp.values())

    def rho(self, T, p):
        self._set_state(T, p)
        return self._gas.density

    def cp(self, T, p):
        self._set_state(T, p)
        return self._gas.cp_mass

    def mu(self, T, p):
        self._set_state(T, p)
        return self._gas.viscosity

    def k(self, T, p):
        self._set_state(T, p)
        return self._gas.thermal_conductivity

    def Pr(self, T, p):
        """Prandtl number = cp·μ / k"""
        return self.cp(T,p)*self.mu(T,p)/self.k(T,p)

    def gamma(self, T, p):
        """If needed: cp/cv for compressible coolant."""
        self._set_state(T, p)
        return self._gas.cp_mass/self._gas.cv_mass
