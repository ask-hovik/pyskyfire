# Regenerative Cooling

This page explains the regenerative cooling model implemented in `pyskyfire.regen`. It traces the path through `solver.py`, shows how `solver.py` calls the equations in `physics.py`, and explains how the hot-gas aerothermodynamic properties are prepared by `pyskyfire.skycea.aerothermodynamics`.

The model is a quasi-one-dimensional steady heat-exchanger calculation. At each axial station, it balances three heat-transfer paths:

1. heat transfer from the combustion gas to the hot wall,
2. conduction through one or more wall layers, and
3. heat transfer from the cold wall into the coolant flowing in the cooling channels.

The result is a streamwise solution for wall temperatures, coolant bulk temperature, coolant pressure, heat flux, heat-transfer coefficients, coolant velocity, and residuals.

---

## 1. Code path overview

The public entry point is:

```python
steady_heating_analysis(
    thrust_chamber,
    boundary_conditions,
    n_nodes=100,
    circuit_index=0,
    solver="newton",
    output=True,
)
```

For the current implementation, the only accepted solver name is `"newton"`. This dispatches to:

```python
solve_heat_exchanger(
    thrust_chamber,
    boundary_conditions,
    n_nodes,
    circuit_index,
    output,
)
```

The coolant temperature and pressure are marched explicitly, but the two wall temperatures are obtained from a local nonlinear heat-balance solve using `scipy.optimize.least_squares`.

The helper class:

```python
HeatExchangerPhysics(thrust_chamber, boundary_conditions, circuit_index)
```

collects the local physics calculations used by the marching solver:

- `hot_side_coefficients(x, T_hw)`
- `cold_side_coefficients(x, T_cw, T_cool)`
- `dQ_hot_dx(x, T_hw)`
- `dQ_cond_dx(x, T_hw, T_cw)`
- `dQ_cold_dx(x, T_cw, T_cool)`
- `coolant_temperature_rate(T_cool, p_cool, dQ_cold_dx)`
- `coolant_pressure_rate(x, T_cool, p_cool)`
- `interface_temperatures(x, T_hw, T_cw)`

These methods call lower-level correlations from `physics.py`, mainly Bartz-type hot-gas heat transfer, Colburn coolant-side heat transfer, Reynolds number, Darcy friction factor, coolant velocity, curvature factor, and adiabatic-wall temperature.

The hot-gas property calls such as:

```python
combustion_transport.get_T(x)
combustion_transport.get_p(x)
combustion_transport.get_h(x)
combustion_transport.get_cp(x)
combustion_transport.get_mu(x)
combustion_transport.get_k(x)
combustion_transport.get_Pr(x)
combustion_transport.get_M(x)
combustion_transport.get_a(x)
combustion_transport.get_gamma(x)
```

are normally supplied by `skycea.aerothermodynamics.Aerothermodynamics`, which precomputes CEA-based equilibrium and temperature-pressure maps along the chamber/nozzle contour.

---

## 2. Symbols and sign conventions

The most important local variables are:

| Symbol | Meaning | Units |
|---|---:|---:|
| `x` | axial coordinate | m |
| `dx` | axial node spacing used by the marching solver | m |
| `T_hw` | hot-side wall temperature | K |
| `T_cw` | coolant-side wall temperature | K |
| `T_cool` | coolant bulk static temperature | K |
| `p_static` | coolant static pressure | Pa |
| `p_stagnation` | coolant stagnation pressure | Pa |
| `dA_hot/dx` | hot-side thermal area per unit axial length | m |
| `A_cool` | flow area of one coolant channel | m² |
| `D_h` | coolant hydraulic diameter | m |
| `D_hyd` | hot-gas-side hydraulic diameter approximation | m |
| `mdot_c` | coolant mass flow through one channel | kg/s |
| `mdot_g` | hot-gas mass flow | kg/s |
| `h_hot` | effective temperature-based hot-side coefficient | W/m²/K |
| `h_g` | enthalpy-based hot-side coefficient | kg/m²/s |
| `h_cold` | coolant-side heat-transfer coefficient | W/m²/K |
| `qpp_hot` | hot-side wall heat flux | W/m² |
| `dQ_hot_dx` | hot-side heat input per unit axial length | W/m |
| `dQ_cond_dx` | wall conduction per unit axial length | W/m |
| `dQ_cold_dx` | coolant heat pickup per unit axial length | W/m |

Heat flow is positive from the combustion gas into the wall and then into the coolant.

---

## 3. Boundary conditions

The `BoundaryConditions` class stores the inlet state of one cooling circuit:

```python
BoundaryConditions(T_coolant_in, p_coolant_in, mdot_coolant)
```

where:

$$
T_{c,in} = \texttt{T\_coolant\_in}
$$

$$
p_{c,in} = \texttt{p\_coolant\_in}
$$

$$
\dot{m}_{c,total} = \texttt{mdot\_coolant}
$$

The mass flow supplied here is the total mass flow through the selected cooling circuit. The solver divides it evenly across all geometric channel instances:

$$
N_{chan} = N_{positions} N_{channels/leaf}
$$

$$
\dot{m}_{c,chan} = \frac{\dot{m}_{c,total}}{N_{chan}}
$$

This single-channel mass flow is used for local velocity, Reynolds number, heat-transfer coefficient, coolant temperature rise, and pressure drop.

---

## 4. Aerothermodynamics from `skycea.aerothermodynamics`

The regenerative cooling solver needs the hot-gas state along the contour. In the CEA-based workflow, this is provided by the `Aerothermodynamics` class.

### 4.1 Initialization from thrust, area ratio, and characteristic length

One constructor is:

```python
Aerothermodynamics.from_F_eps_Lstar(
    fu, ox, MR, p_c, F, eps, L_star,
    T_fu_in=298.15,
    T_ox_in=298.15,
    p_amb=1.013e5,
    npts=15,
)
```

The oxidizer-to-fuel mixture ratio is:

$$
MR = \frac{\dot{m}_{ox}}{\dot{m}_{fu}}
$$

The method builds CEA fuel and oxidizer objects at the specified inlet temperatures and runs a `CEA_Wrap.RocketProblem` at chamber pressure and exit area ratio. Internally, chamber pressure is converted from pascal to psi using:

$$
p[\mathrm{psi}] = 0.000145038\,p[\mathrm{Pa}]
$$

CEA returns design-point quantities such as characteristic velocity, vacuum specific impulse, chamber density, chamber temperature, throat temperature, and throat pressure. These are then used to derive engine-level quantities.

The total mass flow is calculated from the requested thrust and vacuum specific impulse:

$$
\dot{m} = \frac{F}{I_{sp,vac} g_0}
$$

using:

$$
g_0 = 9.81~\mathrm{m/s^2}
$$

The fuel and oxidizer mass flows are:

$$
\dot{m}_{fu} = \frac{\dot{m}}{1 + MR}
$$

$$
\dot{m}_{ox} = \dot{m} - \dot{m}_{fu}
$$

The throat area follows from the definition of characteristic velocity:

$$
c^* = \frac{p_c A_t}{\dot{m}}
$$

so that:

$$
A_t = \frac{c^* \dot{m}}{p_c}
$$

The throat radius is:

$$
r_t = \sqrt{\frac{A_t}{\pi}}
$$

The exit area and exit radius are:

$$
A_e = \varepsilon A_t
$$

$$
r_e = \sqrt{\frac{A_e}{\pi}}
$$

The chamber residence time is computed from the characteristic length, throat area, chamber density, and total mass flow:

$$
t_{stay} = \frac{L^* A_t \rho_c}{\dot{m}}
$$

The chamber volume estimate is then:

$$
V_c = \frac{\dot{m} t_{stay}}{\rho_c}
$$

which reduces to:

$$
V_c = L^* A_t
$$

when substituting the previous expression.

The vacuum thrust coefficient is derived from vacuum specific impulse:

$$
C_{F,vac} = \frac{I_{sp,vac} g_0}{c^*}
$$

Ambient thrust coefficient is modeled by subtracting the pressure-thrust penalty:

$$
C_{F,amb} = C_{F,vac} - \frac{p_{amb}}{p_c}\frac{A_e}{A_t}
$$

and the corresponding ambient specific impulse is:

$$
I_{sp,amb} = \frac{C_{F,amb} c^*}{g_0}
$$

The sea-level thrust coefficient and sea-level specific impulse use:

$$
p_{SL} = 1.01325\times10^5~\mathrm{Pa}
$$

$$
C_{F,SL} = C_{F,vac} - \frac{p_{SL}}{p_c}\frac{A_e}{A_t}
$$

$$
I_{sp,SL} = \frac{C_{F,SL} c^*}{g_0}
$$

### 4.2 Initialization from thrust and exit pressure

The alternative constructor is:

```python
Aerothermodynamics.from_F_pe_Lstar(
    fu, ox, MR, p_c, F, p_e, L_star,
    T_fu_in=298.15,
    T_ox_in=298.15,
    p_amb=1.013e5,
    npts=15,
)
```

Instead of specifying exit area ratio directly, it specifies exit pressure. CEA is called with the pressure ratio:

$$
\frac{p_c}{p_e}
$$

The area ratio is then obtained from the CEA result:

$$
\varepsilon = \frac{A_e}{A_t}
$$

After that, the same mass-flow, area, radius, residence-time, and thrust-coefficient equations are used.

### 4.3 Precomputing the property maps

After construction, the method:

```python
compute_aerothermodynamics(contour, Nt=64)
```

builds a two-dimensional property table along the contour.

First, the axial grid is created:

$$
x_i = \mathrm{linspace}(x_{min}, x_{max}, n_{pts})
$$

where:

$$
x_{min} = \texttt{contour.xs[0]}
$$

$$
x_{max} = \texttt{contour.xs[-1]}
$$

The local area ratio is:

$$
\varepsilon_i = \frac{A(x_i)}{A_t}
$$

At each station, CEA is called in subsonic or supersonic mode depending on the sign of `x`:

$$
x_i < 0 \Rightarrow \text{subsonic CEA solve at } \varepsilon_i
$$

$$
x_i \ge 0 \Rightarrow \text{supersonic CEA solve at } \varepsilon_i
$$

The resulting equilibrium state is stored in column zero of every property map:

$$
M_i,\ T_i,\ p_i,\ \rho_i,\ c_{p,i},\ \gamma_i,\ h_i,\ a_i,\ \mu_i,\ k_i,\ Pr_i,\ MW_i
$$

CEA pressure is stored in bar in `p_map`; the getter later converts it to pascal:

$$
p[\mathrm{Pa}] = 10^5 p[\mathrm{bar}]
$$

For each axial row, a local temperature grid is built from the equilibrium temperature down to 200 K:

$$
T_{i,j} = \mathrm{linspace}(T_{eq,i},\ 200~\mathrm{K},\ N_T)
$$

For each nonzero temperature-grid column, a CEA temperature-pressure problem is solved at fixed local pressure:

$$
\text{TP solve at } \left(T_{i,j},\ p_i\right)
$$

The same property maps are filled with the TP result. The maps therefore contain:

- column 0: equilibrium nozzle/chamber solution at local area ratio,
- columns 1 through `Nt - 1`: TP solutions at the same local static pressure but different temperatures.

Implementation note: in the current code, the `Nt` argument is overwritten internally by:

```python
Nt = len(self.x_nodes)
```

so the number of temperature columns becomes equal to the number of axial stations, not necessarily the value passed as `Nt`.

### 4.4 Getter behavior

The getter methods all call the same internal logic. With no temperature or enthalpy argument, the getter interpolates the equilibrium column along `x`:

$$
Z(x) = (1-w_x) Z_i + w_x Z_{i+1}
$$

where:

$$
w_x = \frac{x - x_i}{x_{i+1} - x_i}
$$

This is used for calls such as:

```python
get_T(x)
get_p(x)
get_h(x)
get_cp(x)
get_gamma(x)
get_M(x)
get_a(x)
get_mu(x)
get_k(x)
get_Pr(x)
```

When a temperature is supplied, for example:

```python
get_h(x, T=T_wall)
```

bilinear interpolation is first tried over the precomputed `(x, T)` map. For each of the two bracketing axial rows, it interpolates in temperature:

$$
Z_i(T) = \mathrm{interp}\left(T;\ T_{i,:},\ Z_{i,:}\right)
$$

$$
Z_{i+1}(T) = \mathrm{interp}\left(T;\ T_{i+1,:},\ Z_{i+1,:}\right)
$$

Then it interpolates in `x`:

$$
Z(x,T) = (1-w_x)Z_i(T) + w_x Z_{i+1}(T)
$$

If the requested temperature lies outside the local table range, a live CEA TP solve is performed at:

$$
(T,\ p(x))
$$

where `p(x)` is taken from the equilibrium pressure column.

Pressure is special: `get_p(x, T=..., h=...)` ignores `T` and `h` and always returns the interpolated equilibrium pressure:

$$
p(x) = 10^5\,p_{map,0}(x)
$$

### 4.5 Composition interpolation

Equilibrium product mole fractions are stored as species dictionaries. When composition is requested at a point between two CEA stations, the code forms the union of species in the two bracketing dictionaries and linearly blends each species:

$$
X_s(x) = (1-w_x)X_{s,i} + w_x X_{s,i+1}
$$

Missing species are treated as zero. Negative values are clipped to zero, and the final dictionary is renormalized:

$$
X_s^{norm}(x) = \frac{\max(X_s(x),0)}{\sum_r \max(X_r(x),0)}
$$

### 4.6 Enthalpy-pressure lookup path

The HP lookup path is marked as under construction in the source code. It is used when a getter is called with `h=...`, for example:

```python
get_T(x, h=h_target)
```

The intended target is a CEA HP solve at:

$$
(h_{target},\ p(x))
$$

To do this, the implementation assigns the mixture enthalpy to one selected “adjuster” fuel species. First, fuel and oxidizer component weights are normalized within their streams:

$$
\hat{w}_{fu,m} = \frac{w_{fu,m}}{\sum_r w_{fu,r}}
$$

$$
\hat{w}_{ox,m} = \frac{w_{ox,m}}{\sum_r w_{ox,r}}
$$

The stream mass fractions implied by mixture ratio are:

$$
w_{fu,stream} = \frac{1}{1+MR}
$$

$$
w_{ox,stream} = \frac{MR}{1+MR}
$$

The overall mass fraction of each fuel component is:

$$
w_m = w_{fu,stream}\hat{w}_{fu,m}
$$

and the overall mass fraction of each oxidizer component is:

$$
w_m = w_{ox,stream}\hat{w}_{ox,m}
$$

The first fuel species is chosen as the adjuster. If its molecular mass is `MW_adj` and its overall mass fraction is `w_adj`, the assigned molar enthalpy is chosen such that:

$$
h_{target} = w_{adj}\frac{H_{adj}}{MW_{adj}}
$$

so:

$$
H_{adj} = \frac{h_{target} MW_{adj}}{w_{adj}}
$$

All other species receive zero assigned formation enthalpy in this artificial HP setup. CEA is then run at:

$$
(p(x),\ H_{adj})
$$

The molecular mass used for the adjuster is computed from an exploded chemical formula such as `"C 2 H 6 O 1"`:

$$
MW = \sum_e n_e MW_e
$$

where `n_e` is the integer count of element `e` in the formula.

---

## 5. Marching grid and unknowns in `solver.py`

The cooling solver reads the selected cooling circuit:

```python
circuit = thrust_chamber.cooling_circuits[circuit_index]
```

and builds a marching grid from `circuit.x_domain`. If the circuit direction is positive:

$$
x_i = \mathrm{linspace}(x_0, x_N, n_{nodes})
$$

If the circuit direction is negative, the order is reversed:

$$
x_i = \mathrm{linspace}(x_N, x_0, n_{nodes})
$$

The axial spacing used in the solver is:

$$
\Delta x = \left|\frac{x_{N-1}-x_0}{n_{nodes}-1}\right|
$$

The stored solution arrays are:

- `T_hw_arr`: hot-side wall temperature,
- `T_cw_arr`: coolant-side wall temperature,
- `T_cool_arr`: coolant bulk static temperature,
- `p_static_arr`: coolant static pressure,
- `p_stagnation_arr`: coolant stagnation pressure,
- `dQ_dA_arr`: local heat flux,
- `velocity_arr`: coolant velocity,
- `T_stagnation_arr`: coolant stagnation temperature.

The inlet values are:

$$
T_{cool,0} = T_{c,in}
$$

$$
p_{static,0} = p_{c,in}
$$

$$
p_{0,0} = p_{c,in}
$$

The initial guesses for both wall temperatures at the first station are:

$$
T_{hw,guess} = \frac{T_g(x_0)+T_{c,in}}{2}
$$

$$
T_{cw,guess} = \frac{T_g(x_0)+T_{c,in}}{2}
$$

At each station, the unknowns solved by least squares are:

$$
T_{cw}
$$

and:

$$
\Delta T_w = T_{hw} - T_{cw}
$$

so that:

$$
T_{hw} = T_{cw} + \Delta T_w
$$

The lower bound enforces:

$$
T_{cw} \ge T_{cool}
$$

and:

$$
\Delta T_w \ge 0
$$

so the wall ordering remains:

$$
T_{hw} \ge T_{cw} \ge T_{cool}
$$

The hot-wall upper bound is based on the local gas temperature:

$$
T_{hw,max} = \max(T_{cool,i}+1,\ T_{g,i}-10^{-3})
$$

---

## 6. Hot-gas-side heat transfer

The hot-side calculation is performed by:

```python
HeatExchangerPhysics.hot_side_coefficients(x, T_hw)
```

### 6.1 Local gas-side geometry and state

The hot-gas-side hydraulic diameter is approximated as twice the local chamber radius:

$$
D_{hyd,g}(x) = 2r(x)
$$

The local hot-gas flow area is:

$$
A_g(x) = A(x)
$$

The hot-gas mass flow is:

$$
\dot{m}_g = \texttt{combustion\_transport.mdot}
$$

The local gas static temperature is:

$$
T_g = T_g(x)
$$

The temperature used in the current Bartz-property correction is the arithmetic mean of wall and gas temperature:

$$
T_{gr} = \frac{T_{hw}+T_g}{2}
$$

The gas enthalpy at the equilibrium state is:

$$
H_g = H_g(x)
$$

The hot-wall enthalpy is requested from the combustion-transport model at the wall temperature:

$$
H_{hw} = H(x,T_{hw})
$$

If this lookup fails, the implementation falls back to:

$$
H_{hw} = H_g
$$

The local Mach number and speed of sound are:

$$
M_g = M(x)
$$

$$
a_g = a(x)
$$

A reference enthalpy is computed:

$$
H_{gr} = \frac{1}{2}(H_{hw}+H_g) + 0.18\left(\frac{1}{2}M_g^2 a_g^2\right)
$$

Implementation note: 

Currently `H_gr`is calculated, but not used to query reference condition gas properties. The actual properties used in the Bartz correlation are currently retrieved at `x` as shown below. The reference-condition properties should ideally be evaluated at `H_gr`, but the current aerothermodynamics module does not yet support that workflow robustly. This is not trivial to solve, as it generally requires us to move away from NASA CEA as a reference lookup. There are currently no other programs that support as wide a library of propellants, so ditching it would hurt propellant compatibility strongly.

$$
c_{p,gr} \leftarrow c_p(x)
$$

$$
\mu_{gr} \leftarrow \mu(x)
$$

$$
k_{gr} \leftarrow k(x)
$$

$$
Pr_{gr} \leftarrow Pr(x)
$$

### 6.2 Bartz-type enthalpy-driven coefficient

The hot-side coefficient is calculated by `physics.h_gas_bartz_enthalpy_driven`:

$$
h_{gr} = 0.026\frac{k_{gr}}{D_{hyd,g}}
\left(\frac{c_{p,gr}}{k_{gr}\mu_{gr}}\right)^{0.4}
\left(\frac{\dot{m}_g D_{hyd,g}}{A_g}\right)^{0.8}
\left(\frac{T_g}{T_{gr}}\right)^{0.8}
$$

The solver multiplies this by a user-supplied hot-side correction factor:

$$
h_{gr,eff} = h_{gr}\,C_{hot}
$$

where:

$$
C_{hot} = \texttt{thrust\_chamber.h\_gas\_corr}
$$

The enthalpy-based heat-transfer coefficient used in the heat-flux equation is:

$$
h_g = \frac{h_{gr,eff}}{c_{p,gr}}
$$

Because `h_gr` has units of W/m²/K and `c_p` has units of J/kg/K, `h_g` has units of kg/m²/s. This is why the hot-side heat flux is written in terms of enthalpy difference rather than temperature difference.

### 6.3 Adiabatic-wall temperature and adiabatic-wall enthalpy

Adiabatic-wall temperature is computed using `physics.T_aw`:

$$
r = Pr^{1/3}
$$

$$
T_{aw} = T_\infty\left(1 + r\frac{\gamma-1}{2}M_\infty^2\right)
$$

In the hot-side call, this becomes:

$$
T_{aw} = T_g\left(1 + Pr_{gr}^{1/3}\frac{\gamma-1}{2}M_g^2\right)
$$

The heat flux itself is driven by adiabatic-wall enthalpy, not directly by `T_aw`. We compute:

$$
H_{aw} = H_g + \frac{1}{2}Pr_{gr}^{1/3}M_g^2 a_g^2
$$

### 6.4 Hot-side heat flux and heat per unit length

The hot-side wall heat flux is:

$$
q''_{hot} = h_g\left(H_{aw}-H_{hw}\right)
$$

The effective temperature-based hot-side coefficient reported for plotting and diagnostics is:

$$
h_{hot} = \frac{q''_{hot}}{T_{aw}-T_{hw}}
$$

The heat transfer per unit axial length is then:

$$
\frac{d\dot{Q}_{hot}}{dx} = q''_{hot}\frac{dA_{hot}}{dx}
$$

where `dA_hot/dx` comes from:

```python
cooling_circuit.dA_dx_thermal_exhaust(x)
```

---

## 7. Wall conduction

Wall conduction is calculated by:

```python
dQ_cond_dx(x, T_hw, T_cw)
```

The solver supports a stack of wall layers. Each wall layer has local thickness:

$$
\delta_j(x)
$$

and thermal conductivity evaluated at the mean wall temperature:

$$
k_j = k_j\left(\frac{T_{hw}+T_{cw}}{2}\right)
$$

The hot-side thermal area per unit axial length is:

$$
\frac{dA_{hot}}{dx}
$$

The thermal resistance per unit axial length of wall layer `j` is:

$$
R_j(x) = \frac{\delta_j(x)}{k_j\left(\frac{T_{hw}+T_{cw}}{2}\right)\frac{dA_{hot}}{dx}}
$$

The total wall resistance per unit length is the series sum:

$$
R_{tot}(x) = \sum_j R_j(x)
$$

The conduction heat flow per unit axial length is:

$$
\frac{d\dot{Q}_{cond}}{dx} = \frac{T_{hw}-T_{cw}}{R_{tot}}
$$

Interface temperatures are reconstructed after solving. Starting with:

$$
T_0 = T_{hw}
$$

and using:

$$
q'_x = \frac{d\dot{Q}_{cond}}{dx}
$$

for each wall layer:

$$
T_{j+1} = T_j - q'_x R_j
$$

The resulting list is:

$$
[T_{hw},\ T_1,\ T_2,\ \ldots,\ T_{cw}]
$$

The output array reverses this order so that the saved temperature vector at each station is:

$$
[T_{cool},\ T_{cw},\ \ldots,\ T_{hw}]
$$

---

## 8. Coolant-side heat transfer

The coolant-side coefficient is calculated by:

```python
cold_side_coefficients(x, T_cw, T_cool)
```

### 8.1 Coolant film and bulk properties

The film temperature is:

$$
T_{cf} = \frac{T_{cool}+T_{cw}}{2}
$$

The current implementation obtains a pressure for these property calls from the combustion-transport model:

$$
p \leftarrow p_g(x)
$$

and then evaluates:

$$
k_{cf} = k_c(T_{cf},p)
$$

$$
c_{p,cr} = c_{p,c}(T_{cf},p)
$$

$$
\mu_{cf} = \mu_c(T_{cf},p)
$$

The coolant bulk density and viscosity are evaluated at the coolant bulk temperature:

$$
\rho_c = \rho_c(T_{cool},p)
$$

$$
\mu_c = \mu_c(T_{cool},p)
$$

Implementation note: this pressure choice is specific to `cold_side_coefficients`. The coolant temperature and pressure marching functions use the coolant pressure array. If coolant properties are strongly pressure-dependent in a particular case, this distinction should be checked.

### 8.2 Coolant velocity and Reynolds number

The local single-channel flow area is:

$$
A_c(x) = \texttt{cooling\_circuit.A\_coolant(x)}
$$

The local hydraulic diameter is:

$$
D_c(x) = \texttt{cooling\_circuit.Dh\_coolant(x)}
$$

Coolant velocity is computed by `physics.u_coolant`:

$$
u_c = \frac{\dot{m}_{c,chan}}{\rho_c A_c}
$$

The Reynolds number is computed by `physics.reynolds`:

$$
Re_c = \frac{\rho_c u_c D_c}{\mu_c}
$$

### 8.3 Curvature factor

The local channel radius of curvature is:

$$
R_{curv}(x) = \texttt{cooling\_circuit.radius\_of\_curvature(x)}
$$

The curvature factor computed by `physics.phi_curv` is:

$$
\phi_{curv} = \left[Re_c\left(\frac{0.5D_c}{R_{curv}}\right)^2\right]^{0.05}
$$

For a straight section:

$$
R_{curv}=\infty \Rightarrow \phi_{curv}=1
$$

Implementation note: the current `cold_side_coefficients` method computes `phi_curv`, but then calls the Colburn correlation with `phi_curv=1`. Therefore, curvature is reported in the returned dictionary but is not currently applied to the coolant-side heat-transfer coefficient.

### 8.4 Colburn coolant-side heat-transfer coefficient

The coolant-side heat-transfer coefficient is calculated by `physics.h_coolant_colburn`:

$$
h_c = 0.023\frac{k_{cf}}{D_c}
\left(\frac{c_{p,cr}}{k_{cf}\mu_{cf}}\right)^{0.4}
\left(\frac{\dot{m}_{c,chan}D_c}{A_c}\right)^{0.8}
\phi_{curv}
$$

In the current solver call:

$$
\phi_{curv}=1
$$

The result is multiplied by a user-supplied coolant-side correction factor:

$$
h_{c,eff} = h_c C_{cold}
$$

where:

$$
C_{cold} = \texttt{thrust\_chamber.h\_cold\_corr}
$$

### 8.5 Coolant heat pickup per unit length

The coolant-side thermal resistance per unit length is delegated to the cooling-circuit geometry object:

```python
R_coolant_per_len(x, h_c=h_c, T_wall_rep=T_rep)
```

with representative wall/coolant temperature:

$$
T_{rep} = \frac{T_{cw}+T_{cool}}{2}
$$

The heat pickup per unit axial length is:

$$
\frac{d\dot{Q}_{cold}}{dx} = \frac{T_{cw}-T_{cool}}{R_{cool/len}}
$$

This lets the geometry object decide the effective coolant-side area and fin efficiency details, instead of hard-coding them directly in `solver.py`.

---

## 9. Local wall-temperature solve

At each axial station `i`, the solver forms three per-cell heat flows:

$$
Q_{hot,i} = \left(\frac{d\dot{Q}_{hot}}{dx}\right)_i \Delta x
$$

$$
Q_{cond,i} = \left(\frac{d\dot{Q}_{cond}}{dx}\right)_i \Delta x
$$

$$
Q_{cold,i} = \left(\frac{d\dot{Q}_{cold}}{dx}\right)_i \Delta x
$$

The local steady heat balance is:

$$
Q_{hot,i} = Q_{cond,i} = Q_{cold,i}
$$

The nonlinear residual vector is:

$$
R_1 = Q_{hot,i} - Q_{cond,i}
$$

$$
R_2 = Q_{cond,i} - Q_{cold,i}
$$

The residuals are scaled by:

$$
Q_{ref} = \max\left(|Q_{hot,i}|, |Q_{cond,i}|, |Q_{cold,i}|, 1\right)
$$

so the least-squares residual vector is:

$$
\mathbf{r} =
\begin{bmatrix}
R_1/Q_{ref}\\
R_2/Q_{ref}
\end{bmatrix}
$$

The solver uses `least_squares` with:

- trust-region reflective method, `method="trf"`,
- robust `soft_l1` loss,
- bounds enforcing physical wall-temperature ordering,
- `xtol = ftol = gtol = 1e-10`,
- `max_nfev = 200`.

After convergence, the solved values are:

$$
T_{cw,i} = T_{cw,sol}
$$

$$
T_{hw,i} = T_{cw,sol}+\Delta T_{w,sol}
$$

The converged wall temperatures become the initial guesses for the next axial station.

---

## 10. Coolant temperature marching

After the wall temperatures are solved at station `i`, the coolant temperature at station `i+1` is updated using the heat absorbed by the coolant in the current cell.

The continuous energy equation is:

$$
\frac{dT_{cool}}{dx} = \frac{1}{\dot{m}_{c,chan}c_p}\frac{d\dot{Q}_{cold}}{dx}
$$

In the implemented marching step, the solver first computes the per-cell coolant heat pickup:

$$
Q_{cold,i}=\left(\frac{d\dot{Q}_{cold}}{dx}\right)_i \Delta x
$$

and then applies:

$$
\Delta T_{cool,i} = \frac{Q_{cold,i}}{\dot{m}_{c,chan}c_p(T_{cool,i},p_i)}
$$

so:

$$
T_{cool,i+1} = T_{cool,i} + \Delta T_{cool,i}
$$

In the source, the function name and argument name still refer to `dQ_cold_dx`, but in the marching call the quantity passed has already been multiplied by `dx`. The implemented update is therefore a finite-volume energy update using heat per cell.

---

## 11. Coolant pressure marching

The pressure update is calculated by:

```python
coolant_pressure_rate(x, T_cool, p_cool)
```

At the current station, the coolant density is:

$$
\rho_c = \rho_c(T_{cool},p_{cool})
$$

The single-channel velocity is:

$$
u_c = \frac{\dot{m}_{c,chan}}{\rho_c A_c}
$$

The Reynolds number is:

$$
Re_{D_h} = \frac{\rho_c u_c D_h}{\mu_c}
$$

The Darcy friction factor is then computed by `physics.f_darcy`.

### 11.1 Darcy friction factor

The following laminar threshold is used:

$$
Re_{lam}=2300
$$

and turbulent threshold:

$$
Re_{turb}=3500
$$

For laminar flow:

$$
f = \frac{64}{Re_{D_h}}
$$

For turbulent smooth-wall flow, a Petukhov-type expression is used:

$$
f = \left(0.79\ln Re_{D_h} - 1.64\right)^{-2}
$$

If roughness is supplied, the Colebrook-White equation is solved iteratively:

$$
\frac{1}{\sqrt{f}} + 2\log_{10}\left(
\frac{\epsilon(x)}{3.71D_h}
+ \frac{2.51}{Re_{D_h}\sqrt{f}}
\right)=0
$$

where:

$$
\epsilon(x) = \texttt{roughness(x)}
$$

In the transitional regime, the solver linearly blends the laminar and turbulent friction factors:

$$
f = \mathrm{interp}\left(Re_{D_h};\ [2300,3500],\ [f_{lam},f_{turb}]\right)
$$

### 11.2 Stagnation pressure gradient

A geometric path-length factor from is computed:

```python
circuit.ds_dx(x)
```

The stagnation pressure gradient is:

$$
\frac{dp_0}{dx} = -\frac{f}{D_h}\frac{\rho_c u_c^2}{2}\frac{ds}{dx}
$$

An equivalent-length curvature term is also computed:

$$
\frac{dL_{eq}}{dx} = \frac{2KD_h}{\pi R_{curv} f}
$$

and:

$$
\text{curvature factor} = 1 + \frac{dL_{eq}}{dx}
$$

but this factor is not applied in the returned pressure gradient in the current implementation.

### 11.3 Static pressure gradient with area change

The coolant channel area derivative is:

$$
\frac{dA_c}{dx}
$$

The static pressure gradient is computed as:

$$
\frac{dp}{dx} = \frac{dp_0}{dx} - \frac{\rho_c u_c^2}{A_c}\frac{dA_c}{dx}
$$

The explicit pressure updates are:

$$
p_{i+1} = p_i + \left(\frac{dp}{dx}\right)_i\Delta x
$$

$$
p_{0,i+1} = p_{0,i} + \left(\frac{dp_0}{dx}\right)_i\Delta x
$$

After the full march, the solver recomputes static pressure from stagnation pressure and dynamic pressure:

$$
q_i = \frac{1}{2}\rho_i u_i^2
$$

$$
p_i = p_{0,i} - q_i
$$

This corrected static pressure overwrites the previously marched static-pressure array in the returned results.

---

## 12. Derived outputs

After all nodes have been marched, the solver calculates several reporting quantities.

### 12.1 Heat flux

At each station:

$$
Q_{hot,i} = \left(\frac{d\dot{Q}_{hot}}{dx}\right)_i\Delta x
$$

$$
A_{hot,i} = \left(\frac{dA_{hot}}{dx}\right)_i\Delta x
$$

The reported heat flux is:

$$
q''_i = \frac{Q_{hot,i}}{A_{hot,i}}
$$

If `A_hot` is zero, the code stores zero to avoid division by zero.

### 12.2 Coolant velocity

The velocity reported at each node is:

$$
u_i = \frac{\dot{m}_{c,chan}}{\rho_i A_{c,i}}
$$

where density is evaluated using the solved coolant temperature and static pressure.

### 12.3 Stagnation temperature

The stagnation temperature is computed from the static coolant temperature and kinetic-energy term:

$$
T_{0,i} = T_i + \frac{u_i^2}{2c_{p,i}}
$$

where:

$$
c_{p,i}=c_p(T_i,p_i)
$$

### 12.4 Wall-interface temperature array

The returned `T` array has shape:

$$
(n_{nodes},\ 1+n_{walls}+1)
$$

The first column is coolant bulk temperature. The remaining columns are wall-interface temperatures ordered from cold side to hot side:

$$
T[i,:] = [T_{cool},\ T_{cw},\ T_{interface,1},\ldots,\ T_{hw}]
$$

### 12.5 Heat-transfer coefficients and adiabatic-wall temperature

For each node, the solver recomputes and stores:

$$
h_{hot}
$$

$$
h_g
$$

$$
h_{cold}
$$

$$
T_{aw}
$$

These are returned as:

```python
"h_hot"
"h_hot_enthalpy"
"h_cold"
"T_aw_hot"
```

---

## 13. Residual logging

When residual logging is enabled, every local wall-temperature solve records:

$$
(cell,\ iteration,\ R_1,\ R_2)
$$

The residual magnitude is:

$$
|R| = \sqrt{R_1^2+R_2^2}
$$

The global residual history is aggregated iteration-by-iteration. For finite `p`, the implemented norm is:

$$
R_{global,k} = \left(\mathrm{mean}\left(|R|^p\right)\right)^{1/p}
$$

For `p = \infty`, it is:

$$
R_{global,k} = \max(|R|)
$$

The final per-cell residual is the last recorded residual magnitude for each cell.

---


## 14. Assumptions and implementation details

The current regenerative cooling solver is an engineering heat-exchanger model. Its main assumptions are:

- The solution is steady in time.
- The solver marches along one spatial coordinate.
- Axial conduction in the wall is neglected.
- Wall conduction is treated as one-dimensional through the wall stack.
- Each coolant channel in the selected circuit receives an equal share of total circuit mass flow.
- Hot-gas heat transfer is modeled by an enthalpy-driven Bartz-style correlation.
- Coolant-side heat transfer is modeled by a Colburn-style turbulent internal-flow correlation.
- Coolant pressure drop is modeled with Darcy friction and a separate area-change correction.
- Hot-gas properties are supplied by equilibrium and TP CEA tables, then interpolated.

There are also a few important implementation details to keep in mind:

- `H_gr` is calculated in the hot-side model, but the current gas properties used in the Bartz equation are still queried at `x`, not at the computed reference enthalpy.
- `phi_curv` is computed in the coolant-side model, but the current Colburn call passes `phi_curv=1`, so curvature does not currently modify `h_cold`.
- `coolant_pressure_rate` computes an equivalent-length curvature factor, but does not apply it to the returned pressure gradients.
- `compute_aerothermodynamics(contour, Nt=...)` currently overwrites `Nt` with `len(self.x_nodes)`.
- `get_p(x, T=..., h=...)` always returns the equilibrium-column pressure and ignores the supplied temperature or enthalpy.
- In `cold_side_coefficients`, coolant properties are evaluated using `combustion_transport.get_p(x)` rather than the marched coolant pressure. Other coolant calculations use the coolant pressure array.

These are limitations in the current implementation that is being worked on for future versions. Nevertheless, the exact behavior is reported here so results using the program are interpreted correctly. 

---

## 15. Returned data structure

`steady_heating_analysis` returns the dictionary produced by `solve_heat_exchanger`:

```python
{
    "x": x_domain,
    "T": T_full,
    "T_static": T_cool_arr,
    "T_stagnation": T_stagnation_arr,
    "p_static": p_static_arr,
    "p_stagnation": p_stagnation_arr,
    "dQ_dA": dQ_dA_arr,
    "velocity": velocity_arr,
    "h_hot": h_hot_arr,
    "h_hot_enthalpy": h_hot_enthalpy_arr,
    "h_cold": h_cold_arr,
    "T_aw_hot": T_aw_hot_arr,
    "residuals": (global_R, final_R),
}
```
---

## 16. Summary

The regenerative cooling solver couples a CEA-based hot-gas property model to a one-dimensional cooling-channel heat-exchanger calculation. The solver walks along the selected cooling circuit, solves a local nonlinear balance for the hot- and cold-side wall temperatures, and then marches coolant temperature and pressure downstream.

At each station, the central balance is:

$$
Q_{hot}=Q_{cond}=Q_{cold}
$$

with:

$$
Q_{hot}=q''_{hot}A_{hot}
$$

$$
Q_{cond}=\frac{T_{hw}-T_{cw}}{R_{wall}}
$$

$$
Q_{cold}=\frac{T_{cw}-T_{cool}}{R_{coolant}}
$$

The hot side is enthalpy-driven, the wall is modeled as a stack of thermal resistances, and the coolant side is temperature-driven. This makes the solver fast enough for design iteration while retaining the main physical couplings needed for regenerative cooling analysis.
