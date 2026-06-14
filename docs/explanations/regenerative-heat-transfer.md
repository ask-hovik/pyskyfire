# Regenerative Cooling Heat Transfer in Pyskyfire

The regenerative cooling solver is one of the core capabilities of pyskyfire. It predicts the heat transfer from the hot combustion gases through the thrust chamber wall into the coolant flowing in the cooling channels. The solver is implemented as a quasi-1D, steady-state model that balances heat transfer on the hot gas side, conduction through the wall (including multiple layers), and heat transfer into the coolant.

This explanation focuses on the physics and assumptions implemented in the `regen` module (primarily `physics.py` and `solver.py`).

## 1. Overall Approach and Assumptions

The regenerative cooling model follows a **quasi-1D steady-state** approach, consistent with established tools such as those described by Binder et al. (1997). The main assumptions are:

- Flow and heat transfer are treated as one-dimensional along the axial direction of the thrust chamber.
- No axial heat conduction in the wall (temperature gradients in the radial direction dominate).
- Propellants are treated as incompressible while in the liquid phase.
- Chemical equilibrium is assumed at every station in the combustion gas.
- Heat transfer coefficients are calculated using established engineering correlations (Bartz on the hot gas side, Colburn on the coolant side) rather than full CFD.
- The model is enthalpy-driven on the hot gas side (better correlation with experimental data for combustion flows) and temperature-driven on the coolant side.

The solver marches along the thrust chamber contour, solving a system of differential equations for wall temperatures, coolant temperature, and coolant pressure at each station.

## 2. Combustion Gas Aerothermodynamics (Hot Gas Properties)

Before solving the cooling problem, the hot gas properties along the thrust chamber must be known. These are obtained from chemical equilibrium calculations (via the `skycea` module, based on NASA CEA-style methods).

At each axial station $x$:

- The local area ratio $\varepsilon(x) = A_x / A_t$ is calculated from the contour geometry.
- The local Mach number $M(x)$ is found from the isentropic area-Mach relation (subsonic before the throat, supersonic after).
- Local temperature $T(x)$ and pressure $p(x)$ are computed using isentropic relations.
- Gas composition is equilibrated at the local $T$ and $p$.
- Thermodynamic and transport properties are then interpolated:

$$
\rho_g, \mu_g, \Pr_g, c_{p,g}, k_g = f(T, p)
$$

The isentropic exponent is:

$$
\gamma = \frac{c_p}{c_v}
$$

These properties (especially those evaluated at reference enthalpy conditions averaged between the freestream gas and wall) are used as inputs to the hot gas heat transfer calculation.

## 3. Hot Gas Side Heat Transfer

pyskyfire uses an **enthalpy-driven** heat transfer formulation on the hot gas side, following Binder et al. (1997). This approach generally shows better agreement with experimental data than purely temperature-driven models for combustion gases.

The heat transfer rate to the wall is:

$$
\dot{Q}_{hw} = h_g A_h (H_{aw} - H_{hw})
$$

where:
- $H_{aw}$ = adiabatic wall enthalpy
- $H_{hw}$ = hot wall enthalpy
- $A_h$ = heat transfer area on the hot gas side

The enthalpy-based heat transfer coefficient is defined as:

$$
h_g = \frac{h_{gr}}{c_{p,gr}}
$$

Here the subscript $r$ denotes properties evaluated at a reference enthalpy condition (average between freestream gas and wall surface conditions).

The hot gas heat transfer coefficient $h_{gr}$ is calculated using the **Bartz correlation** (as presented by Binder):

$$
h_{\text{Bartz}} = 0.026 \frac{k_{gr}}{D_h} \left( \frac{c_{p,gr}}{k_{gr} \mu_{gr}} \right)^{0.4} \left( \frac{\dot{m}_g D_h}{A_{\text{chmb}}} \right)^{0.8} \frac{T_g}{T_{gr}}
$$

A curvature correction is sometimes applied in the code for regions of high contour curvature.

## 4. Conduction Through the Wall

Heat conduction through the thrust chamber wall (which may consist of multiple layers, e.g., liner + thermal barrier coating) is modeled using Fourier's law in one dimension (radial direction).

For a wall layer of thickness $\delta_j$ and conductivity $k_j$:

$$
R_j(x) = \frac{\delta_j(x)}{k_j} \frac{dA_h}{dx}
$$

The temperature drop across each layer is:

$$
\Delta T_j = \frac{d\dot{Q}_{\text{cond}}}{dx} R_j
$$

The temperature at the interface between layers is updated iteratively:

$$
T_{j+1} = T_j - \Delta T_j
$$

This allows pyskyfire to compute both the hot-side wall temperature ($T_{hw}$) and cold-side wall temperature ($T_{cw}$) even with multi-layer walls.

## 5. Cold Side (Coolant) Heat Transfer

On the coolant side, a **temperature-driven** formulation is used (Newton's law of cooling):

$$
\dot{Q}_{cw} = h_c A_c (T_{cw} - T_{\text{cool}})
$$

The coolant-side heat transfer coefficient is calculated using a modified **Colburn correlation** with a curvature correction factor:

$$
h_{\text{Colburn}} = 0.023 \frac{k_{cf}}{D_{h,c}} \left( \frac{c_{p,cf}}{k_{cf} \mu_{cf}} \right)^{0.4} \left( \frac{\dot{m}_c D_{h,c}}{A_{\text{chn}}} \right)^{0.8} \varphi_{\text{curv}}
$$

where the subscript $cf$ denotes coolant properties evaluated at film conditions (average of bulk coolant and cold wall temperature).

The **curvature correction factor** $\varphi_{\text{curv}}$ accounts for the effect of channel curvature on the boundary layer:

$$
\varphi_{\text{curv}} = \left( \text{Re}_{D_h} \frac{D_h}{2 R_{\text{curv}}} \right)^{0.05}
$$

The local radius of curvature $R_{\text{curv}}(x)$ is computed from the second derivative of the contour:

$$
R_{\text{curv}}(x) = \frac{\left(1 + \left( \frac{dr}{dx} \right)^2 \right)^{3/2}}{\frac{d^2 r}{dx^2}}
$$

## 6. Overall Solver and Coupling (solver.py)

The three heat transfer mechanisms (hot gas, conduction, coolant) are coupled through energy balance at steady state:

$$
\dot{Q}_{hw} = \dot{Q}_{\text{cond}} = \dot{Q}_{cw}
$$

The solver marches along the thrust chamber and solves the following system of differential equations:

**Hot gas side:**
$$
\frac{d\dot{Q}_{hw}}{dx} = h_g(x) \frac{dA_h}{dx} (H_{aw}(x) - H_{hw}(x))
$$

**Wall conduction:**
$$
\frac{d\dot{Q}_{\text{cond}}}{dx} = k(x) \frac{dA_h}{dx} \frac{T_{hw}(x) - T_{cw}(x)}{L(x)}
$$

**Coolant side:**
$$
\frac{d\dot{Q}_{cw}}{dx} = h_c(x) \frac{dA_c}{dx} (T_{cw}(x) - T_{\text{cool}}(x))
$$

From energy balance on the coolant:

$$
\frac{dT_{\text{cool}}}{dx} = \frac{1}{\dot{m} c_p} \frac{d\dot{Q}_{cw}}{dx}
$$

Coolant pressure drop (including friction):

$$
\frac{dp_{0,\text{cool}}}{dx} = -\frac{f \rho v^2}{D_h \cdot 2}
$$

The friction factor $f$ is determined from the Reynolds number using the laminar Darcy formula or the Colebrook-White equation (with linear blending in the transitional regime $2300 < \text{Re}_{D_h} < 3500$).

The solver iterates until the residuals (typically update residuals between iterations) converge to a satisfactory tolerance. pyskyfire tracks both local cell residuals and a global norm (usually $L_\infty$ or averaged $L_2$).

## Summary of Key Physics

- **Hot gas side**: Enthalpy-driven (Binder/Bartz)
- **Wall**: Multi-layer conduction with temperature continuity
- **Coolant side**: Temperature-driven (Colburn) + curvature correction
- **Coupling**: Strict energy balance at steady state
- **Numerical method**: Forward marching 1D finite volume / finite difference approach

This model strikes a practical balance between physical fidelity and computational speed, making it suitable for design space exploration while still producing results that compare well with RL10 data (as shown in the validation section of the thesis).

---

**References** (from the original thesis)
- Binder et al. (1997) — Enthalpy-based heat transfer approach
- Huzel & Huang — Classic rocket propulsion design reference
- Bartz and Colburn correlations (standard engineering practice)T