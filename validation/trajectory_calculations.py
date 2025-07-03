import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Constants
g0 = 9.81
Isp = 311
T = 845e3 * 9
Cd = 0.4
d = 6  # Rocket diameter
Sref = np.pi * (d/2)**2
pe = 0.8e5  # Exit pressure
mdot = T / (g0 * Isp)
Ae = Sref / 2  # Exit area for rocket engine
p0 = 101325  # Pascal, sea level pressure
h_scale = 7640  # Scale height
Rp = 6387e3  # Earth radius
m0 = 566600
rho0 = 1.225

# Functions
def find_m(t):
    return m0 - mdot * t

def find_D(rho, v):
    return 0.5 * rho * v**2 * Cd * Sref

def find_T(pinf):
    return mdot * Isp * g0 + (pe - pinf) * Ae

def find_pinf(h):
    return p0 * np.exp(-h / h_scale)

def find_rho(h):
    return rho0 * np.exp(-h / h_scale)

def find_g(h):
    return g0 / (1 + h / Rp)**2

def find_h(v, gamma, t):
    return v * t * np.sin(gamma)

def find_x(dxdt, h, t):
    return dxdt/(1 + h/Rp)*t

def find_gamma(t):
    end_time = 300
    if t < end_time:
        gamma = np.pi/2 - np.pi/2 * t / end_time
    else:
        gamma = 0
    return gamma

# Differential equations
def derivatives(t, y):
    v = y[0]
    gamma = find_gamma(t)
    h = find_h(v, gamma, t)
    pinf = find_pinf(h)
    rho = find_rho(h)
    T = find_T(pinf)
    D = find_D(rho, v)
    m = find_m(t)
    g = find_g(h)
    
    dvdt = (T - D) / m - g * np.sin(gamma)
    return [dvdt]

# Time points and storage for parameters
t_span = (0, 225)
t_eval = np.linspace(0, 225, num=50)

# Initial condition
v0 = 1e-5
y0 = [v0]

# Solve ODE
solution = solve_ivp(derivatives, t_span, y0, t_eval=t_eval, method='RK45')

# Extract results
v = solution.y.flatten()

# Use results to update functions
gamma = np.array([find_gamma(t) for t in t_eval])
h = np.array([find_h(vi, gi, ti) for vi, gi, ti in zip(v, gamma, t_eval)])
dxdt = v * np.cos(gamma)
x = np.array([find_x(dxi, hi, ti) for dxi, hi, ti in zip(dxdt, h, t_eval)])
pinf = find_pinf(h)
rho = find_rho(h)
T = find_T(pinf)
D = find_D(rho, v)
m = find_m(t_eval)
g = find_g(h)

# Begin plot setup
plt.figure(figsize=(18, 24))

# Plot altitude
plt.subplot(5, 2, 1)
plt.plot(t_eval, h)
plt.title('Altitude')

# Plot pressure at infinity
plt.subplot(5, 2, 2)
plt.plot(t_eval, pinf)
plt.title('Pressure at Infinity')

# Plot density
plt.subplot(5, 2, 3)
plt.plot(t_eval, rho)
plt.title('Density')

# Plot thrust
plt.subplot(5, 2, 4)
plt.plot(t_eval, T)
plt.title('Thrust')

# Plot drag
plt.subplot(5, 2, 5)
plt.plot(t_eval, D)
plt.title('Drag')

# Plot mass
plt.subplot(5, 2, 6)
plt.plot(t_eval, m)
plt.title('Mass')

# Plot gravity
plt.subplot(5, 2, 7)
plt.plot(t_eval, g)
plt.title('Gravity')

# Plot flight path angle
plt.subplot(5, 2, 8)
plt.plot(t_eval, np.degrees(gamma))  # Convert radians to degrees for plotting
plt.title('Flight Path Angle')

# Plot velocity
plt.subplot(5, 2, 9)
plt.plot(t_eval, v)  
plt.title('Velocity')

# Plot downrange position
plt.subplot(5, 2, 10)
plt.plot(t_eval, x)
plt.title('Position')

plt.tight_layout()
plt.show()

# trajectory visualisation
plt.figure(figsize=(10, 6))  # Set the size of the figure (width, height in inches)
plt.plot(x, h, label='Rocket Trajectory')  
plt.xlabel('Downrange Position (meters)')  
plt.ylabel('Altitude (meters)')  
plt.title('Rocket Trajectory: Altitude vs. Downrange Position')  
plt.grid(True)  
plt.legend()  
plt.show()
