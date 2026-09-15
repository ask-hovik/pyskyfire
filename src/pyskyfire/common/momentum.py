"""One-dimensional static-pressure changes shared by heated duct solvers."""


def acceleration_pressure_drop(rho_up, velocity_up, rho_down, velocity_down):
    """Trapezoidal integral of rho*u*du [Pa]; positive for acceleration.

    Exact for constant density (Bernoulli) and constant area (mass-flux form).
    This is reversible static-pressure change, separate from wall friction.
    """
    return (
        0.5
        * (rho_up * velocity_up + rho_down * velocity_down)
        * (velocity_down - velocity_up)
    )
