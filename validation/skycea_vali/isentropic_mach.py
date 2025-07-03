def M_from_A_subsonic(A, A_t, gamma):
    """Get the Mach number from the local flow area, assuming subsonic flow.

    Args:
        A (float): Local area (m2)
        A_t (float): Throat area (m2)
        gamma (float): Ratio of specific heats cp/cv
    """

    def func_to_solve(Mach):
        return  A/A_t - A_At(M = Mach, gamma = gamma)
        
    return scipy.optimize.root_scalar(func_to_solve, bracket = [1e-10, 1], x0 = 0.5).root

def M_from_A_supersonic(A, A_t, gamma):
    """Get the Mach number from the local flow area, assuming supersonic flow.

    Args:
        A (float): Local area (m2)
        A_t (float): Throat area (m2)
        gamma (float): Ratio of specific heats cp/cv
    """
    def func_to_solve(Mach):
        return  A/A_t - A_At(M = Mach, gamma = gamma)

    return scipy.optimize.root_scalar(func_to_solve, bracket = [1, 500], x0 = 1).root

def A_At(M, gamma):
    """Ratio of local flow area to the throat area, for a given Mach number [1].

    Args:
        M (float): Mach number
        gamma (float): Ratio of specific heats cp/cv
    """
    return 1/M * ( (2 + (gamma-1) * M**2) / (gamma + 1) )**( (gamma+1) / (2 * (gamma-1) ) )
