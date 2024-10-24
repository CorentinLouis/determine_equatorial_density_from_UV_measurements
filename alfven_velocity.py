from numpy import sqrt

# Defined the Alfven speed
def v_alfven_calculation(B, rho):
    mu_0 = 1.2566e-6  # Vacuum permittivity [kg.m.A^-2.s^-2]
    v_A = sqrt(B ** 2 / (mu_0 * rho))

    return v_A
