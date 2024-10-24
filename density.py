from numpy import exp, sqrt

# Defined the density function
def rho_z(x_i, y_i, z_i, r_0, H, rho_0, torus = True, disk = False, n=1):
    
    r_i = sqrt(x_i**2 + y_i**2)
    if disk == True:
        density = rho_0 * exp(-(r_i/r_0)**n)*exp(-(z_i/H)**n)
        torus == False
    if torus == True:
        density = rho_0 * exp(- (sqrt((r_i - r_0)**2 + z_i**2) / H)**n)
    
    return density


