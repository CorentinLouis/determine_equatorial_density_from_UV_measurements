import numpy

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle
from density import rho_z

def plot_MFL_and_rho_CS(x_N, y_N, z_N,
                        x_S, y_S, z_S, 
                        mask_CS_N, mask_CS_S,
                        r_0, H, rho_0, disk = False,
                        cb_title = r'$\rho$ (kg.cm$^{-3}$)',
                        ):



    equatorial_vals = numpy.concatenate([numpy.sqrt(x_N[mask_CS_N]**2 + y_N[mask_CS_N]**2), numpy.sqrt(x_S[mask_CS_S]**2 + y_S[mask_CS_S]**2)])
    z_vals_CS = numpy.concatenate([z_N[mask_CS_N], z_S[mask_CS_S]])
    
    
    plt.figure(figsize=(8, 6))

    plt.plot(numpy.sqrt(x_N**2+y_N**2), z_N, 'k')
    plt.plot(numpy.sqrt(x_S**2+y_S**2), z_S, 'k', linestyle = '--')
    plt.axhline(y=0, xmin = 0, xmax = 25, color = 'k')



    rho_CS_N = rho_z(x_N[mask_CS_N], y_N[mask_CS_N], z_N[mask_CS_N], r_0, H, rho_0, disk=disk)
    rho_CS_S = rho_z(x_S[mask_CS_S], y_S[mask_CS_S], z_S[mask_CS_S], r_0, H, rho_0, disk=disk)
    rho_CS_vals = numpy.concatenate([rho_CS_N, rho_CS_S])
    
    
    # Northern points with color corresponding to rho_N_CS
    sc1 = plt.scatter(equatorial_vals, z_vals_CS, c=rho_CS_vals, cmap='viridis', marker='o',  norm=LogNorm())


    # Add colorbar for rho_N_CS and rho_S_CS
    plt.colorbar(sc1, label=cb_title)

    circle1 = Circle((0, 0), 1, color='black', fill=True, linewidth=2)
    plt.gca().add_patch(circle1)
    plt.xlim(0,28)
    plt.ylim(-14,14)

    circle2 = Circle((numpy.sqrt(x_N[0]**2+y_N[0]**2),z_N[0]),.2, color='black', fill=True, linewidth=2)
    plt.gca().add_patch(circle2)

    # Label axes
    plt.xlabel(r'equatorial radius ($R_J)$')
    plt.ylabel(r'z ($R_J)$')