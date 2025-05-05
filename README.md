# Determining the Jovian equatorial density from UV measurements of the Galilean moons' footprint

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/876014209.svg)](https://zenodo.org/badge/latestdoi/876014209)


## Overview
This pipeline allows to determine the density $\rho$ and the height scale $H$ of the Jovian equatorial density, from the Trans-hemispheric Electron Beam (TEB) and Main Alfvén Wing (MAW) UV spots and the separation in longitude separation between them.


## Exemple of Use
The script requires specifying key parameters such as:
  -* a magnetic field model that gives the magnetic field value along the magnetic field line in the North (r_N, theta_N, phi_N, B_total_N) and South (r_S, theta_S, phi_S, B_total_S). Distance values need to be in Jovian radii. Magnetic field values need to be in Tesla;
  -* an array of electron energy to be tested (E_array);
  -* specific scale height values H for the Current Sheet to be tested (z_lim_N_optimized);
  -* longitude of the observed MAW and TEB spots (longitude_observed_MAW,longitude_observed_TEB).
  
Some optionnal parameters can be added, such as:
  -* a first guess of the value of $\rho_0$ to help the code to obtain the right values for the longitude of the spots and the separation in longitude beteween them (rho_0_volumetric_first_guess; default is 3.6436e-21 kg/m^3, i.e., 4000 particules/cm^3);
  -* the rotation rate of the moon (default is 0.0098, Callisto's rotation rate;
  -* the equator of reference (jovicentric_equator, centrifugal_equator or magnetic_equator; default is centrifugal_equator = True);
  -* verbose mode can be de-activated using verbose = False;
  -* density model can be choosen between disk and torus (default is torus = True).

Rabia_2024_Nature.ipynb file is a jupyter notebook that give an example to reproduce the results in Jonas et al. (2025).
In this Jupyter notebook, two methods are given:
  -* The first one takes different values of H as inputs and give the best $\rho$ results for these different H that minimise the longitudes values of the TEB, MAW, and the longitude separation between them
  -* The second one minimizes the longitude values of TEB and MAW spots and the longitude separation between them for different H and $\rho$ values and give the best result.

The first method gives the result showing in Table 2 in Rabia et al. (2025) for values H = 0.75, 1.00, 1.50, 2.00 and 3.00 $R_J$ while the second method gives the best results highlited in red in Table to for H = 0.94 $R_J$.

## Application
We give here an example for Callisto, but this method can be applied to any other Galilean moons, as long as you have a longitude value for the TEB and MAW spots, and therefore the separation between them.
A second version of this code will contain the same method for different pair of spots, including the first Reflected Alfvén Wing (RAW) spot.

A general approach will be the following (with the required key parameters defined by users):

```python
# Outer minimization function to find the optimal z_lim_N
result = minimize(
    objective_function_zlimN_from_longitude,
    z_lim_N_first_guess,
    args=(r_0, rho_0_volumetric_first_guess,
          r_N, theta_N, phi_N, B_N,
          r_S, theta_S, phi_S, B_S,
          E, rotation_rate_callisto,
          360-phi_N[0], longitude_observed_MAW, longitude_observed_TEB, 
          10*71492e3, # equatorial_radius_lim
          71492e3, # Rp
          MAW, TEB, # MAW=True, TEB=False
          disk, torus,  # disk=False, torus=True
          jovicentric_equator, centrifugal_equator, magnetic_equator,
          delta_L),  # jovicentric_equator=False, centrifugal_equator=True, magnetic_equator=False
    method='Nelder-Mead',
    tol=1e-6
)

z_lim_N_optimized = result.x


(H_array, # Height Scale in Rp
            rho_0_array, #total mass density in the center of the PS in kg.cm^-3
            n_0_ions_array, # ions density in the PS cm^-3
            n_0_electrons, # electrons density in the PS cm^-3
            longitude_TEB, # TEB longitude
            longitude_MAW, # MAW longitude
            rho_at_Moon_array, # total mass density at the moon
            n_at_Moon_ions, # ions density at the moon orbit cm^-3
            n_at_Moon_electrons, # electrons density at the moon orbit cm^-3 
            units
            ) = density_and_scale_height_of_CS_calculation(r_N, theta_N, phi_N, B_total_N,
                                   r_S, theta_S, phi_S, B_total_S,
                                   E_array,
                                   z_lim_N_optimized, r_0,
                                   longitude_observed_MAW,
                                   longitude_observed_TEB,
                                   rho_0_volumetric_first_guess = rho_0_volumetric_first_guess,
                                   rotation_rate = rotation_rate_callisto,
                                   jovicentric_equator = False,
                                   centrifugal_equator = True,
                                   magnetic_equator = False,
                                   verbose = True,
                                   disk = False,
                                   torus = True,
                                   delta_L = delta_L)

```
