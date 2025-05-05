# Determining the Jovian equatorial density from UV measurements of the Galilean moons' footprint

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/876014209.svg)](https://zenodo.org/badge/latestdoi/876014209)


## Overview
This pipeline allows to determine the density $\rho$ and the height scale $H$ of the Jovian equatorial density, from the Trans-hemispheric Electron Beam (TEB) and Main Alfvén Wing (MAW) UV spots and the separation in longitude separation between them.

The main steps of the pipeline are the following:

1. We first assume the Plasma Sheet (LS) to be aligned with one of the equator (jovicentric, centrifugal or magnetic). We then assume a plasma sheet height scale $H$, which gives the values magnetic field line length value $L$ in the Northern ($L_N$) and in the Southern ($L_S$) hemispheres inside and outside of the PS: $L_\text{N in PS}$, $L_\text{N out PS}$, $L_\text{S in PS}$ and $L_\text{S out PS}$

2. We assume a plasma density $\rho_0$  at the center of the PS and we determine $\rho$ the density along the magnetic field line in the PS using the following density profile equation:
$\rho_i = \rho_0 exp \left( -\sqrt{\frac{(r_i - r_0)^2 + z_i^2}{H}} \right)$
with $r_0$ the equatorial diameter set at the moons' orbital distance, $r_i = x_i^2+ y_i^2$ the equatorial radial distance, and $z_i$ the altitude above the equator of the position of the measurement point i. 

4. We determine the Alfvén speed velocity $v_A$, based on the calculated $\rho_i$ and the magnetic field amplitude $B_i$ using a magnetic field model. From that, we obtain the values of $t_\text{TEB}$, $t_\text{MAW}$, $\lambda_\text{TEB}$, $\lambda_\text{MAW}$ and therefore $\Delta \lambda_\text{calculated}$.

5. By applying this method on different values of $\rho_0$, we minimize $|\Delta \lambda_\text{calculated} - \Delta \lambda_\text{observed}|$ .

6. We then run the same above calculations for different values of the scale height $H$ to minimize $|\lambda_\text{MAW calculated} - \lambda_\text{MAW observed}|$  and $|\lambda_\text{TEB calculated} - \lambda_{TEB observed}|$ .




## Exemple of Use
The script requires specifying key parameters such as:
  - a magnetic field model that gives the magnetic field value along the magnetic field line in the North (r_N, theta_N, phi_N, B_total_N) and South (r_S, theta_S, phi_S, B_total_S). Distance values need to be in Jovian radii. Magnetic field values need to be in Tesla;
  - an array of electron energy to be tested (E_array);
  - specific scale height values H for the Current Sheet to be tested (z_lim_N_optimized);
  - longitude of the observed MAW and TEB spots (longitude_observed_MAW,longitude_observed_TEB).
  - a first guess of the value of $\rho_0$ to help the code to obtain the right values for the longitude of the spots and the separation in longitude beteween them (rho_0_volumetric_first_guess; default is 3.6436e-21 kg/m^3, i.e., 4000 particules/cm^3);
  - the rotation rate of the moon (default is 0.0098, Callisto's rotation rate;
  - the equator of reference (jovicentric_equator, centrifugal_equator or magnetic_equator; default is centrifugal_equator = True);
  - density model can be choosen between disk and torus (default is torus = True).
  - - verbose mode can be de-activated using verbose = False;

Rabia_2024_Nature.ipynb file is a jupyter notebook that give an example to reproduce the results in Jonas et al. (2025).
In this Jupyter notebook, two methods are given:
  - The first one takes different values of H as inputs and give the best $\rho$ results for these different H that minimise the longitudes values of the TEB, MAW, and the longitude separation between them
  - The second one minimizes the longitude values of TEB and MAW spots and the longitude separation between them for different H and $\rho$ values and give the best result.

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
