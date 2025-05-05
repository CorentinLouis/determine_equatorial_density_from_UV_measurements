# Determining the Jovian equatorial density from UV measurements of the Galilean moons' footprint

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/876014209.svg)](https://zenodo.org/badge/latestdoi/876014209)


## Overview
This pipeline allows to determine the density $\rho$ and the height scale $H$ of the Jovian equatorial density, from the Trans-hemispheric Electron Beam (TEB) and Main Alfvén Wing (MAW) UV spots and the separation in longitude separation between them.


## Exemple of Use
The script requires specifying key parameters such as:
  -* a magnetic field model that gives the magnetic field value along the magnetic field line in the North (r_N, theta_N, phi_N, B_total_N) and South (r_S, theta_S, phi_S, B_total_S);
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
