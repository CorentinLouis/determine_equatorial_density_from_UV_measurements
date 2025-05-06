# Determining the Jovian equatorial density from UV measurements of the Galilean moons' footprints

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/876014209.svg)](https://zenodo.org/badge/latestdoi/876014209)


## Overview
This pipeline allows the determination of the density $\rho$ and the height scale $H$ of the Jovian equatorial density from the Trans-hemispheric Electron Beam (TEB) and Main Alfvén Wing (MAW) UV spots and the longitudinal separation between them.

The main steps of the pipeline are the following:

1. We first assume the Plasma Sheet (PS) is aligned with one of the equators (jovicentric, centrifugal, or magnetic). We then assume a plasma sheet height scale  $H$, which gives the magnetic field line length values $L$ in the Northern ($L_N$) and Southern ($L_S$) hemispheres, inside and outside of the PS: $L_\text{N in PS}$, $L_\text{N out PS}$, $L_\text{S in PS}$ and $L_\text{S out PS}$.

2. We assume a plasma density  $\rho_0$  at the center of the PS and determine $\rho$, the density along the magnetic field line in the PS, using the following density profile equation:  

   $\rho_i = \rho_0 exp \left( -\sqrt{\frac{(r_i - r_0)^2 + z_i^2}{H}} \right)$
   
   where $r_0$  is the equatorial radius set at the moon's orbital distance, $r_i = x_i^2+ y_i^2$ is the equatorial radial distance, and  $z_i$  is the altitude above the equator of the position of the measurement point $i$.

3. We determine the Alfvén speed $v_A$, based on the calculated  $\rho_i$  and the magnetic field amplitude $B_i$, using a magnetic field model. From that, we obtain the values of $t_\text{TEB}$, $t_\text{MAW}$, $\lambda_\text{TEB}$, $\lambda_\text{MAW}$, and therefore $\Delta \lambda_\text{calculated}$.

4. By applying this method to different values of $\rho_0$, we minimize $|\Delta \lambda_\text{calculated} - \Delta \lambda_\text{observed}|$ .

5. We then repeat the calculations for different values of the scale height $H$ to minimize $|\lambda_\text{MAW calculated} - \lambda_\text{MAW observed}|$  and $|\lambda_\text{TEB calculated} - \lambda_{TEB observed}|$.

## How to use it

The main function to be used is ```density_and_scale_height_of_CS_calculation```.

The script requires specifying key parameters such as:
  - a magnetic field model that gives the magnetic field values along the magnetic field line connected to the Galilean moon, in the North (`r_N`, `theta_N`, `phi_N`, `B_total_N`) and South (`r_S`, `theta_S`, `phi_S`, `B_total_S`). Distance values need to be in Jovian radii. Magnetic field values need to be in Tesla;
  - an array of electron energies to be tested (`E_array`);
  - specific scale height values $H$ for the Current Sheet to be tested (`z_lim_N_array`);
  - the distance along the chosen equator that will be taken as reference (`r_0`, in planetary Radius)
  - longitudes of the observed MAW and TEB spots (`longitude_observed_MAW`, `longitude_observed_TEB`);
  - a first guess of the value of $\rho_0$ to help the code obtain the correct values for the longitudes of the spots and their separation (`rho_0_volumetric_first_guess`; default is $3.6436e-21 \text{kg/m}^3$, i.e., $4000 \text{particles/cm}^3$);
  - the rotation rate of the moon (`rotation_rate`; default is $0.0098$, Callisto's rotation rate);
  - the equator of reference (`jovicentric_equator`, `centrifugal_equator`, or `magnetic_equator`; default is `centrifugal_equator = True`);
  - density model can be chosen between `disk` and `torus` (default is `torus = True`);
  - Uncertainties can be added using the $\Delta L$ (`delta_L`) input. These uncertainties are based on the tilt angle of the Alfvén wings, which increases the length of the magnetic field line in the plasma sheet as $L = L_\text{in PS} \times \Delta L$, with $\Delta L= \frac{1}{cos \theta}, with $\theta = atan(M_A)$, with $M_A$ the Mach Alfvén number.

  - verbose mode can be deactivated using `verbose = False`.

```python
(H_array,                  # Height Scale in Rp
 rho_0_array,              # Total mass density at the center of the PS in kg/cm^3
 n_0_ions_array,           # Ion density in the PS in cm^{-3}
 n_0_electrons,            # Electron density in the PS in cm^{-3}
 longitude_TEB,            # TEB longitude
 longitude_MAW,            # MAW longitude
 rho_at_Moon_array,        # Total mass density at the moon
 n_at_Moon_ions,           # Ion density at the moon's orbit in cm^{-3}
 n_at_Moon_electrons,      # Electron density at the moon's orbit in cm^{-3}
 units
) = density_and_scale_height_of_CS_calculation(
        r_N, theta_N, phi_N, B_total_N,
        r_S, theta_S, phi_S, B_total_S,
        E_array,
        z_lim_N_array, r_0,
        longitude_observed_MAW,
        longitude_observed_TEB,
        rho_0_volumetric_first_guess=rho_0_volumetric_first_guess,
        rotation_rate=rotation_rate,
        jovicentric_equator=False,
        centrifugal_equator=True,
        magnetic_equator=False,
        verbose=True,
        disk=False,
        torus=True,
        delta_L=delta_L
)
```


## Example of Use: Reproducing Rabia et al. (2025) Results

The `Rabia_2025_Nature.ipynb` file is a Jupyter notebook that provides an example of how to reproduce the results in Jonas et al. (2025).
In this notebook, two methods are provided:
1. The first one takes different values of $H$ as inputs and provides the best $\rho$ results for each $H$ that minimize the TEB and MAW longitudes, and the longitudinal separation between them.
2. The second one minimizes the longitudes of the TEB and MAW spots, as well as their separation, across different $H$ and $\rho$ values to find the best fit.

The first method reproduces the results shown in Table 2 of Rabia et al. (2025) for $H=0.75$, $1.00$, $1.50$, $2.00$, and $3.00$ $R_\text{J}$​, while the second method yields the best results highlighted in red in the table for $H=0.94 \text{R}_\text{J}$​.

