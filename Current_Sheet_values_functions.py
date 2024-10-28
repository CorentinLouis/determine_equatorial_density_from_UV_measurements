import numpy
from scipy.optimize import minimize

from coordinates_system_transformation import spherical_to_cartesian, cartesian_to_spherical

from ion_populations import single_ion_population_mass
from longitude_functions import delta_longitude_calculation_no_CS, objective_function_rho_0_from_delta_longitude, longitude_MAW_TEB_calculation
from density import rho_z

def density_and_scale_height_of_CS_calculation(r_N, theta_N, phi_N, B_total_N,
                                   r_S, theta_S, phi_S, B_total_S,
                                   E_array,
                                   z_lim_N_array, r_0, 
                                   lambda_observed_MAW,
                                   lambda_observed_TEB,
                                   jovicentric_equator = False,
                                   centrifugal_equator = False,
                                   magnetic_equator = False,
                                   verbose = True,
                                   disk = False,
                                   torus = False,
                                   rho_0_volumetric_first_guess = 3.6436e-21, #kg/m^3 --> 4000 particules/cm^3
                                   rotation_rate =  0.0098,
                                   equatorial_radius_lim = 10*71492e3,
                                   R_p = 71492e3,
                                   no_CS_calculation = False
                                   ):

    m_e = 9.109e-31 # kg
    keV = 1.602e-16 # J == kg.m^2/s^2

    delta_longitude_observed = numpy.abs(lambda_observed_TEB - lambda_observed_MAW)

    if jovicentric_equator:
        (x_N,y_N,z_N) = spherical_to_cartesian(r_N, theta_N, phi_N)
        (x_S,y_S,z_S) = spherical_to_cartesian(r_S, theta_S, phi_S)
    else:
        if centrifugal_equator:
            theta_d = 9.3*1/3 # Centrifugal Equator tilt
        if magnetic_equator:
            theta_d = 9.3 
        phi_d = 204.2

        theta_N_transformed = 90-(theta_d*numpy.cos((phi_d-(360-phi_N))*numpy.pi/180.)+(90-theta_N))
        theta_S_transformed = 90-(theta_d*numpy.cos((phi_d-(360-phi_S))*numpy.pi/180.)+(90-theta_S))

        (x_N,y_N,z_N) = spherical_to_cartesian(r_N, theta_N_transformed, phi_N)
        (x_S,y_S,z_S) = spherical_to_cartesian(r_S, theta_S_transformed, phi_S)

    v_array = numpy.sqrt(2*E_array / m_e)

    # Outputs
    H_array = numpy.zeros((v_array.shape[0], z_lim_N_array.shape[0]))
    rho_0_array = numpy.zeros((v_array.shape[0], z_lim_N_array.shape[0]))
    n_0_ions_array = numpy.zeros((v_array.shape[0], z_lim_N_array.shape[0]))
    n_0_electrons_array = numpy.zeros((v_array.shape[0], z_lim_N_array.shape[0]))
    longitude_TEB_array = numpy.zeros((v_array.shape[0], z_lim_N_array.shape[0]))
    longitude_MAW_array = numpy.zeros((v_array.shape[0], z_lim_N_array.shape[0]))
    n_at_Moon_ions_array = numpy.zeros((v_array.shape[0], z_lim_N_array.shape[0]))
    n_at_Moon_electrons_array = numpy.zeros((v_array.shape[0], z_lim_N_array.shape[0]))
    
    
    for i_v, v in enumerate(v_array):
        for i_zlim_N, z_lim_N in enumerate(z_lim_N_array): 
                
            # determination of the MFL section in the CS and outside the CS
            equatorial_radius_N = numpy.sqrt(x_N**2 + y_N**2)
            equatorial_radius_S = numpy.sqrt(x_S**2 + y_S**2)
            mask_CS_N = (z_N <= z_lim_N) & (equatorial_radius_N > equatorial_radius_lim) # Second condition is to be sure no points close to the planet are taken
            mask_CS_S = equatorial_radius_S >= equatorial_radius_N[mask_CS_N].min()
    
            x_vals_CS = numpy.concatenate((x_N[mask_CS_N], x_S[mask_CS_S]))
            y_vals_CS = numpy.concatenate((y_N[mask_CS_N], y_S[mask_CS_S]))
            z_vals_CS = numpy.concatenate((z_N[mask_CS_N], z_S[mask_CS_S]))
            equatorial_vals = numpy.concatenate((equatorial_radius_N[mask_CS_N], equatorial_radius_S[mask_CS_S]))
            H =  (z_vals_CS.max() - z_vals_CS.min())/2 # Scale height in m
    
            # longitude of TEB and MAW if no CS
            if no_CS_calculation == True:
                (t_MAW_no_CS, t_TEB_no_CS,
                    Delta_t_no_CS,
                    longitude_MAW_no_CS, longitude_TEB_no_CS,
                    delta_longitude_no_CS) = delta_longitude_calculation_no_CS(r_N, theta_N, phi_N,
                                                                               r_S, theta_S, phi_S,
                                                                               rotation_rate,
                                                                               v)
                if verbose:
                    print(f'Associated values without Current Sheet: \n'
                            +r'λ$_{theorique}$ = ' + f'{360-phi_N[0]}°'+ '\n'
                            +r't$_{MAW}$ = ' + f'{t_MAW_no_CS} seconds'+'\n' 
                            +r't$_{TEB}$ = ' + f'{t_TEB_no_CS} seconds'+ '\n'
                            +r'λ$_{MAW}$ = ' + f'{longitude_MAW_no_CS}°'+ '\n'
                            +r'λ$_{TEB}$ = ' + f'{longitude_TEB_no_CS}°'+ '\n'
                            +r'$\Delta t$ = ' + f'{Delta_t_no_CS} seconds' + '\n'
                            +r'$\Delta \lambda$ = ' + f'{delta_longitude_no_CS} degrees')
    
                    
            # determination of rho_0 at the center of the CS based on minimzation of |Δλ_calculated - Δλ_observed|
            results = minimize(
                        objective_function_rho_0_from_delta_longitude,
                        rho_0_volumetric_first_guess, 
                        args=(r_0, H,
                                x_N, y_N, z_N, 
                                x_S, y_S, z_S,
                                B_total_N, B_total_S,
                                mask_CS_N, mask_CS_S,
                                v, rotation_rate, delta_longitude_observed,
                                disk, torus),
                        method='Nelder-Mead',                         # Optimization method
                        tol=1e-6                                      # Tolerance for convergence
                        )
        
            rho_0_volumetric_optimized = results.x[0]
            if verbose:
                print(f'({v:.2E} m/s (E = {v**2/2*m_e/keV:.2f} keV)')
                print(f'H: {H/R_p:.02f}')
                print(f'ρ_0_optimized: {rho_0_volumetric_optimized/1e6:.2E} kg/cm^3')
        
            # Determination of λ_MAW and λ_TEB based from rho_0_volumetric_optimized 
        
    
            longitude_MAW, longitude_TEB = longitude_MAW_TEB_calculation(x_N, y_N, z_N, B_total_N,
                                                                        x_S, y_S, z_S, B_total_S,
                                                                        mask_CS_N, mask_CS_S,
                                                                        r_0, H, rho_0_volumetric_optimized,
                                                                        v,
                                                                        rotation_rate,
                                                                        360-phi_N[0],
                                                                        disk = disk, torus = torus)
            delta_longitude_CS_MAW_TEB = numpy.abs(longitude_TEB-longitude_MAW)

            m_amu_ions, m_kg_ions, q_mean_ions = (single_ion_population_mass(r_0/R_p))
        
            ions_density = (rho_0_volumetric_optimized/1e6)/(m_kg_ions)  # as m_i >> m_e we consider that rho_0_volumetric_optimized is ~∑(ions_density*m_ions)
            electrons_density = ions_density*q_mean_ions # needs to have quasi-neutrality between - and + charges
            
            rho_at_Moon = rho_z(x_N[0], y_N[0], z_N[0], r_0, H, rho_0_volumetric_optimized, disk=disk, torus = torus)
            ions_density_at_Moon = (rho_at_Moon/1e6)/(m_kg_ions)
            electrons_density_at_Moon = ions_density_at_Moon*q_mean_ions
            if verbose:  
                print(f'ρ0: {rho_0_volumetric_optimized/1e6:.2E} kg.cm^-3') 
                print(f'n_0_ions: {ions_density:.2E} cm^-3')
                print(f'n_0_electrons: {electrons_density:.2E} cm^-3')
                print(f'λ (TEB) from UV measurement: {lambda_observed_TEB:.02f}\n'
                        +f'λ (MAW) from UV measurement: {lambda_observed_MAW:.02f}\n'
                        +f'Δλ observed: {numpy.abs(lambda_observed_TEB-lambda_observed_MAW):.02f}°')
                print(f'λ TEB: {longitude_TEB:.02f}°, \n'
                        +f'λ MAW: {longitude_MAW:.02f}°\n'
                        +f'Δλ: {delta_longitude_CS_MAW_TEB:.02f}°')
                
                print(f'n_ions @ Moon: {ions_density_at_Moon:.2E} cm^-3')
                print(f'n_electrons @ Moon: {electrons_density_at_Moon:.2E} cm^-3')
                print('\n \n \n')

            H_array[i_v, i_zlim_N] = H/R_p
            rho_0_array[i_v, i_zlim_N] = rho_0_volumetric_optimized/1e6
            n_0_ions_array[i_v, i_zlim_N] = electrons_density
            n_0_electrons_array[i_v, i_zlim_N] = electrons_density
            longitude_TEB_array[i_v, i_zlim_N] = longitude_TEB
            longitude_MAW_array[i_v, i_zlim_N] = longitude_MAW
            n_at_Moon_ions_array[i_v, i_zlim_N] = ions_density_at_Moon
            n_at_Moon_electrons_array[i_v, i_zlim_N] = electrons_density_at_Moon

        units_dict = {
            "H": f"Rp: {R_p/1e3} km", # Height Scale in Rp
            "rho_0": "kg.cm^-3", # total mass density in the center of the PS
            "n_0_ions": "cm^-3", # ions density in the PS
            "n_0_electrons": "cm^-3", # electrons density in the PS
            "longitude_TEB": "degrees", # TEB longitude
            "longitude_MAW": "degrees", # MAW longitude
            "n_at_Moon_ions": "cm^-3", # ions density at the moon orbit
            "n_at_Moon_electrons": "cm^-3" # electrons density at the moon orbit
        }
                
    
    return (H_array, # Height Scale in Rp
            rho_0_array, #total mass density in the center of the PS in kg.cm^-3
            n_0_ions_array, # ions density in the PS cm^-3
            n_0_electrons_array, # electrons density in the PS cm^-3
            longitude_TEB_array, # TEB longitude
            longitude_MAW_array, # MAW longitude
            n_at_Moon_ions_array, # ions density at the moon orbit cm^-3
            n_at_Moon_electrons_array, # electrons density at the moon orbit cm^-3
            units_dict)