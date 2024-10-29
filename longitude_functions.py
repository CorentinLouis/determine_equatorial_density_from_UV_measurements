import numpy
from scipy.optimize import minimize


from density import rho_z
from alfven_velocity import v_alfven_calculation
from magnetic_field_functions import length_magnetic_field_line, delta_length_magnetic_field_line
from coordinates_system_transformation import spherical_to_cartesian, cartesian_to_spherical

def delta_longitude_calculation_no_CS(r_N, theta_N, phi_N,
                                      r_S, theta_S, phi_S,
                                        rotation_rate_callisto,
                                        v_electrons
                                        ):
    c = 3e8 # m/s
    longitude_theorique_callisto = 360-phi_N[-1]

    length_mfl_B_N = length_magnetic_field_line(r_N, theta_N, phi_N, spherical = True, cartesian = False)
    length_mfl_B_S = length_magnetic_field_line(r_S, theta_S, phi_S, spherical = True, cartesian = False)


    t_TEB_no_CS = length_mfl_B_S/c + (length_mfl_B_S+length_mfl_B_N)/v_electrons
    t_MAW_no_CS = length_mfl_B_N/c
    Delta_t_no_CS = t_TEB_no_CS - t_MAW_no_CS

    delta_longitude_no_CS = Delta_t_no_CS * rotation_rate_callisto

    longitude_TEB_no_CS = (360-phi_N[0])-t_TEB_no_CS*rotation_rate_callisto
    longitude_MAW_no_CS = (360-phi_N[0])-t_MAW_no_CS*rotation_rate_callisto
    
    return t_MAW_no_CS, t_TEB_no_CS, Delta_t_no_CS, longitude_MAW_no_CS, longitude_TEB_no_CS, delta_longitude_no_CS

# Function to calculate delta longitude
def delta_longitude_calculation(x_N, y_N, z_N, 
                                  x_S, y_S, z_S,
                                  mask_CS_N, mask_CS_S,
                                  v,
                                  v_A_N,
                                  v_A_S,
                                  rotation_rate_callisto,
                                  delta_L = 1):
    """
    All input values needs to be in SI units
    delta_L (default is 1: no modification): modification factor to take into account that the Length of the magnetic field line in the plasma sheet is in fact longer, due to the tilt (theta = atan(Ma)) of the Alfvén wing 
    """
    c=3e8 # m/s 
    
    delta_length_mfl_B_N_CS = delta_length_magnetic_field_line(x_N[mask_CS_N], y_N[mask_CS_N], z_N[mask_CS_N], cartesian = True, spherical = False)
    #delta_length_mfl_B_N_out_CS = delta_length_magnetic_field_line(x_N[~mask_CS_N], y_N[~mask_CS_N], z_N[~mask_CS_N], cartesian = True, spherical = False)*R_J
    delta_length_mfl_B_S_CS = delta_length_magnetic_field_line(x_S[mask_CS_S], y_S[mask_CS_S], z_S[mask_CS_S], cartesian = True, spherical = False)
    #delta_length_mfl_B_S_out_CS = delta_length_magnetic_field_line(x_S[~mask_CS_S], y_S[~mask_CS_S], z_S[~mask_CS_S], cartesian = True, spherical = False)*R_J


    length_mfl_B_N_CS = length_magnetic_field_line(x_N[mask_CS_N], y_N[mask_CS_N], z_N[mask_CS_N], cartesian = True, spherical = False)
    length_mfl_B_N_out_CS = length_magnetic_field_line(x_N[~mask_CS_N], y_N[~mask_CS_N], z_N[~mask_CS_N], cartesian = True, spherical = False)
    length_mfl_B_S_CS = length_magnetic_field_line(x_S[mask_CS_S], y_S[mask_CS_S], z_S[mask_CS_S], cartesian = True, spherical = False)
    length_mfl_B_S_out_CS = length_magnetic_field_line(x_S[~mask_CS_S], y_S[~mask_CS_S], z_S[~mask_CS_S], cartesian = True, spherical = False)


    length_mfl_B_S = length_mfl_B_N_CS + length_mfl_B_N_out_CS
    length_mfl_B_N = length_mfl_B_S_CS + length_mfl_B_S_out_CS


    # Time calculation (t_TEB_CS and t_MAW_CS)
    t_TEB_CS = numpy.sum(delta_length_mfl_B_S_CS / v_A_S[:-1])*delta_L + length_mfl_B_S_out_CS / c  + (length_mfl_B_S + length_mfl_B_N)/ v


    t_MAW_CS = numpy.sum(delta_length_mfl_B_N_CS / v_A_N[:-1])*delta_L + length_mfl_B_N_out_CS / c


    # Delta time and delta longitude
    Delta_t_CS = numpy.abs(t_TEB_CS - t_MAW_CS)
    delta_longitude_MAW = t_MAW_CS * rotation_rate_callisto
    delta_longitude_TEB = t_TEB_CS * rotation_rate_callisto
    delta_longitude_CS_MAW_TEB = Delta_t_CS*rotation_rate_callisto
    
    
    return t_MAW_CS, t_TEB_CS, Delta_t_CS, delta_longitude_MAW, delta_longitude_TEB, delta_longitude_CS_MAW_TEB



c = 3e8 #m/s
# Define the nested function for optimizing rho_0 given a fixed z_lim_N
def optimize_rho_0(z_lim_N, r_0, H, rho_0_volumetric_first_guess,
                   x_N, y_N, z_N, B_N,
                   x_S, y_S, z_S, B_S,
                   mask_CS_N, mask_CS_S,
                   v, rotation_rate_callisto, 
                   delta_longitude_observed, disk=False, torus=False,
                   delta_L = 1):
    # Define the inner objective function for rho_0 minimization
    def objective_function_rho_0_from_delta_longitude_inner(rho_0, r_0, H,
                                       x_N, y_N, z_N, 
                                       x_S, y_S, z_S,
                                       B_N, B_S,
                                       mask_CS_N, mask_CS_S,
                                       v, rotation_rate_callisto, delta_longitude_observed,
                                       disk = False, torus = False,
                                       delta_L = 1):

        #if numpy.any(mask_CS_N):
        rho_CS_N = rho_z(x_N[mask_CS_N], y_N[mask_CS_N], z_N[mask_CS_N], r_0, H, rho_0, disk=disk, torus = torus)
        rho_CS_N_volumetric = rho_CS_N
        v_A_N = v_alfven_calculation(B_N[mask_CS_N], rho_CS_N_volumetric) 
        #else:
        #    v_A_N = c
        #if numpy.any(mask_CS_S):
        rho_CS_S = rho_z(x_S[mask_CS_S], y_S[mask_CS_S], z_S[mask_CS_S], r_0, H, rho_0, disk=disk, torus = torus)
        rho_CS_S_volumetric = rho_CS_S
        v_A_S = v_alfven_calculation(B_S[mask_CS_S], rho_CS_S_volumetric)
        #else:
        #    v_A_S = c

        
        
        (_, _, _,_ ,_ , delta_longitude_CS_MAW_TEB) = delta_longitude_calculation(x_N, y_N, z_N, 
                                                                                    x_S, y_S, z_S,
                                                                                    mask_CS_N, mask_CS_S,
                                                                                    v,
                                                                                    v_A_N,
                                                                                    v_A_S,
                                                                                    rotation_rate_callisto,
                                                                                    delta_L = delta_L)



        
        # Minimize the difference between calculated and observed delta longitude
        #print(f'|Δλ_calculated - Δλ_observed|: {(delta_longitude_CS_MAW_TEB - delta_longitude_observed) ** 2:.02f}')
        return (delta_longitude_CS_MAW_TEB - delta_longitude_observed) ** 2

    # Perform the rho_0 optimization
    result = minimize(
        objective_function_rho_0_from_delta_longitude_inner,
        rho_0_volumetric_first_guess,
        args=(r_0, H, 
                x_N, y_N, z_N,
                x_S, y_S, z_S,
                B_N, B_S,
                mask_CS_N, mask_CS_S,
                v, rotation_rate_callisto,
                delta_longitude_observed, disk, torus, delta_L),
        method='Nelder-Mead',
        tol=1e-6
    )

    # Extract the optimized rho_0 value
    return result.x[0]


# Objective function to minimize
def objective_function_rho_0_from_delta_longitude(rho_0, r_0, H,
                                       x_N, y_N, z_N, 
                                       x_S, y_S, z_S,
                                       B_N, B_S,
                                       mask_CS_N, mask_CS_S,
                                       v, rotation_rate_callisto, delta_longitude_observed,
                                       disk = False, torus = False,
                                       delta_L = 1):



    rho_CS_N = rho_z(x_N[mask_CS_N], y_N[mask_CS_N], z_N[mask_CS_N], r_0, H, rho_0, disk=disk, torus = torus)
    rho_CS_S = rho_z(x_S[mask_CS_S], y_S[mask_CS_S], z_S[mask_CS_S], r_0, H, rho_0, disk=disk, torus = torus)

    rho_CS_N_volumetric = rho_CS_N
    rho_CS_S_volumetric = rho_CS_S

    v_A_N = v_alfven_calculation(B_N[mask_CS_N], rho_CS_N_volumetric) 
    v_A_S = v_alfven_calculation(B_S[mask_CS_S], rho_CS_S_volumetric)

    
    
    (_, _, _,_ ,_ , delta_longitude_CS_MAW_TEB) = delta_longitude_calculation(x_N, y_N, z_N, 
                                                                                x_S, y_S, z_S,
                                                                                mask_CS_N, mask_CS_S,
                                                                                v,
                                                                                v_A_N,
                                                                                v_A_S,
                                                                                rotation_rate_callisto,
                                                                                delta_L=delta_L)



    
    # Minimize the difference between calculated and observed delta longitude
    #print(f'|Δλ_calculated - Δλ_observed|: {(delta_longitude_CS_MAW_TEB - delta_longitude_observed) ** 2:.02f}')
    return (delta_longitude_CS_MAW_TEB - delta_longitude_observed) ** 2


# Define the main objective function for z_lim_N optimization
def objective_function_zlimN_from_longitude(z_lim_N, r_0, rho_0_volumetric_first_guess,
                                           r_N, theta_N, phi_N, B_N,
                                           r_S, theta_S, phi_S, B_S,
                                           E, rotation_rate_callisto,
                                           longitude_instantaneous, longitude_observed_MAW, longitude_observed_TEB,
                                           equatorial_radius_lim = 10*71492e3,
                                           R_p = 71492e3,
                                           MAW=False, TEB=False, disk=False, torus=False,
                                           jovicentric_equator=False, centrifugal_equator=False, magnetic_equator=False,
                                           delta_L = 1):

    m_e = 9.109e-31 # kg
    delta_longitude_observed = numpy.abs(longitude_observed_MAW-longitude_observed_TEB)
    v = numpy.sqrt(2*E / m_e)
    if jovicentric_equator:
        (x_N,y_N,z_N) = spherical_to_cartesian(r_N, theta_N, phi_N)
        (x_S,y_S,z_S) = spherical_to_cartesian(r_S, theta_S, phi_S)
    else:
        if centrifugal_equator:
            theta_d = 9.3*1/3 # Centrifugal Equator tilt
        elif magnetic_equator:
            theta_d = 9.3 
        
        else:
            raise ValueError("Equator type not specified. Choose between jovicentric, centrifugal, or magnetic.")
        phi_d = 204.2

        theta_N_transformed = 90-(theta_d*numpy.cos((phi_d-(360-phi_N))*numpy.pi/180.)+(90-theta_N))
        theta_S_transformed = 90-(theta_d*numpy.cos((phi_d-(360-phi_S))*numpy.pi/180.)+(90-theta_S))

        (x_N,y_N,z_N) = spherical_to_cartesian(r_N, theta_N_transformed, phi_N)
        (x_S,y_S,z_S) = spherical_to_cartesian(r_S, theta_S_transformed, phi_S)

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

    # Optimize rho_0 for the current z_lim_N
    
    
    rho_0_volumetric_optimized = optimize_rho_0(z_lim_N, r_0, H, rho_0_volumetric_first_guess,
                                                x_N, y_N, z_N, B_N,
                                                x_S, y_S, z_S, B_S,
                                                mask_CS_N, mask_CS_S,
                                                v, rotation_rate_callisto,
                                                delta_longitude_observed, disk=disk, torus=torus, delta_L = delta_L)

    # Calculate the longitudes using the optimized rho_0
    longitude_MAW, longitude_TEB = longitude_MAW_TEB_calculation(
        x_N, y_N, z_N, B_N, x_S, y_S, z_S, B_S, mask_CS_N, mask_CS_S, 
        r_0, H, rho_0_volumetric_optimized, v, rotation_rate_callisto, 
        longitude_instantaneous, disk=disk, torus=torus, delta_L= delta_L
    )

    # Determine which longitude to compare
    if MAW:
        longitude = longitude_MAW
        longitude_observed = longitude_observed_MAW
    if TEB:
        longitude = longitude_TEB
        longitude_observed = longitude_observed_TEB

    # Compute the squared difference as the objective function value
    return (longitude - longitude_observed) ** 2

def longitude_MAW_TEB_calculation(x_N, y_N, z_N, B_total_N,
                                      x_S, y_S, z_S, B_total_S,
                                      mask_CS_N, mask_CS_S,
                                      r_0, H, rho_0_volumetric,
                                      v_electrons,
                                      rotation_rate_callisto,
                                      longitude_callisto,
                                      disk = False,
                                      torus = False,
                                      delta_L = 1):

    rho_CS_N_volumetric = rho_z(x_N[mask_CS_N], y_N[mask_CS_N], z_N[mask_CS_N], r_0, H, rho_0_volumetric, disk=disk, torus = torus)
    rho_CS_S_volumetric = rho_z(x_S[mask_CS_S], y_S[mask_CS_S], z_S[mask_CS_S], r_0, H, rho_0_volumetric, disk=disk, torus = torus)

    v_A_N = v_alfven_calculation(B_total_N[mask_CS_N], rho_CS_N_volumetric) 
    v_A_S = v_alfven_calculation(B_total_S[mask_CS_S], rho_CS_S_volumetric)


    (_, _, _, 
            delta_longitude_MAW, delta_longitude_TEB, _
            ) = delta_longitude_calculation(x_N, y_N, z_N, 
                                    x_S, y_S, z_S,
                                    mask_CS_N, mask_CS_S,
                                    v_electrons,
                                    v_A_N,
                                    v_A_S,
                                    rotation_rate_callisto, 
                                    delta_L = delta_L)


    # Delta time and delta longitude
    
    longitude_MAW = longitude_callisto - delta_longitude_MAW
    longitude_TEB = longitude_callisto - delta_longitude_TEB
    return longitude_MAW, longitude_TEB
