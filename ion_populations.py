import pandas as pd 
import numpy as np 

def ion_equatorial_density(specie,rho):
    # Density at the centrifugal equator
    # From Dougherty et al.(2017), for R>= 6 
    
    # Equation & coefficient from Bagenal & Delamere (2011)
    a1,a2,a3 = 1987,14,0.05
    b1,b2,b3 = -8.2,-3.2,-0.65

    # Ion & electron equatorial centrifugal density constant from Dougherty et al.(2017) 
    coeff_density = pd.DataFrame(columns=['e-','O+','O++','S+','S++','S+++','H+','Na+','O+_hot'], index=range(4))
    coeff_density.iloc[:,:] = np.transpose(
            np.array([[1,2451,-6.27,4.21],   [0.24,592,-7.36,0.368],[0.03,76.3,-6.73,0.086],
                      [0.07,163,-6.81,0.169],[0.22,538,-6.74,0.598],[0.04,90.7,-6.21,0.165],
                      [0.02,50.6,-5.31,0.212],[0.04,97.2,-6.75,0.106],[0.06,134,-4.63,1.057]]))
    
    
    ni,a,b,c = tuple(coeff_density[specie])
    
    n_2011 = a1*(rho/6)**b1 + a2*(rho/6)**b2 + a3*(rho/6)**b3
    if rho >= 6 and rho<= 15.2  : 
        n = a*(rho/6)**b 
        
    if rho > 15.2 : 
        n = c*n_2011

    
    return n


def single_ion_population_mass(rho):
    
    species_model = ['S+','S++','S+++','O+','O++','H+','Na+']
    m_ion = {'O+':16,'O++':16,'S+':32,'S++':32,'S+++':32,'H+':1,'Na+':23}
    
    z_ion = {'O+':1,'O++':2,'S+':1,'S++':2,'S+++':3,'H+':1,'Na+':1}

    # Density [cm-3] of each species
    n_i = [ion_equatorial_density(sp, rho) for sp in species_model] # cm-3
    n_tot = sum(n_i) # cm-3
    n_relative = np.array(n_i)/n_tot # % 
    
    # Mean mass [amu]
    m_amu = 0 
    # Mean charge
    q_mean = 0
    for i in range(len(species_model)):
        m_amu+=n_relative[i]*m_ion[species_model[i]]
        q_mean+=n_relative[i]*z_ion[species_model[i]]
    
    

    # Convert amu to kg 
    m_kg = m_amu*1.66054e-27  
    
    return m_amu, m_kg, q_mean
