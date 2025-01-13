#Import libraries
import numpy as np

import Ln_functions as Ln
from importlib import reload
reload(Ln)

#Set the values for the important constants:

#Set the Boltzmann constant:
kb = 0.08617 # meV * K^{-1}

#Set the Planck Constant:
hbar = 6.582*(10**(-13)) #meV.s

#Set the electron charge constant as the unity:
e = 1 

# Spin thermoelectric coefficient when \Delta V^s = 0

def S(mu_array, Temperature_parameter, energy_array, polarization_left, polarization_right, gamma_parameter, Hamiltonian):

    result = -(1/(e*Temperature_parameter))*( Ln.L_function_up(1, mu_array, Temperature_parameter, energy_array,polarization_left, polarization_right, gamma_parameter, Hamiltonian) +\
                                             Ln.L_function_down(1, mu_array, Temperature_parameter, energy_array,polarization_left, polarization_right, gamma_parameter, Hamiltonian))/\
                                            ( Ln.L_function_up(0, mu_array, Temperature_parameter, energy_array,polarization_left, polarization_right, gamma_parameter, Hamiltonian) +\
                                              Ln.L_function_down(0, mu_array, Temperature_parameter, energy_array,polarization_left, polarization_right, gamma_parameter, Hamiltonian) )
    
    # Reemplazar NaNs con ceros
    result_no_nans = np.nan_to_num(result)    
    return result_no_nans

# Spin thermoelectric coefficients when \Delta V^s \neq 0:

# Seebeck UP:

def S_up(mu_array, Temperature_parameter, energy_array, polarization_left, polarization_right, gamma_parameter, Hamiltonian):

    result = (1/(e*Temperature_parameter))*( Ln.L_function_up(1, mu_array, Temperature_parameter, energy_array,\
                                            polarization_left, polarization_right, gamma_parameter, Hamiltonian)/\
                                            Ln.L_function_up(0, mu_array, Temperature_parameter, energy_array,\
                                                           polarization_left, polarization_right, gamma_parameter, Hamiltonian) )
    
    # Reemplazar NaNs con ceros
    result_no_nans = np.nan_to_num(result)    
    return result_no_nans

# Seebeck DOWN:

def S_down(mu_array, Temperature_parameter, energy_array, polarization_left, polarization_right, gamma_parameter, Hamiltonian):

    result = (1/(e*Temperature_parameter))*( Ln.L_function_down(1, mu_array, Temperature_parameter, energy_array, \
                                                             polarization_left, polarization_right, gamma_parameter, Hamiltonian)/\
                                            Ln.L_function_down(0, mu_array, Temperature_parameter, energy_array,\
                                                             polarization_left, polarization_right, gamma_parameter, Hamiltonian) )
    
    # Reemplazar NaNs con ceros
    result_no_nans = np.nan_to_num(result)    
    return result_no_nans