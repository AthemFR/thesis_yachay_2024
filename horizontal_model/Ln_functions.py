#Import important libraries
import numpy as np
from importlib import reload

import functions_spin_p as fy
from importlib import reload

reload(fy)

#Set the values for the important constants:

#Set the Boltzmann constant:
kb = 0.08617 # meV * K^{-1}

#Set the Planck Constant:
hbar = 6.582*(10**(-13)) #meV.s

#Set the electron charge constant as the unity:
e = 1 

#The spin dependent Ln function is defined based on the integration of the following function:

# Fermi derivative respect to the Energy

def fermi_derivative(energy_array, mu_parameter, Temperature_parameter):
    
    #Convert the temperature in to Energy value by product with Boltzman constant
    temperature_kb = Temperature_parameter*kb

    result= np.zeros(len(energy_array))
    
    for i in range(len(energy_array)):
        result[i] = ( np.exp((energy_array[i]-mu_parameter)/(temperature_kb)) )/(\
             (1+np.exp((energy_array[i]-mu_parameter)/(temperature_kb)) )**2*temperature_kb)

    return result


# Ln function spin UP:

def L_function_up(order, mu_array, Temperature_parameter, energy_array,\
                  polarization_left, polarization_right, gamma_parameter, Hamiltonian):
    
    #Array of transmissions up:
    unique_transmission = fy.Transmission(energy_array, gamma_parameter, \
                                          Hamiltonian, polarization_left, polarization_right)[0]
    transmissions=[]

    for i in range(len(mu_array)):
        transmissions.append(unique_transmission)

    transmissions= np.array(transmissions)

    #Array of Fermi Derivative arrays
    fermi_derivatives = []

    for i in range(len(mu_array)):
        fermi_derivatives.append(fermi_derivative(energy_array,\
             mu_array[i], Temperature_parameter))

    fermi_derivatives= np.array(fermi_derivatives)
    #fermi_derivatives.shape

    #Array of (Energy-mu)^{n} arrays:

    ener_min_mu = []

    for i in range(len(mu_array)):
        ener_min_mu.append(np.power((energy_array-mu_array[i]),order))

    ener_min_mu= np.array(ener_min_mu)
    ener_min_mu.shape

    #Product of all arrays:

    product = fermi_derivatives*ener_min_mu*transmissions

    #Integration over the Energy array:
    dE = energy_array[1]-energy_array[0]

    result =(dE/(2*np.pi))*np.sum(product, axis=1)

    #We avoid the using of hbar constant and we set it just as constant!

    result.shape #meV^{n}

    return result


# Ln function spin DOWN:

def L_function_down(order, mu_array, Temperature_parameter, energy_array,\
                  polarization_left, polarization_right, gamma_parameter, Hamiltonian):
    
    #Array of transmissions up:
    unique_transmission = fy.Transmission(energy_array, gamma_parameter,\
                                          Hamiltonian, polarization_left, polarization_right)[1]
    transmissions=[]

    for i in range(len(mu_array)):
        transmissions.append(unique_transmission)

    transmissions= np.array(transmissions)

    #Array of Fermi Derivative arrays
    fermi_derivatives = []

    for i in range(len(mu_array)):
        fermi_derivatives.append(fermi_derivative(energy_array,\
             mu_array[i], Temperature_parameter))

    fermi_derivatives= np.array(fermi_derivatives)
    fermi_derivatives.shape

    #Array of (Energy-mu)^{n} arrays:

    ener_min_mu = []

    for i in range(len(mu_array)):
        ener_min_mu.append(np.power((energy_array-mu_array[i]),order))

    ener_min_mu= np.array(ener_min_mu)
    ener_min_mu.shape

    #Product of all arrays:

    product = fermi_derivatives*ener_min_mu*transmissions

    #Integration over the Energy array:
    dE = energy_array[1]-energy_array[0]

    result =(dE/(2*np.pi))*np.sum(product, axis=1)

    #We avoid the using of hbar constant and we set it just as constant!

    result.shape #meV^{n}

    return result