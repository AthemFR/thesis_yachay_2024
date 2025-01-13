# Update:
# Solving some problems related to the Hamiltonian renormalized
# En la anterior version la diagonal solo era los elementos de 
# las energias renormalizadas y nada mas.
# esto en el caso que no sea el modelo de Geyer funciona ya que en la
# diagonal no habria nada, pero en el modelo de Geyer si que afecta


#Import libraries:
import numpy as np

# Define the function for renormalize each site energy:

def ren(array, energy):
    result = array[0] + ( (array[2]**2) /(energy - array[1]) )
    return result
    
#Implement the previous function in the Hamiltonian for renormalization:

def hamiltonian_ren(Hamiltonian, sites_ren, energy):
    n = len(Hamiltonian)
    solution_1 = np.zeros((n,n), dtype=complex)

    for i in range(n):
        solution_1[i][i] = ren(sites_ren[i],energy)

    result = solution_1 + Hamiltonian

    return result