#Previous version:
#
##############################################################################
#New update:
# The Sigma Matrix is for the horizontal model

import numpy as np

from numpy.linalg import inv

#Function to solve the Dyson Equation based on Green Functions

def green_matrix(energy, eta_parameter, H):

    #Extract the dimension of the Hamiltonian matrix:
    n = len(H)
    #print("The size is:", n)

    # Matrix M which is the Hamiltonian
    M = np.zeros((n,n), dtype=complex)

    for i in range(n):
        for j in range(n):
            M[i, j] = H[i, j]
    #print("The matrix Hcen is:")
    #print(M)

    #Now write the matrix E which is E + imaginary*eta:

    E = np.zeros((n,n), dtype=complex)
    #print(E)

    for i in range(0, n):
        #print(i)
        E[i,i] = energy + (1j)*eta_parameter
        #print(E)

    #print(E)

    #Get the inverse i.e. the Green function Matrix

    Green_r = inv(E- M)
    #print(E-M-Sigma)

    #print(Green_r)

    return Green_r

#Function for the DOS

def DOS(energy_array, eta_parameter, H):
    #DOS axis:
    dos_axis = np.zeros(len(energy_array))
    #Run in the array:
    for i in range(len(energy_array)):
        #print(i)
        Green_retarded = green_matrix(energy_array[i], eta_parameter, H)
        #DOS as the axis, remember that we are making DOS = -Im(Tr(G)) for each energy in the energy array
        dos_axis[i] = -np.sum(np.diag(Green_retarded).imag)
    return dos_axis



