#Previous version:
#
##############################################################################
#New update:
# The Sigma Matrix is for the horizontal model

import numpy as np

from numpy.linalg import inv

#Function to solve the Dyson Equation based on Green Functions

def green_matrix(energy, gamma_parameter, H, polarization_left, polarization_right):

    #Extract the dimension of the Hamiltonian matrix:
    n = len(H)

    # Matrix M which is the Hamiltonian + omega + imaginaty*eta

    M = np.zeros((n,n), dtype=complex)

    for i in range(n):
        for j in range(n):
            M[i, j] = H[i, j]

    # print("The matrix Hcen is:")
    # print(M)

    #Now write the matrix E:

    E = np.zeros((n,n), dtype=complex)

    for i in range(0, n):
        E[i,i] = energy 
        #print(E)

    #print(E)

    #Definition of the Gamma Matrices:
    Gamma_left = np.zeros((n,n),dtype=complex)
    Gamma_left[0][0] = (gamma_parameter/2)*((1+polarization_left))
    Gamma_left[1][1] = (gamma_parameter/2)*((1-polarization_left))
    #print(Gamma_left)

    Gamma_right = np.zeros((n,n),dtype=complex)
    Gamma_right[n-2][n-2]= (gamma_parameter/2)*((1+polarization_right))
    Gamma_right[n-1][n-1] = (gamma_parameter/2)*((1-polarization_right))
    #print(Gamma_right)

    #Sigma i.e. Self energy matrix definition
    Sigma = np.zeros((n,n),dtype=complex)

    #Sigma matrix:
    Sigma = -(1j)*(1/2)*(Gamma_left+Gamma_right)

    #Get the inverse i.e. the Green function Matrix

    Green_r = inv(E- M - Sigma)
    #print(E-M-Sigma)

    #print(Green_r)

    #Now get the transmission:
    #Define Green Advanced Function, Level Width Function (Gamma) for left and right
    Green_a=Green_r.conj().T

    #print(Green_a)

    A = np.matmul(Green_r,Gamma_left)
    B = np.matmul(Green_a,Gamma_right)

    C = np.matmul(A,B)

    # print("The matrix G^r Gamma_L G^a Gamma_R:")
    # print(C)

    #transm = np.sum(np.diag(C))
    transm = np.abs(np.diag(C))

    #print(transm)

    return transm[n-2], transm[n-1]



#Function for the Landauer Transmission
#def Transmission(energy_array, gamma_parameter, eta_parameter, H, dim):

def Transmission(energy_array, gamma_parameter, H, polarization_left, polarization_right):

    #Transmission axis:
    transmission_axis_up = np.zeros(len(energy_array))
    transmission_axis_down = np.zeros(len(energy_array))

    #Run in the array:
    for i in range(len(energy_array)):
        #print(i)
        #Transmission as the trace of the part of the matrix of Green and Wide band function:
        transmission_axis_up[i] = green_matrix(energy_array[i], gamma_parameter, H, polarization_left, polarization_right)[0]
        transmission_axis_down[i] = green_matrix(energy_array[i], gamma_parameter, H, polarization_left, polarization_right)[1]

    #print(transmission_axis_up)

    return transmission_axis_up, transmission_axis_down

