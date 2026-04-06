# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 23:35:31 2026

@author: DELL
"""

import numpy as np
from functools import reduce
from scipy.linalg import expm
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from collections import Counter
import pandas as pd
from pathlib import Path

#%% code up winding number due to bulk hamiltonian 

def J_tilte_a(t, t_a, pulse_length, l, T, pulse_amp): # t_a define the START POINT in each period T for pulse a
    t_p = np.mod (t, T) #for periodic changings, but there's problem for Jz needed to be figured out

    t_start_a = l * T + t_a
    
    t_end_a = t_start_a + pulse_length
    
    if t_end_a < T:
    
#    if t_start_a < t < t_end:        
#        return  pulse_amp   
#    else:
#        return 0 -- this if else is not suitable when t is an array

        return np.where((t_p >= t_start_a) & (t_p < t_end_a), pulse_amp, 0)
    #np.where(condition, x, y) if condition is true : take x(pulse amp here) ; 
    #                          if condition is false: take y (0)3.+1
    if t_end_a > T:
        
        t_end_a_new = t_end_a - T
        
        return np.where((t_p >= t_start_a) & (t_p < T) | (t_p >= 0) & (t_p < t_end_a_new), pulse_amp, 0)
    
#parameter adjustion:
#(here always take T = 1, l = 0, set the numerical platform of time, also unify the y axis as epsilon in fig_1)
T = 1
l = 0
#(as said in the paragraph above (16), t_x = 0, t_y = T/3, t_z = 2*T/3)
t_x = 0
t_y = T/3
t_z = T * (2/3)
    
delta_t = T/2

J_a = 2* ((0.9 * np.pi ))/ ((T * 0.5) * 4)#  2* ((0.9 * np.pi ))/ ((T * 0.5) * 4)-- this is 2x the value in paper 

                                       # Q: Can I vary this J_a to perform perturb on Hamiltonian?
                                       # Ans: I don't think so. I think this might work: 1.change the pulse amplitude by 1%
                                       #                                                 2.shift the timing by 0.01T
                                       #                                                 3.add a tiny disorder term


#%% (here change H_test to real H, which is parameterised by k:
#(firstly try to use the [H] build a list of 6 H; 
#(another method is that to firstly write H1/U1 to H6/U6, then U_full = U1U2U3...U6)


# QUESTION: should I define a function: array of Hamiltonians? known H_i depend on Jx, Jy, Jz, k
# -- firstly we try only do as quantities
# here t_x,y,z are setted, delta_t are setted, J_x,y,z are setted 
def J_tilte_x(t):
    return J_tilte_a(t, t_x, delta_t, l, T, J_a)
def J_tilte_y(t):
    return J_tilte_a(t, t_y, delta_t, l, T, J_a)
def J_tilte_z(t):
    return J_tilte_a(t, t_z, delta_t, l, T, J_a) 
# didn't write a loop here because J_a and J_start are also variables, can be different/change 
# we defined fixed l = 1, T_length = T = 1 
#%%
# try plot: t is from 0 to 1, J_tilte_x, y, z 

t_test = np.linspace(-2, 4, 100) 
# function of J_a need to be adjusted as in that t can be out of range of one period(0,T)

plt.plot(t_test, J_tilte_x(t_test), label='J_tilte_x')
plt.plot(t_test, J_tilte_y(t_test), label='J_tilte_y')
plt.plot(t_test, J_tilte_z(t_test), label='J_tilte_z')

plt.xlabel('t')
plt.ylabel('J_tilte_a')
plt.title('Plot of J_tilta_a(pulse of hopping amp) vs t')
plt.legend()
plt.grid(True)
plt.show()

#%% some test try figure properties of J func
t_specific_time_point_test = np.pi - 3

call_func_J_tilte_x = J_tilte_z(t_specific_time_point_test)

print (call_func_J_tilte_x, ' see if J_tilte_a = x return a scalar multiply the Pauli matrices')
#plot succeed

#%%

# **a_0 is lattice spacing**, which is the quantity ensures translational invariance(in x direction):
# H(r+ a_0 *x_hat, r' + a_0 *x_hat, t) = H(r, r', t) -- here H is Hamiltonian and x_hat is unit vector in x direction

# Q:1.  what is r, r'(both vectors) specifically here? do they refer to interactions between r and r'? 
# Q:2.  does (r-r' )'s modulus<l shows some kind of interaction range between two sites/positions? 
#(just like the interaction in SSH model, but in SSH model H11 doesn't represent real position(as in meter or...) - )

#(- it represent onsite interaction: hopping amplitude between site 1 and site 1 )
# a_0 is also called 'lattice_length' in the subsequent codes, (coz a_0 was used directly in parameter in functions?

#lattice_length = 1
a_0 = 1
stage_length = T/6 #(numerical calculations)
# due to equation (28)
#stage_number = 6 #(just N afterwards)
#Pauli matrices:
pauli_x = np.array([[0.0, 1.0], [1.0, 0.0]])
pauli_y = np.array([[0.0, -1j], [1j, 0.0]])
pauli_z = np.array([[1.0, 0.0], [0.0, -1.0]])

# Pulse setting Q: where is J_0? Please find it later
    
#%%
#-----K BASIS HAMILTONIAN MATRIX H 2X2-----
#Here is the 2by2 Hamiltonian that describes the [sublattices coupling] in momentum space
#            -- we will refer this as '2by2 Hamiltonian' /'2by2 H' in the future

#delta_t is the actual T/2, pulse length of our case 
def H_matrix_2by2(k_x, k_y, t, a_0): 
    
    def J_tilte_x(t):
        return J_tilte_a(t, t_x, delta_t, l, T, J_a)
    def J_tilte_y(t):
        return J_tilte_a(t, t_y, delta_t, l, T, J_a)
    def J_tilte_z(t):
        return J_tilte_a(t, t_z, delta_t, l, T, J_a) 
    
    H = (J_tilte_x(t) + J_tilte_y(t) * np.cos (k_x * a_0) + J_tilte_z(t) * np.cos(k_y * a_0)) * pauli_y + (J_tilte_y(t) * np.sin (k_x * a_0) + J_tilte_z(t) * np.sin(k_y * a_0)) * pauli_x
    #(28)
    
    return H

# generate a list of t which takes 6 samples of static H(1 to 6)

#t = T/6 * i + T/12 (T/12 is for if there is any confusion in the boundary points of the six stages)


# run a small test for the future 2N 2N matrix
H_test_for_entrices = H_matrix_2by2(np.pi/2, np.pi/2, 5/12, 1)
a_11 = H_test_for_entrices[0][0]
a_12 = H_test_for_entrices[0][1]
#(1.413716694115407-1.413716694115407j) (1.413716694115407-1.413716694115407j)

#%% Guess: bulk-defect as in bulk hamiltonian's winding number decide what will happen when we introduce defect
#hence we don't need to calculate winding number through the one with time vortex; we can try: without time vortex
# write H_eff -- H_eff is calculated from U(k, T)

def U_full_discretisation(k_x, k_y, num_time_stages): # which means along lhs to rhs, from (1) earlier/smaller t to (6) later/larger t

    U_list_N_matrices = []
    t_desicretisation= T /num_time_stages
#    time_stage = num_time_stages
    
    for i in range (num_time_stages):
        
    
        t_i = T - t_desicretisation * (i + 1)
    
        H_i = H_matrix_2by2(k_x, k_y, t_i, a_0)
#expm(-1j* H_i_test* delta_t_test)
        
        U_i = expm(-1j* H_i* (t_desicretisation)) #!!!!! # here using exponential form construct each of U_1 to U_6
        
        U_list_N_matrices.append(U_i)
        
#        print(U_list_N_matrices)
        
    U_full = reduce(np.matmul, U_list_N_matrices) # here adjoint the full time evolution 
        # here U_list_N_matrices = [e^{-iH(T-dt)}]
    return U_full   

#%%
#TEST : if U_full is unitary? 
U_full_6TS =  U_full_discretisation(np.pi/2, np.pi/5, 6) # these kx, ky can be changed 
#adjoint of U_full_6TS(U_adgga)
U_full_6TS_adj = U_full_6TS.conj().T
print("\nU† =", U_full_6TS_adj)
# check if U is really 'unitary'
check_unitariness = U_full_6TS @ U_full_6TS_adj
print('unitary?', check_unitariness) # -- we should see this is always identity(U is indeed Unitary)
#ANS： Yes. U_full (U(k, T)) is unitary
#But:  It doesn't eqaul to identity at all; so due to RLBL article III-A, it needs to be modified 
#Q: is U diagonalised? If so it would be much more better! Since we don't need to include 'V' as in U_dia = V U V^-1 subsequently
#%% the logarithm is extracting the eigenphases of U when U is unitary 
# H_eff = template code -- 
def heff_from_U(U, T, epsilon, tol=1e-10):
    """
    Compute H_eff = (i/T) log_epsilon(U)
    with the logarithm branch cut along exp(-i * epsilon * T).

    Parameters
    ----------
    U : (N,N) complex ndarray
        Floquet/unitary operator.
    T : float
        Period.
    epsilon : float
        Quasienergy value defining the branch cut.
    tol : float
        Numerical tolerance for unitarity check.

    Returns
    -------
    Heff : (N,N) complex ndarray
        Effective Hamiltonian.
    evals : (N,) complex ndarray
        Eigenvalues of U.
    quasienergies : (N,) float ndarray
        Branch-chosen quasienergies used in Heff.
    """
    U = np.asarray(U, dtype=complex)
    N = U.shape[0]

    # Optional: check unitarity
    I = np.eye(N, dtype=complex)
    if not np.allclose(U.conj().T @ U, I, atol=tol):
        print("Warning: U is not exactly unitary within tolerance.")

    # Diagonalize U
    evals, V = np.linalg.eig(U)

    # For a unitary matrix, eigenvalues should lie on the unit circle:
    # evals_n = exp(i * alpha_n), where alpha_n is an angle.
    alpha = np.angle(evals)   # principal angles in (-pi, pi]

    # Branch cut direction is exp(-i * epsilon * T)
    cut = -epsilon * T

    # Put alpha on the interval (cut - 2*pi, cut]
    alpha_branch = ((alpha - cut) % (2 * np.pi)) + cut
    alpha_branch = np.where(alpha_branch > cut, alpha_branch - 2 * np.pi, alpha_branch)

    # Since log(lambda_n) = i * alpha_branch_n for |lambda_n|=1,
    # Heff eigenvalues are:
    quasienergies = -alpha_branch / T

    # Reconstruct Heff
    V_inv = np.linalg.inv(V)
    Heff = V @ np.diag(quasienergies) @ V_inv

    # Clean up tiny non-Hermitian numerical noise
    Heff = 0.5 * (Heff + Heff.conj().T)

    return Heff, evals, quasienergies

#%% write my own Heff code

# find V as similarity trans matrix. U_diagonalised = V U V^-1 = (diagonal: eigenval = e^(-j Phi_1)... e^(-j n))

#U_full_bulk is the 2x2 full time evolution opt difined by different kx, ky, but not t as in they are U_full = U (k_vec, T)

# later maybe can embed the U_full_bulk(T) = [e ^ (-j *H(k,t_n_last) *delta_t)] ... [e ^ (-j *H(k,t_n_first) *delta_t)] 
# where H2x2 calculated by function H_matrix_2by2, with input k_x, k_y, t_n, a_0(as fixed lattice length)
# from INPUT: kx, ky, t 

# Q1: what is the upper and lower limitation for the winding number integeral(equation (4)) when U is general case?
# Q2: does general case U_epsilon(k, t) we setted, use the same definition as U(k, t) as in t_vor article, eqn 21 (U(...t) does mean to multiply all H up to t)
#     - Ans: the first half in one period (T) have U(k, 2t) used same; hence should use func H(k,t) >> times up to U(t)
#            double check the maths behind why we can (ex:U(T) = U6U5...U1) break U(t)= U(t)...U(1)
#     - Ans: the second half in one period (T) = V_epsilon, which is from U(T)'s corresponding effective Ham Heff 


#LATER ADAPT U_full to: when calculating Heff, INPUT needed now: U(T); but we can also break U(full) into that we define H(k, t) and t descritisation(summed over)
# INPUT will be: kx, ky, num of discritisation(for T/n = delta t, H(t) in expoential of U_n...U_1 is H( T/n * i)); subQ: can we introduce time vortex in bulk k space Hamiltonian? 
#                                                                                                                        If we cannot, why? 

def Calculate_Heff_from_Ufull (U_full_bulk, epsilon): 
    
    eigvals_Ufull, eigvecs_Ufull = np.linalg.eig(U_full_bulk) #this matrix is diagonalised; what is the order of eigenval arrangement on diagonal?
    # confirm one to one correspondence (somehow? you can write a bit on paper) of the eigvals and eigvecs which compose V 
    U_dia_full = np.diag(eigvals_Ufull) # the way we define eigval is same as how the eigvals returned from np.linalg.eig 
    
    V = eigvecs_Ufull # need to be checked that if they are natrually returned in columns 
    
    
    # Assume V Heff V^-1 from log (with defined branch cut) (U_diag_full) for log with chosen branchcut; now is trying to write branch cut 'into' log
    #branch cut only kicks in for how we extract Heff from U_diag -- next step is: How? 
    
    #Assume the log(matrix) is log(with branchcut) of each of the elements on diag
    log_U_conv_branchcut =  - np.angle(U_dia_full) # np.angle doesn't return phi, but return minus phi (-phi); we return phi for fitting branchcut change
    # if  - np.angle = phi_0 ( a specific )
    # do it eigval-wise 
    
    phi_0_adjusted = []
    for eigval_Ufull in eigvals_Ufull:
        
        phi = -np.angle(eigval_Ufull) # range: (-pi, pi]
        #as 'the angle' in log function that defines the branch cut for the function(which shows together with minus sign in the function )
#----------------- this is the old way, for only epsilon T larger then 0 (< pi)
#        if phi < epsilon * T:     #this part is for phi in range(returned by np.angle) 0 to pi
#            value = phi + 2 * np.pi

#        else:                        #this part is for phi in range(returned by np.angle) -pi to 0
#            if 0 < phi <= np.pi:
#                value = phi
#            elif -np.pi < phi < 0:
#                value = phi + 2 * np.pi
#            else:
#                value = phi
#------------------ this is the old way, for only epsilon T larger then 0 (< pi)  
    # generally below:
        if 0 <= phi <= np.pi:
            phi_0 = phi
        else:
            phi_0 = phi + 2 * np.pi #phi_0 is from 0 to 2pi, including 0 but not 2pi 
        
        # here from this line, in the above indentation, will phi_0 keep going to the next if function and return appendable value ? 
        if  epsilon* T > 0:
            if phi_0 <= epsilon * T:
                phi_0_adjusted.append(phi_0 + 2 * np.pi)
            if phi_0 > epsilon * T :
                phi_0_adjusted.append(phi_0)
        elif epsilon * T <= 0:
            if phi_0 >= epsilon * T + 2 * np.pi :
                phi_0_adjusted.append(phi_0 - 2 * np.pi)
            if phi_0 < epsilon * T + 2 * np.pi:
                phi_0_adjusted.append(phi_0)
    
    phi_newbranch_array = np.array(phi_0_adjusted) # this array has same order with U_diag (as in the arg's corresponding elements, or saying: eigval)
    # this is the branchcut
    H_eff = (1j * (np.diag(phi_newbranch_array)))/T
    
    return H_eff


