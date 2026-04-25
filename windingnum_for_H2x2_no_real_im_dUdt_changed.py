# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 23:35:31 2026

@author: DELL
"""

import numpy as np
#import sympy as sp
import scipy.integrate as spint
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

#print (J_tilte_x(0))
#print (J_tilte_x(0.5))

#print (J_tilte_y(1/3))
#print (J_tilte_x(5/6))
# here is clear that each pulse is defined as start point(of driving) = Ja, end point (of driving) = 0
#%%
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
        
    
        t_i = T - t_desicretisation * (i + 1) + t_desicretisation/2 # here shift sample t_i to midpoint 
                                                                    # i = 0 (leftest in muti) t_i = T - dt + dt/2
                                                                    # i = N -1 t_i = T - (T/N) * N + T/2N
                                                                    # this will be safe for any case
        H_i = H_matrix_2by2(k_x, k_y, t_i, a_0)
#expm(-1j* H_i_test* delta_t_test) 
        
        U_i = expm(-1j* H_i* (t_desicretisation)) #!!!!! # here using exponential form construct each of U_1 to U_6
        
        U_list_N_matrices.append(U_i)
        
#        print(U_list_N_matrices)
        
    U_full = reduce(np.matmul, U_list_N_matrices) # here adjoint the full time evolution 
        # here U_list_N_matrices = [e^{-iH(T-dt)}]
    return U_full   
# this will be for U(k, t)    

#%% this is for not U(k, T)/U(r, T), but U(k/r, t)
def U_t_discretization(k_x, k_y, num_time_stages, t): # which means along lhs to rhs, from (1) earlier/smaller t to (6) later/larger t

    U_list_N_matrices = []
    t_discretization= t /num_time_stages
#    time_stage = num_time_stages
    
    for i in range (num_time_stages):
        
    
        t_i = t - t_discretization * (i + 1)  + t_discretization/2   #t_discretization = delta_t here, **take midpoint  
        # here we shifted to midpoint 
        # note: here t_discretization << T/6
        # because if t is any value, and no classification of t range is done ..|.....|.....|.....| 
        # t/ N_dis doesn't need to not across the '|' and blur the change from Hi to Hi+1 at the boundaries 
        # double check here: i = 0 (leftest), t_i = t - delta_t
     #                       i = num - 1(rightest) , t_i = t - num * dis = t - t = 0 
        H_i = H_matrix_2by2(k_x, k_y, t_i, a_0)
#expm(-1j* H_i_test* delta_t_test) 
        
        U_i = expm(-1j* H_i* (t_discretization)) #!!!!! # here using exponential form construct each of U_1 to U_6
        
        U_list_N_matrices.append(U_i)
        
#        print(U_list_N_matrices)
        
    U_full = reduce(np.matmul, U_list_N_matrices) # here adjoint the full time evolution 
        # here U_list_N_matrices = [e^{-iH(T-dt)}]
    return U_full   
# this will be for U(k, t)
#%%
#test_reduce()
#test_x = np.array([[1j +1, 3j -2], [1, 3]])
#test_y = np.array([[1j + 0.5, 0.2j + 5], [5, 1]])

#list_test = [test_x, test_y]
#list fruit = [apple, pear, peach]
#U_test = reduce(np.matmul, list_test) 
#print(U_test)
# all verified 
#%% here is a simple version can be applied to U_full_k as in 6 constant stages

def U_full_descending(k_x, k_y, num_time_stages): # U_full =  U6 U5 U4 U3 U2 U1 as in U1 is from 0 to T/6

# difference from discretisation func: H(t) sample take t in the midpoint of H const stage 
    U_list_N_matrices = []
    t_desicretisation= T /num_time_stages
#    time_stage = num_time_stages
    
    for i in range (num_time_stages):
        
    
        t_i = T - t_desicretisation * (i + 1/2) # from last one to first one
    
        H_i = H_matrix_2by2(k_x, k_y, t_i, a_0)
#expm(-1j* H_i_test* delta_t_test)
        
        U_i = expm(-1j* H_i* (t_desicretisation)) # H_0 will be H (t = T-(T/6)*1/2), with is 6th stage
        
        U_list_N_matrices.append(U_i)
        
    U_full = reduce(np.matmul, U_list_N_matrices) # here adjoint the full time evolution 
        # here U_list_N_matrices = [e^{-iH(T-dt)}]
    return U_full 

#%%
#TEST : if U_full is unitary? 
U_full_6TS_wrong =  U_full_discretisation(np.pi/2, np.pi/5, 6) # these kx, ky can be changed 
# shouldn't risk to use discretisation func to calculate the 6 const H partern, for that the start/end point of a pulse 
# might be defined as 0 or Ja, uncertain

#adjoint of U_full_6TS(U_adgga)
U_full_6TS_wrong_adj = U_full_6TS_wrong.conj().T
print("\nU† =", U_full_6TS_wrong_adj)
# check if U is really 'unitary'
check_unitariness = U_full_6TS_wrong @ U_full_6TS_wrong_adj
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


#%% Rewrite the function above: 
#   INPUT: U_full_bulk = U(k, T) epislon -- which is used to define the branchcut

def Heff_from_Ufull (U_full_bulk, epsilon): # U_full_bulk = U(k, T)
    
    eigvals_Ufull, eigvecs_Ufull = np.linalg.eig(U_full_bulk)  #Each eigval_alpha corresponding to each eigenvector |phi_alpha>
    
    U_dia_full = np.diag(eigvals_Ufull) # diagonalised U, now each element on the diagonal is lamda_alpha = eigval_alpha corresponding to eigenvector |phi_alpha>
    
    V = eigvecs_Ufull # need to be checked that if they are natrually returned in columns  # ***VERIFICATION 1
    
    PHI_dia = np.angle(U_dia_full) #np.angle works elementwise on each entry and returns the phase angle of each complex number, in the range (-pi, pi]
    #take angle as an array maintian the order of alpha(s) perfectly
    #**VERIFICATION2**: U_diag = np.exp(1j*PHI_diag) instead of **elementwise** definition of phi_alpha = - PHI_alpha

    #Q: expm(1j*M) is matrix exponential -- how does it work diffrently from elementwise exp? 
    PHI_alpha_array = np.diagonal(PHI_dia) # branch range: (-pi, pi]



    PHI_modified_list = []
    
    # possible modified branchrange: 
        #epsilonT range: [-pi, pi] PHI_modified_min > -epsilonT_min - 2pi = -3pi; PHI_modified_max <= -espilonT_max = -(-pi) = pi 
    # hence (rough) PHI_modified range (with above information): (-3pi, pi]
    # we only test PHI_alpha minus 2npi, because only the smallest PHI, (-pi), + 2pi will be in the range of PHI_modified;
    # That need epsilon T to be -pi, hence - epsilon * T = -(-pi) = pi; - epsilon * T - 2 * np.pi = -pi, PHI_alpha = -pi is out of range--> modified to pi
    # but this is a special case, 
    
    for PHI_alpha in PHI_alpha_array:
        
        if  - epsilon * T - 2 * np.pi < PHI_alpha <= - epsilon * T:   # here used (] for PHI, as we take -epsilonT but not  - epsilonT - 2* pi
            
            PHI_modified = PHI_alpha
        
        elif  - epsilon * T - 2 * np.pi <  PHI_alpha - 2 * np.pi <=  - epsilon * T:
            
            PHI_modified = PHI_alpha - 2 * np.pi
        
        elif - epsilon * T - 2 * np.pi <  PHI_alpha - 4 * np.pi <=  - epsilon * T:
            
            PHI_modified = PHI_alpha - 4 * np.pi
        
        elif - epsilon * T - 2 * np.pi <  PHI_alpha + 2 * np.pi <=  - epsilon * T:
            
            PHI_modified  =  PHI_alpha + 2 * np.pi
            print(PHI_alpha, PHI_modified) # this step is to check if the guess, only case here should be -pi + 2pi = pi, is correct
            
        else:
            raise ValueError(f"No condition matched for PHI_alpha = {PHI_alpha}")

        
        # also remember: in one branch each PHI_alpha_modified only has one value choice
        PHI_modified_list.append(PHI_modified)
        
    phi_alpha_array = - np.array(PHI_modified_list)
    
    H_eff_dia = (1/ T) * np.diag(phi_alpha_array)

    H_eff = V @ H_eff_dia @ np.linalg.inv(V)

    return H_eff
#VERIFICATION 2: within the code above 
#VERIFICATION 3: pen and paper: what if we go from phi_diag = - np.angle (U_dia_full)? Done
#VERIFICATION 4: input a real matrix Phi to calculate Heff, see if there's any error.Also compare both function Heff_from_Ufull and Heff_from_Ufull 2 output
#VERIFICATION 5: test if V is as what expected

#%%
def VERIFICATION2(U_full_bulk):
    
    eigvals_Ufull, eigvecs_Ufull = np.linalg.eig(U_full_bulk)  #Each eigval_alpha corresponding to each eigenvector |phi_alpha>
    
    U_dia_full = np.diag(eigvals_Ufull) # diagonalised U, now each element on the diagonal is lamda_alpha = eigval_alpha corresponding to eigenvector |phi_alpha>

    PHI_diag = np.angle(U_dia_full)
    phi_diag = - np.angle(U_dia_full)

    
    if np.allclose(expm(1j * PHI_diag), U_dia_full):
        
        if np.allclose(expm(-1j * phi_diag), U_dia_full):
            
         print ('VERIFIED')
        else:
            
            raise ValueError('error TEST2')
    
    else:
        raise ValueError('error TEST1')
        
    return 

F = VERIFICATION2(U_full_descending((2*np.pi)/8, (2*np.pi)/8, 6))
# VERIFICATION2: verified 

#%%

def Heff_from_Ufull_2 (U_full_bulk, epsilon): # Heff_from_Ufull_2 has no PHI, only phi_alpha 
    
    eigvals_Ufull, eigvecs_Ufull = np.linalg.eig(U_full_bulk)  #Each eigval_alpha corresponding to each eigenvector |phi_alpha>
    
    U_dia_full = np.diag(eigvals_Ufull) # diagonalised U, now each element on the diagonal is lamda_alpha = eigval_alpha corresponding to eigenvector |phi_alpha>
    
    V = eigvecs_Ufull # need to be checked that if they are natrually returned in columns  # ***VERIFICATION 1
    
    phi_dia = - np.angle(U_dia_full) #[-pi, pi)
    #take angle as an array maintian the order of alpha(s) perfectly

    phi_alpha_array = np.diagonal(phi_dia) # branch range: (-pi, pi]
    
    phi_modified_list = []
    
    # possible modified branchrange: 
        #epsilonT range: [-pi, pi] PHI_modified_min > -epsilonT_min - 2pi = -3pi; PHI_modified_max <= -espilonT_max = -(-pi) = pi 
    # hence (rough) PHI_modified range (with above information): (-3pi, pi]
    # we only test PHI_alpha minus 2npi, because only the smallest PHI, (-pi), + 2pi will be in the range of PHI_modified;
    # That need epsilon T to be -pi, hence - epsilon * T = -(-pi) = pi; - epsilon * T - 2 * np.pi = -pi, PHI_alpha = -pi is out of range--> modified to pi
    # but this is a special case, 
    
    for phi_alpha in phi_alpha_array:
        
        if  epsilon * T <= phi_alpha <  epsilon * T + 2 * np.pi :   # here used (] for PHI, as we take -epsilonT but not  - epsilonT - 2* pi
            
            phi_modified = phi_alpha
        
        elif  epsilon * T <= phi_alpha + 2 * np.pi <  epsilon * T + 2 * np.pi:
            
            phi_modified = phi_alpha + 2 * np.pi
        
        elif  epsilon * T <= phi_alpha + 4 * np.pi <  epsilon * T + 2 * np.pi:
            
            phi_modified = phi_alpha + 4 * np.pi
        
        elif epsilon * T <= phi_alpha - 2 * np.pi <  epsilon * T + 2 * np.pi:
            
            phi_modified = phi_alpha - 2 * np.pi
            print(phi_alpha, phi_modified) # this step is to check if the guess, only case here should be -pi + 2pi = pi, is correct
            
        else:
            raise ValueError(f"No condition matched for phi_alpha = {phi_alpha}")

        
        # also remember: in one branch each PHI_alpha_modified only has one value choice
        phi_modified_list.append(phi_modified)

    
    H_eff_dia = (1/ T) * np.diag(np.array(phi_modified_list))

    H_eff = V @ H_eff_dia @ np.linalg.inv(V) # A *B is elementwise multiplication but not matrix multiplication 

    return H_eff

#%%

# Which U equation is recalled? U_full_descending, which is only 6 stages( num_of_stages = 6, fixed, hence t_discretisation = T/6 which is large)
#                                                                          hence here the H_i is taken as the mid-point of each t_discre interval
#                               U_discretisation input with large N, very small t_discretisation; hence H_i is taken as lower edge point (0), delta_t, ...
def Heff_from_kxky (k_x, k_y, num_time_stages, epsilon):
    U_kxky = U_full_descending(k_x, k_y, num_time_stages)   # U(k, T)
    Heff = Heff_from_Ufull (U_kxky, epsilon)
    return Heff

def Heff_from_kxky_2(k_x, k_y, num_time_stages, epsilon):
    U_kxky = U_full_descending(k_x, k_y, num_time_stages)   # U(k, T)
    Heff = Heff_from_Ufull_2 (U_kxky, epsilon)
    return Heff

# here for bulk as in infinite amount of sites, or infinite length of system, no edge, N --> inf; k_x, k_y spacing -->0

def VERIFICATION_3(k_x, k_y, epsilon):
    Heff_1 = Heff_from_kxky (k_x, k_y, 6, epsilon)
    Heff_2 = Heff_from_kxky_2 (k_x, k_y, 6, epsilon)
    
    if np.allclose(Heff_1, Heff_2):

         print ('VERIFIED')
    else:
            
        raise ValueError('error TEST')
        
    return 

#%% can be commented out
Heff_test = Heff_from_kxky ((2*np.pi)/8, (2*np.pi)/7, 6, np.pi/2)

U_full_test = U_full_descending((2*np.pi)/8, (2*np.pi)/7, 6)

Heff_test_exp = expm(-1j * Heff_test * T)
TEST1 = np.allclose(Heff_test_exp, U_full_test)
print(TEST1)

#%% can be commented out
Heff_test_2 = Heff_from_kxky_2 ((2*np.pi)/8, (2*np.pi)/7, 6, np.pi/2)

Heff_test_exp_2 = expm(-1j * Heff_test_2 * T)
TEST2 = np.allclose(Heff_test_exp_2, U_full_test)
print(TEST2)
#%% can be commented out
Heff_verify = VERIFICATION_3((2*np.pi)/8, (2*np.pi)/7, 0)
    
#%%
eigvals_Ufull, eigvecs_Ufull = np.linalg.eig(U_full_test)  #Each eigval_alpha corresponding to each eigenvector |phi_alpha>
    
U_dia_full = np.diag(eigvals_Ufull) # diagonalised U, now each element on the diagonal is lamda_alpha = eigval_alpha corresponding to eigenvector |phi_alpha>
    
V = eigvecs_Ufull

U_return = V @ U_dia_full @ np.linalg.inv(V) #wrong: U_diagonalised = V U V^-1; true:  A=VDV−1
TEST3 = np.allclose(U_return, U_full_test)
print(TEST3)

#%%

#equation 3 as in the definition of U(k, t) Done; U(k, t) = U_k (t, 0)

#I equa 8 finish   

#II integral   

#%% Equation 8

# the U_e function below can only calculate U(t) with t within [0, T]
# only with H(2x2) for kx, ky, without time vor
# -- otherwise : H_i = H2NyNx_xyOBC_time_vor_core (N_x, N_y, t_i, a_0, x_0, y_0) instead of H_matrix_2by2

def U_epsilon_H2by2 (kx, ky, t, N_dis, epsilon): # N_dis stand for time discretisation we choose to calculate U(k,t) here 
                                 #only_bulk_U_involved will be =True or = False
                                 #Physically 'only_bulk_U_involved' means no time vortices 

    U_full_for_Heff = U_full_discretisation(kx, ky, num_time_stages = 6) 
    
    #**if any change on H, H isn't valid for time discretization only 6 (piecewise) anymore, above should be changed 
    
    
    if 0<= t <= T/2:
        
        t_epsilon = 2 * t
    # alternative way record(little notebook p79): we can just not call any U function before, but classify t_epsilon into 6 stages, for each stage we have a certain parttern of 
        U_epsilon =  U_t_discretization(kx, ky, N_dis, t_epsilon)
                                       # here was set as 6 in the could we send to Mimi intially, wrong!
                    #U_t_discretization(k_x, k_y, num_time_stages, t) - these are meaning of variables 
       #U_e (k, t) = U(k, 2t) = U_latest......U_earliest which should be calculated in numerical(infinitesimal) pieces 
       
    elif T/2 <= t <= T:  # here this must be elif, otherwise all t belongs to the first 'branch' 
                         # will finally go to else, which raise error interrupt the function
        
        t_epsilon = 2 *T - 2 * t # 2nd def
        
        H_eff_epsilon = Heff_from_Ufull_2(U_full_for_Heff, epsilon) 
        #Here this H_eff_epsilon has specific k_vec as U_full_for_Heff
        
        U_epsilon =  expm(-1j* H_eff_epsilon * t_epsilon ) # H_eff_epsilon refer to above, t_epsilon as 2nd def
        
    else:
        raise ValueError(f"No condition matched for U(k_vec, t) with time  = {t}")
    return U_epsilon  

# t >T/2 can be considered as diagonalised matrix, if use diagonalised H_eff_epsilon (?)
# t< T/2 has depends on U_t_discretization = U_N-1 ... U_0
# if U_n = diagonalised? H_n = H(2x2), check on eqn (28) NO **


        # U = e ^-i *Heff*t
        # Q: which e, expm or e elementwise we should use? 
        # Ans: expm. expm doesn't equal to e - elementwise for non-diagonalised matrices
        # Heff calculated in eqn(9) is not necessarily diagonalised. 
        #We only used the diagonal elements to set branch cut, will log is used elementwised
        # ** for function of operatorsf(A), power series is f(A) = I+A+2!A2​+3!A3​+; 
        #U(T) = e ^-i Heff T os f(Heff), Heff = scalar * log (U(T)) = function of U, but practically we used U_diagonal to calculate then transfer back to not diagonal matrix 
        # For that epsilon_alpha is set to be Heff's eigenval on each phi_alpha with the branchcut definition
        # ( that's how we 'design' Heff)
# equation 8 finished 
#How to refer the function above: 
#print( U_test_language(only_bulk_U_involved=True))
#print( U_test_language(only_bulk_U_involved=False))
#%% calculate integeral - Q: what U we will input in integral?  Ans: W[U_epsilon]
# Therefore, the number of edge modes of U at quasienergy " must be the same as 
#the number of edge modes of U" at quasienergy =T. The latter quantity is then given by W[U_epsilon] 
# double check about U(0) and U(T) for any k = identity? 
# theoretically 1. 'The first property is that U_epsilon has a trivial Floquet operator:** U(T) = I** all k'
#               2. Does U(0) = U(T)? yes for U_espilon, no for U original; U_epsilon (0) = U(0,0) = I
UT_test_I = U_epsilon_H2by2 (np.pi*2, np.pi/3, T, 50, np.pi/15) #only_bulk_U_involved=True
print(UT_test_I)
#U0_double_check = U_epsilon(np.pi*2, np.pi/3, 0, 50, np.pi/15, only_bulk_U_involved=True)
#print(U0_double_check) # checks out

# how to do tr, how to do /dt or /dk for matrices;

#%% 
def f(x, t, omega = 2.6):
    return(t * np.sin(x + t * omega))

# Here is an example of how to integrate f w.r.t t at x = -1.0

res = spint.quad(lambda t : f(-1.0, t), 0.0, 1.0)
print(res)


#%%
#function needed when calculate integral, **'integral function'
def quad3d(func, a, b, args = ()):
    # calculates the 3d quad integral of func(x, y, z), where
    # a = (x_min, y_min, z_min), b = (x_max, y_max, z_max)
    
    # args is a tuple of all the remaining arguments of func
    
    # Returns value and numerical error
    
    x_min, y_min, z_min = a
    x_max, y_max, z_max = b
    
    if not isinstance(args, tuple):
        # A single argument. We wrap it into a tuple
        args = (args,)
    
    return(spint.tplquad(func, z_min, z_max, y_min, y_max, x_min, x_max, args))

def wassup(x, y, z, lol, lmao):
    return(x+lol*y+lmao*z)

print(quad3d(wassup, (0,0,0), (1,1,1), (2, -3)))

#%%
# Example so you see how it works

def f3d(t, x, y, omega):
    return(x * np.sin(y + omega * t))

# We will integrate for t = (-1.0, 1.0), x = (0.0, 1.0), y = (5.0, 10.0), with omega = 3.0

res, err = quad3d(f3d, (-1.0, 0.0, 5.0), (1.0, 1.0, 10.0), (3.0))
#                       lb t x y          ub t x y 
print(f"Integral is {res:0.5f} pm {err}")

#%%
def matrix_derivative(M, t, dt = 1e-6):
    return((M(t + dt) - M(t - dt)) / (2 * dt))

#matrix derivative: 
    
def lol(t):
    return(np.array([
            [np.cos(t), -np.sin(t)],
            [np.sin(t), np.cos(t)]
        ]))

print(matrix_derivative(lol, 0))


#%%
# construct the thing inside integral 
# N_dis is just as t_smapling in fig 5ab prodece, = 50
 
def Tr_integrated_MOCK(kx, ky, t, N_dis, epsilon):    # only_bulk_U_involved = True -- this will be encoded in U 
    
    Mock_m = np.array([[kx, ky],[N_dis * kx, epsilon* ky ]])
    Mock_trace  = np.trace(Mock_m)
    return Mock_trace 

#BZ range coherent with before: (k_x_TEST = generate_k_x_list(0, 2*np.pi, N_x_TEST))
# as the func in RLBL article: against /  t=0, T /kx = 0, 2pi  / ky = 0, 2pi Tr(...) dt dkx dky
W_res, W_err = quad3d(Tr_integrated_MOCK, ( 0.0, 0.0, 0.0), (2*np.pi, 2*np.pi, 1), (50, np.pi))  # here epsilon = pi
#                      func           lb  x y t         ub x y t
print(f"Integral is {res:0.5f} pm {err}")

#Q : _quad
#    return _quadpack._qagse(func,a,b,args,full_output,epsabs,epsrel,limit)

#TypeError: only size-1 arrays can be converted to Python scalars
#%% Write out Trace func within: 
def matrix_partial_dt(M, kx, ky, t, N_dis, epsilon, dt = 1e-6): # after M and before dt are just arg of M
                                                        # U = U_epsilon_H2by2 (kx, ky, t, N_dis, epsilon) 
                                                        # for example partial against t: kx, ky, N_dis, epsilon are fixed
    dM_dt = ((M(kx, ky,t + dt , N_dis, epsilon) - M( kx, ky,  t - dt, N_dis, epsilon)) / (2 * dt))
    # at kx, ky, t, N_dis, epsilon 
    return dM_dt
    # this function establishes but it will calculate twice for matrix U(t, 0), one with d + dt and one with d -dt 
    # which is 1. less precise 
              #2. call U_epsilon function once, you do once discretization of t which is tricky ( N_dis = 20)
def matrix_partial_dkx(M, kx, ky, t, N_dis, epsilon, dkx = 1e-6): # after M and before dt are just arg of M
                                                        # U = U_epsilon_H2by2 (kx, ky, t, N_dis, epsilon) 
                                                        # for example partial against t: kx, ky, N_dis, epsilon are fixed
    dM_dkx = ((M( kx + dkx, ky, t, N_dis, epsilon) - M( kx - dkx, ky, t,  N_dis, epsilon)) / (2 * dkx))
    # at kx, ky, t, N_dis, epsilon 
    return dM_dkx


def matrix_partial_dky(M, kx, ky, t, N_dis, epsilon, dky = 1e-6): # after M and before dt are just arg of M
                                                        # U = U_epsilon_H2by2 (kx, ky, t, N_dis, epsilon) 
                                                        # for example partial against t: kx, ky, N_dis, epsilon are fixed
    dM_dky = ((M(kx, ky +dky, t, N_dis, epsilon) - M(kx , ky - dky, t,  N_dis, epsilon)) / (2 * dky))
    # at kx, ky, t, N_dis, epsilon 
    return dM_dky
    # this returns dM/dt at t, but also kept kx, ky, N_dis, epsilon as specific input     

# this returns dM/dt at t, but also kept kx, ky, N_dis, epsilon as specific input 



 
def Tr_integrated (kx, ky, t, N_dis, epsilon):    # N_dis and epsilon will be fixed like rhat in mock
                                                  # only kx, ky, t are what intergated against 
    U = U_epsilon_H2by2 (kx, ky, t, N_dis, epsilon)
    U_inverse  = np.linalg.inv(U)
    

#   which is one var ; multiple var: 
    #dU_dt = matrix_partial_dt(U_epsilon_H2by2, kx, ky, t, N_dis, epsilon)
    
    
    
    #--------------below defined a func of dU/dt in more symbolic way ----------------------
    #**if any problem occured, this dU/dt below can be double check for 
    # t = T/2, t = 0 and t = T, for both U_epsilon def1
                                #and  U_epsilon def 2
    
    U_full_for_Heff = U_full_discretisation(kx, ky, num_time_stages = 6)
    # this is for the condition t> T/2, goes to second definition 
    if 0<= t <= T/2:
        
        dU_dt = -2j * H_matrix_2by2(kx, ky, 2*t, a_0) @ U_epsilon_H2by2(kx, ky, t, N_dis, epsilon)
        
    elif T/2 <= t <= T:
        
        Heff_for_div = Heff_from_Ufull_2(U_full_for_Heff, epsilon)
        dU_dt = -1j * Heff_for_div @ expm(-1j* Heff_for_div * (2*T - 2*t))  * (-2)   
                                                                          
   #--------------below defined a func of dU/dt in more symbolic way ----------------------
   
   
    dU_dkx = matrix_partial_dkx(U_epsilon_H2by2, kx, ky, t, N_dis, epsilon)
    dU_dky  = matrix_partial_dky(U_epsilon_H2by2, kx, ky, t, N_dis, epsilon)
    
    part_1  = U_inverse @  dU_dt
    part_comm_1 =  U_inverse @  dU_dkx
    part_comm_2 =  U_inverse @  dU_dky
    comm =  part_comm_1 @ part_comm_2 -  part_comm_2 @ part_comm_1
    Trace  = np.trace(part_1 @ comm)
    # aim for now: 5:45 : adapt matrix div func to be partial div 
    # anyway finish this script of the trace func which is integrated against 
    # mock func 
    return Trace 


#def Tr_integrated_MOCK(kx, ky, t, N_dis, epsilon):    # only_bulk_U_involved = True -- this will be encoded in U 
    
#    Mock_m = np.array([[kx, ky],[N_dis * kx, epsilon* ky ]])
#    Mock_trace  = np.trace(Mock_m)
#    return Mock_trace 

W_res, W_err = quad3d(Tr_integrated, ( 0.0, 0.0, 0.0), (2*np.pi, 2*np.pi, 1), (20, np.pi))  # here epsilon = pi
#                      func           lb t x y          ub t x y 
print(f"Integral is {res:0.5f} pm {err}")
# Q: can we skip the non-digonal matrix U somehow? Heff has a linear version
# cannot; For that from beginning H2x2 isn't digonal matrix 

# Q: do we need sympy to calculate the divU? The point is that each U_epsilon_H2by2 (kx, ky, t, N_dis, epsilon) 
# will call the U_t_discretization which slices U into U_0 ... U_(N-1)
# is numerical and too complicated to be div (as a series of multiplication of e^)
# Ans: We do U_div numerically, not symbolically 



