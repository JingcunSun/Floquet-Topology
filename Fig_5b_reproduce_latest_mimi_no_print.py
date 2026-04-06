

# few general summary of the codes for my CV
# all these codes are working with the Hamiltonian as well as eigen vectors, eigen values of it.
# also include FT (real space, k space)



import numpy as np
from functools import reduce
from scipy.linalg import expm
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from collections import Counter
import pandas as pd
from pathlib import Path


#from collections import Counter
#%%

 # this function if for checking spectrum of two Hamiltonians -(for example real and k space)


# Define function J_tilde_x(t), J_tilde_y, J_tilde_z, they are all step functions with same step length of t. 
    # -Physical meaning：J_x,y,z (u_{I, J})is the strength of (nearest neighbour) hopping. 
    # -Here when J_tilte_a(t) is a function of time, the 'tilte' also encompass a pair (?) of u_{I, J} choice. 


#J_tilte_a(t) are all defined in between 0 and T, and t is independent variable - here you can ask later about 
#-- how to distinguish t and delta_t, J_a



#def J_tilte_a(t, t_start, t_end, J_a): #delta_t = t_end - t_start
# parameter of I_tilte_a: start point t_a (a for x, y, z), pulse length delta_T, (l, T), amp J_a  
  
# this J_tilte_a function is based on the paragraph above eqn(16)
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



#J_x = J_a
#J_y = J_a
#J_z = J_a

#t_a_array = [t_x, t_y, t_z]
#adjusted by the caption under fig_1:(we unify the t_a, J_a, l, T, delta_t) 
#laeve the only parameter adjustable in J_tilte_a is t


#%%

#write six different Hamiltonians -- Figure out how to order H1 to H6, U1...U6, then finish writing the U_full
#H_1 to H_6, here arrange them from latest time to earliest time. 

#J_tilte_x = J_tilte_a(t, t_x, delta_t, l, T, J_x)
#J_tilte_y = J_tilte_a(t, t_y, delta_t, l, T, J_x)
#J_tilte_z = J_tilte_a(t, t_z, delta_t, l, T, J_z)

# these J needed to be written at their corresponding time. 

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

t_test = np.linspace(-2, 1, 100) 
# function of J_a need to be adjusted as in that t can be out of range of one period(0,T)


#%% some test try figure properties of J func
t_specific_time_point_test = np.pi - 3

call_func_J_tilte_x = J_tilte_z(t_specific_time_point_test)


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

#%%  TIME VORTEX
#for induced time vortex H: in H, in (28), we will have the J_tilte_a(t - T* arctan[((H_m_n)y + (H_m_n)y'))/((H_m_n)x + (H_m_n)x'))/2pi] 
# wait? T/2pi is just omega -- no, wrong. The omega part is still denote as epsilon I think as in epsilon_n T = from -pi to pi, 
# also 'density of state' is abs|psi_n| ^2 which can be labelled by epsilon_n - phi_n - psi_n(which is the state)

# here (H_m_n)x(') or y(') denotes the relationship between {the entrices in the Hamiltonian mth row, nth column (which corresponding to x - x' and y- y')}
# and {the x, y, x', y' unit cells they corresponding to}
# Or saying in another way H now is corresponding to both x, x', y, y', and the H_kx_ky(2x2) as in terms of summations of kxky 
#                                                                                      J_tilte denoted by x, y and x' y' 

#'dog'
#%%
#time discretisation:
  
def U_full_discretisation(k_x, k_y, num_time_stages): # which means along lhs to rhs, from (1) earlier/smaller t to (6) later/larger t
    
    H_list_6_matrices =[]
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
# ----- U FROM K BASIS MATRIX H 2X2 (U(T)-H(K))----
#Here is the 2by2 time evolution matrix: (1 to 6) staged time evolution respectively & full time evolution

def U_full_ascending(k_x, k_y, num_time_stages): # which means along lhs to rhs, from (1) earlier/smaller t to (6) later/larger t
    
#    H_list_6_matrices =[]
    U_list_N_matrices =[]
    
#    time_stage = num_time_stages
    
    for i in range (num_time_stages):
    
        t_i = (T/num_time_stages) * i + T/(num_time_stages * 2)
    
        H_i = H_matrix_2by2(k_x, k_y, t_i, a_0)
#expm(-1j* H_i_test* delta_t_test)
        
        U_i = expm(-1j* H_i* (T/num_time_stages)) #!!!!! # here using exponential form construct each of U_1 to U_6
        
        U_list_N_matrices.append(U_i)
        
#        print(U_list_N_matrices)
        
    U_full = reduce(np.matmul, U_list_N_matrices) # here adjoint the full time evolution 
        # here U_list_N_matrices = [U_0, U_1, ..., U_6]
    return U_full 

# Here the order of U_i multiplied is wrong, the actual time - ordered U_full should be -
# - DESCENDING ; we will come back later
# ASCENDING = DESCENDING, SO DOESN'T MATTER 

#%%
#test the order of U_i appended:

#test_order = []
    
#for m in range (6): # also, for (loop) index in range N, it's always 0,1,...,(N-1)
#    M_m = m
#    test_order. append (M_m)

#print(test_order)


#%% Here introduce the 'BASIS' VARIABLE, momentum(k vector) of x and y directions
# we have these functions generate k_y_list and k_x_list, but subsequence functions of plot phase bands/H real basis...so on
#                                                         if they don't use 'k_x_list' or 'k_y_list', then their N_x/y 
#                                                         (which is number of generated k_x and k_y isn't adjusted here)
#*******************************************************************************************
N_y_adjust = 10 #(N_y = L_y / a_0 ) # to distinguish from 'N_y’ embedded in the basis transform functions late
                # Q( also denoted at the test k basis 2Ny x2Ny Ham: is the [number of k_y list](which is N_y_adjust)
                #    here input into the [2x2 Ham], equivalent to
                #    the ['block-diagonal' form of 2Ny x 2Ny k basis Hamiltonian]? 
                
def generate_k_y_list(left_lim, right_lim, N_y):
    
    k_y_generated = np.linspace(left_lim, right_lim, N_y+1) # vary this numbers of k_y we plot here, 
                                         # will vary the desity of bulk phase plotted
    return k_y_generated[1:]

k_y_list = generate_k_y_list(0 *np.pi, 2* np.pi, N_y_adjust)
#k_y_test_list_no_including_2pi = k_y_list[ :-1]     # trail for separate phase bands plot
                          
#print(k_y_list, 'k_y_list')
                                        
# k_x sampling portal (not only for edge! For all kx!!)
 
#*******************************************************************************************

N_x_adjust = 50 #(N_x = L_x / a_0 )
def generate_k_x_list(left_lim, right_lim, N_x):
    
    k_x_generated = np.linspace(left_lim, right_lim, N_x+1) # vary this numbers of k_y we plot here, 
                                         # will vary the desity of bulk phase plotted
    return k_x_generated[1:]

k_x_list = generate_k_x_list(0 *np.pi, 2* np.pi, N_x_adjust)# very important - here we defined the independent var k_x! 
#print(generate_k_x_list(0 *np.pi, 2* np.pi, 3))
#print(k_x_list, 'k_x_list')

#If want to vary list k_x(ex: how many k_x we plot against within 0-2pi), please vary it here 
# (Here we use the same k_x list for bulk plot as well as edge plot
#**********************************************************************************

#%%

#plt.figure()

# Here this function below is designed for how to plot the quasienergy (epsilon_n, which is arg of the eigval of U_full's )-
#- against kx, for specific k_y value 


def plot_phase_vs_kx(k_x, k_y):
    Ufull_matrix_asc_list =[]
    eigval_U_asc_T_list = []
#    eigvec_U_asc_T_list = []
    phase1_T_asc_list = []
    phase2_T_asc_list = []
    quasi1_T_asc_list = []
    quasi2_T_asc_list = []
    
    
    for i in range (len(k_x)):
        
        Ufull_matrix_asc_list.append(U_full_ascending(k_x[i], k_y, 6))
    #    print(i)
#        print(Ufull_matrix_asc_list)
        
        
    for i in range (len(k_x)):
        
        eigvals, eigvecs = np.linalg.eig(Ufull_matrix_asc_list[i])
    #    print(eigvals,'--eigen value;', eigvecs,'--eigen vector', i,'i') # all these eigen values are complex!!

        eigval_U_asc_T_list.append(eigvals)
        phase1_T_asc_list.append(- np.angle(eigvals[0])) #these are epsilon_n*T = Phi_n (n here =1, 2) -- which is what we plotted
        phase2_T_asc_list.append(- np.angle(eigvals[1])) # see above; but anyway we setted T = 1 before, so-- 
        quasi1_T_asc_list.append(- np.angle(eigvals[0]/T))  # quasienergy epsilon_n = phase_n/T (MAYBE WE WILL USE THIS LATER)                                              
        quasi2_T_asc_list.append(- np.angle(eigvals[1]/T))  # inrigorously, phi_n = epsilon_n, phase band is quasiE here
    # here phase1 and phase2 are for 2 (arg of) eigenvals of 2by2 H matrix
    
    #OBSERVATION:
    # ***We can see that: both of phase1 and phase2 (for all kx (and ky) under list), -
    # - the phase pattern is central symmetry against 0 
    

# -- Does here to use a list A filled in by last i loop, it needed to have another new i loop to distribute A[i]?




from matplotlib.lines import Line2D # This is purely for plot legend 

legend_elements = [
    Line2D([0], [0], marker='.', color='tab:blue',
           linestyle='None', markersize=6,
           label='Phase band 1 (1st eigenvalue)'),
    Line2D([0], [0], marker='.', color='tab:purple',
           linestyle='None', markersize=6,
           label='Phase band 2 (2nd eigenvalue)')
]

#%% try U_full_discretisation
def plot_phase_vs_kx(k_x, k_y):
    Ufull_list =[]
    eigval_U_asc_T_list = []
#    eigvec_U_asc_T_list = []
    phase1_T_asc_list = []
    phase2_T_asc_list = []
    quasi1_T_asc_list = []
    quasi2_T_asc_list = []
    
    
    for i in range (len(k_x)):
        
        Ufull_list.append(U_full_discretisation(k_x[i], k_y, 100))
    #    print(i)
#        print(Ufull_matrix_asc_list)
        
        
    for i in range (len(k_x)):
        
        eigvals, eigvecs = np.linalg.eig(Ufull_list[i])
    #    print(eigvals,'--eigen value;', eigvecs,'--eigen vector', i,'i') # all these eigen values are complex!!

        eigval_U_asc_T_list.append(eigvals)
        phase1_T_asc_list.append(- np.angle(eigvals[0])) #these are epsilon_n*T = Phi_n (n here =1, 2) -- which is what we plotted
        phase2_T_asc_list.append(- np.angle(eigvals[1])) # see above; but anyway we setted T = 1 before, so-- 


    

#%%
# The loop(PLOT) below is: 
#try to not plot 2pi, for that I feel the phase1/phase2 separate phase band plot seems having a bit of strange shift -
# - which the same type of that when I plot only specific k_y value(also, for k_y doesn't equal to 0, 
# - phase1 or phase2's plot doesn't need to obey particle-hole symmetry)

# THIS IS JUST A TEST WHICH DOES NOTHING
#for i in range (len(k_y_test_list_no_including_2pi)):
    
#    plot_phase_vs_kx(k_x_list, k_y_test_list_no_including_2pi[i])
    


#%%

#-----REAL BASIS HAMILTONIAN MATRIX H 2NX2N-----
# This below cell shouldn't be put just after 2by2 lattice coupling Hamiltonian
#Please don't confuse the 2by2 sublattices coupling Hamiltonian with the (2N by 2N) REAL SPACE Hamiltonian!!
#(here this  H2N matrix is for specific k_x)
# here y and y' denoted the hopping relation between y and y'. 
#-----------------basically this is H_y y't----------------------

# here we used inverse FT in y direction, involved (1/N_y) (--not L_y here)
def real_H2Ny_PBC_fix_k_x(N_y, k_x, t): # t can be any from [  t_i = T/N * i + T/(N * 2)  ] -- Here N is Ly
                                 #NOTICE: for the subsequence codes, Ly(number of sites)
    
    M = np.zeros((2*N_y, 2*N_y), dtype =complex)
    
#    a_0 = 1 # here is the lattices length
    
    # here we need to sum agaisnt k_y, for corresponding (2pi/Ny)* z; we used integer z instead, so that we can alsways adjust the length of ky_list
    # Which is different from the beginning 
    for z in range (1, N_y+1):
        
        k_y = 2* (np.pi/ (N_y*a_0)) * z  # k_y = 2pi/N -- This **N** has went through the whole process of computing OBC 
                                 # N - real basis converted H_yy'kx (only Nx/y = number of kx/y, density of phase bands are same)
                                 # N - turn off corner
                                 # N - U 
                                 # all these N denote the number of unit cells (num of y or x)
        H_ky = H_matrix_2by2(k_x, k_y, t, a_0) # here is how we introduced 2by2 Hamiltonian 
        
        
        
        for i in range(N_y): #--i in range(N) is for row, j in range N is for column. So y-y' here is y for row y'for column 
            for j in range(N_y):
                M[2*i:2*i+2, 2*j:2*j+2] += (1/N_y) * (np.exp(1j*  k_y * ((i+1) - (j+1)))) * H_ky
    # M is covered by ky rounds and rounds, that's how we sum against  ky
    return M



#create real object to test this real 2N by 2N matrix:
    
#H_specific_t_stage = real_basis_H2N_fix_k_x(4, k_x_list[1], (T/12)*11)

#(here N =10, time stage is at the 6th(midpoint of stage 6)


#%% HERE WE WRITE A TEST FUNCTION FOR PLOT PHASE BAND (equasienergy spectrum) OF real basis PBC in xy H
# list of 6 Hamiltonian, Periodic boundary condition(PBC) in y (fixed PBC kx)
def H_6_2Ny_PBC(N_y, k_x): #here this function should have same spectrum(EIGEN_yVALS)as the [K BASIS H2x2]
                    #6 refers as list of 6 piecewise H 
    
    H2Ny_PBC_6_matrices = []
    
    for x in range (6): #x refers index for t
    
        t_x = T/6 *  x + T/(6 * 2)   # here this x is equivalent to i in the 2x2 H list func
        
        real_H2Ny =  real_H2Ny_PBC_fix_k_x(N_y, k_x, t_x)
        
        H2Ny_PBC_6_matrices.append(real_H2Ny) # generate 6 2N_y real H
        
    return H2Ny_PBC_6_matrices
        
def Ufull_descending_PBC(N_y, k_x): # this is what it should be: U6U5U4U3U2U1 inverse from the bbnotes
# var index: 1        
    U_list_1 = []
#reference :
#    H_list_6_2N_y_edge = H_6_2N_y_OBC(N, k_x)
    
    for i in range (6):

        U_i = expm(-1j* (H_6_2Ny_PBC(N_y, k_x)[i])* (T/6))
        
        U_list_1.append(U_i) # append should be in the loop scope, 
                            # one grid inside the for loop title if you want to append the thing produced once loop (per i)
    
    U_reverse_1 = np.flip(U_list_1)
    
    U_full_1 = reduce(np.matmul,  U_reverse_1)
        
    return U_full_1

def plot_check_phase_vs_kx(k_x, N): # here in this function you can choose N_y (= Ly/a_0) to be any number
                                    # if choose same number as len(k_y_list) --> num of y here = initial N_y generate ky list 
                                    # this replot (plot_check) of phase band will be same as the initial plotfunc [plot_phase_vs_kx]
                            # *here k_x is a list
                            # plot to check the [phase band] of H_y_y'_kx PBC same as H2x2
    Ufull_descending_PBC_kx_list =[]
    phi_n_list_of_kx_list = []
    kx_list_for_plot = [] 
    for i in range (len(k_x)): # This loop's each round run one specific k_x
        
        Ufull_descending_PBC_kx_list.append(Ufull_descending_PBC(N, k_x[i]))

        eigvals, eigvecs = np.linalg.eig(Ufull_descending_PBC_kx_list[i]) # will this really return anything? 
        # YES! YOU CAN EXTRACT THE NEW APPENDED LIST IN SAME SCOPE
        
        eigvals_fixed_kx_list = eigvals # just rename
#        print(len(eigvals_fixed_kx_list), eigvals_fixed_kx_list,'eigvals of fixed k_x of Ufull')
 #       print(eigvals_kx_2Ny_list, 'eigevals')
        
        phi_n_fixed_kx = [] # phi as phase bands labelled by n
        
        for j in range (len(eigvals_fixed_kx_list)): # j labels phase n, epsilon_n(double check at home)
            
            phase_kx_j = - np.angle(eigvals_fixed_kx_list[j])
            
            phi_n_fixed_kx. append (phase_kx_j) # this is list of phase for specific k_x, length 2N
            
#        print(phase_kx_2Ny_list, 'phases')
       
        phi_n_list_of_kx_list. append(phi_n_fixed_kx) # this will be a list of different k_x, for each k_x has j
                                                     # -- layer one: len(k_x) -- layer two(each k_x), j_phase
        
#       print [k_x[i] * len(phi_n_fixed_kx)]
        kx_list_for_plot.append([k_x[i] * len(phi_n_fixed_kx)])
#        print(len(phi_n_fixed_kx))
    

#%%    
#    print (len(kx_list_for_plot))
#    print(len(phi_n_list_of_kx_list))

#CHIRAL EDGE STATE PLOT
Ny_check = 10 # this is the Ly = N*a0 before Ly defined 
              # if choose same number as len(k_y_list) --> num of y here = initial N_y generate ky list 
              # this replot (plot_check) of phase band will be same as the initial plotfunc [plot_phase_vs_kx]

test_phi_along_kx = plot_check_phase_vs_kx(k_x_list, Ny_check)

#%%

#Followed by: def real_basis_H2N_fix_k_x(N, k_x, t): # t can be any from [  t_i = T/N * i + T/(N * 2)  ]

# this below function will produce a list of 6 H
#----------------------------------------------------

# HERE IS A LIST OF 6 PIECEWISE H FOR Y OBC(H_y, y', kx), OBC

# here we used inverse FT in y direction, involved (1/N_y) (--not L_y here)
def H_6_2Ny_OBC(N_y, k_x): # here num of stages are fixed to be 6; lattice length a_0 fixed to be 1
                           # here k_x is not a list but [a number], a specific kx
    H2Ny_list_6_matrices =[]
    
    H2Ny_list_edge_6 = []
    
    for x in range (6):
    
        t_x = T/6 *  x + T/(6 * 2)   # here this x is equivalent to i in the 2x2 H list func
        
        real_H2Ny =  real_H2Ny_PBC_fix_k_x(N_y, k_x, t_x)
        
        H2Ny_list_6_matrices.append(real_H2Ny) # generate 6 2N real H
        
        # from here, we are trying to turn the 2N matrix to edge condition
        real_H2Ny_edge = real_H2Ny. copy()
        rows, cols = real_H2Ny_edge .shape
    
        top_right = [(i, j ) for i in range(rows) for j in range(cols) if i < 2 and j >= cols - 2]
        top_right_indices = top_right[:4]
        
 

# Bottom-left 4 entries
        bottom_left = [(i, j) for i in range(rows) for j in range(cols) if i >= rows - 2 and j < 2]
        bottom_left_indices = bottom_left[:4]

# Turn off corner entries (we have set each of the Hamiltonians for each stage has corner terms 0)
        for i, j in top_right_indices + bottom_left_indices:
            real_H2Ny_edge[i, j] = 0
            
        H2Ny_list_edge_6. append(real_H2Ny_edge)
    
    return H2Ny_list_edge_6 # this will return the list of 6 edge functions which turned off corner already 


test_func_H_6_2Ny_OBC = H_6_2Ny_OBC(10, k_x_list[3])


#%%
#continue by 2N by 2N for edges, the full time (expm) matrix for 2N by 2N dimension:
    
#k_x_list = np.linspace(0,  2 * np.pi, 100)

##IMPORTANT: Here this function provide U_full, which is the core of 2N x 2N matrix of y y' (real basis); from H_real_edge( no corner term)

    
def U_full_ascending_edge(N_y, k_x): # which means along lhs to rhs, from (1) earlier/smaller t to (6) later/larger t
    
# here still, the num of stages are fixed as 6

    U_list_6_matrices =[]
    
    H_list_6_2Ny_edge = H_6_2Ny_OBC(N_y, k_x)
    
    for i in range (6):
    
   #     t_i = (T/N) * i + T/(N * 2)
    
        H_i = H_list_6_2Ny_edge[i]
#expm(-1j* H_i_test* delta_t_test)
        
        U_i = expm(-1j* H_i* (T/6)) #!!!!!
        
        U_list_6_matrices.append(U_i)
        
#        print(U_list_N_matrices)
        
    U_full = reduce(np.matmul, U_list_6_matrices)
        
    return U_full

#%%
#***************************************

 # adjust portal, this will be input into function through 'N'
#Ly = 10 # here is the strip in real space y direction of our model -- kind of useless in this section
 
Ny_plot_yOBC_edge = 10 # which is for chiral edge state plot, but not for edge density diagnosis --  


#%%
#*********************************************

# FOLLOWED BY [real_basis_H2N_fix_k_x] FUNCTION, or saying: H_y_y'_kx for y OBC, 
# here in the plot function including the full time evolution operator of it [U_full_ascending_edge] 

def plot_edge_vs_kx(k_x, N_y): # here k_x is a list(taken from  # here the plot function 
    Ufull_kx_list =[]
    phase_list_of_list = []
#    phase2_T_asc_list = []
    
    kx_list_for_plot = [] 
    for i in range (len(k_x)):
    
        
        # This loop's each round run one specific k_x
        
        Ufull_kx_list.append(U_full_ascending_edge(N_y, k_x[i]))
    #    print(i)
#        print(Ufull_matrix_asc_list)
        
        
    for i in range (len(k_x)): # this is just for using the U_full as a list, it is actually still in 
                               # -- specify only k_x layer
                               
         # for each fixed k_x, we have one Ufull_kx_list[i], corresponding to 2N eigvals 
        
        eigvals_edge, eigvecs_edge = np.linalg.eig(Ufull_kx_list[i]) # here should be 2N eigen vals, here this eig fun
                                                           # -- will produce eigvals as a list 
    #    print(eigvals,'--eigen value;', eigvecs,'--eigen vector', i,'i') # all these eigen values are complex!!
        eigvals_kx_2Ny_list = eigvals_edge  # just rename
        
 #       print(eigvals_kx_2Ny_list, 'eigevals')
        
        phase_kx_2Ny_list = []
        
        for j in range (len(eigvals_kx_2Ny_list)):
            
            phase_kx_j = - np.angle(eigvals_kx_2Ny_list[j])
            
            phase_kx_2Ny_list. append (phase_kx_j) # this is list of phase for specific k_x, length 2N
            
#        print(phase_kx_2Ny_list, 'phases')
       
        phase_list_of_list. append(phase_kx_2Ny_list) # this will be a list of different k_x, for each k_x has j
                                                     # -- layer one: len(k_x) -- layer two(each k_x), j_phase
        
    
        
        kx_list_for_plot. append([k_x[i]] * len(phase_kx_2Ny_list))
#          for kx_val in k_x:
#        U = U_full_ascending_edge(N, kx_val)
#        eigvals, _ = np.linalg.eig(U)
#        phases = -np.angle(eigvals)  # shape: (2N,)
        
#        phase_list.extend(phases)
   #     kx_expanded.extend([kx_val] * len(phases))  # repeat kx_val for each phase

#    plt.figure(figsize=(8, 5))
#    print (len(kx_list_for_plot))
#    print(len( phase_list_of_list))

            
  #  plt.plot(k_x, phase_list_of_list,  'o', label='edge phase')
            

        
 #   for y in range (len(eigval_kx_2N_list)):
        
        
        
#        phase_kx_2Ny_list. append (phase_kx_2N)
        
#        phase1_T_asc_list.append(- np.angle(eigvals[0]))
#        phase2_T_asc_list.append(- np.angle(eigvals[1]))

# -- Does here to use a list A filled in by last i loop, it needed to have another new i loop to distribute A[i]?

#    plt.plot(k_x, phase1_T_asc_list,  'o', label='phase 1(1st eigval), phi_n = phi_0')
#    plt.plot(k_x, phase2_T_asc_list,  'o', label='phase 2(2nd eigval), phi_n = phi_1')

    
#for i in range (len(k_y_list)):
    

#    plot_phase_vs_kx(k_x_list, k_y_list[i])

#CHIRAL EDGE STATE PLOT
    


#%%
#now try eigen vectors --

#-----------------------------------------------------------------------------
# IMPORTANT-- FROM HERE, WE ARE BEGINNING PLOTTING THE INTENSITY OF STATES(EPSILON_N)

#kx_list_edge has 10 arguments

def find_Phi_n_vec_vs_kx(k_x, N_y): # here k_x is a list  N is number of y  = 1...N_Y
 #   Ufull_kx_list =[]
    eigen_vec_listed_by_kx = []
#    phase2_T_asc_list = []
    
#    kx_list_for_plot = [] 
    for i in range (len(k_x)):
    
        
        # This loop's each round run one specific k_x
        
#        Ufull_kx_list.append(U_full_ascending_edge(N, k_x[i]))
    #    print(i)
#        print(Ufull_matrix_asc_list)
        
        
#   for i in range (len(k_x)): # this is just for using the U_full as a list, it is actually still in 
                               # -- specify only k_x layer
                               
         # for each fixed k_x, we have one Ufull_kx_list[i], corresponding to 2N eigvals 
        
        eigvals_edge, eigvecs_edge = np.linalg.eig(U_full_ascending_edge(N_y, k_x[i])) # here should be 2N eigen vals, here this eig fun
                                                           # -- will produce eigvals as a list 
                                                           #this is for the edged -- obc, not pbc in y 
    #    print(eigvals,'--eigen value;', eigvecs,'--eigen vector', i,'i') # all these eigen values are complex!!
        eigen_vec_listed_by_kx.append(eigvecs_edge)
        
        print(i, 'eigen vectors shape:', eigvecs_edge.shape)
        print(i, 'eigen vectors shape:', U_full_ascending_edge(N_y, k_x[i]).shape)
        # here we checked by printing: U_full_ascending_edge(N, k_x[i]) is 20x20 matrix -- why? wait, yes it should be
        # U_full_ascending_edge(N, k_x[i]) should be 2N by 2N matrix. Here we have N = Ly = 10
        #Eigenvectors matrix: (N x N, here N =3)
        #[[1. 0. 0.]
        #[0. 1. 0.]
        #[0. 0. 1.]]
        #First eigenvector: [1. 0. 0.]
        #so for our 20x20 matrix, each of the Phi_n is a eigen vector with 20 component. 

    return eigen_vec_listed_by_kx
    
# here this function print out eigen vecs for each kx has 20x20, why? 

#print (k_x_list)
        
Phi_n_vec_list = find_Phi_n_vec_vs_kx(k_x_list, Ny_plot_yOBC_edge ) # there are 20 eigvecs for 2N x 2N (20) matrix 

                                    #calculate the sum of amp of them for a specific pair: first pair: y =1, alpha =1 or 2

# how could phi_list be a 20x20 matrix?? -- this should be 10 x 20, for 10 kx??
#find the amp :
#Rho_y = []
#for i in range (10):
#    Phi_y = (Phi_list[i])**2 +   (Phi_list[i+1])**2     
#    Rho_y.append(Phi_y)
#%%
# here this function will find all the eigvalues

def find_epsilon_n_vs_kx(k_x, N_y): # here k_x is a list  N is number of y I guess, 1...LY=20 if N =20 (see previously
                                                                   # Ly = ?
 #   Ufull_kx_list =[]
    eigvals_listed_by_kx = []
    phase_ns_list_listed_by_kx =[] # this will contain the arg of U_full matrix
#    phase2_T_asc_list = []
    
#    kx_list_for_plot = [] 
    for i in range (len(k_x)):
        
        eigvals_edge, eigvecs_edge = np.linalg.eig(U_full_ascending_edge(N_y, k_x[i])) # here should be 2N eigen vals, here this eig fun
                                                           # -- will produce eigvals as a list 
                                                           #this is for the edged -- obc, not pbc in y 
    #    print(eigvals,'--eigen value;', eigvecs,'--eigen vector', i,'i') # all these eigen values are complex!!
        eigvals_listed_by_kx.append(eigvals_edge)
        # here we gain a list of list if stay outside, in this [inside] level, we are still have only 
        # list len 20 for n eigvals, specific kx, but can we use it all for angle(phase)？ 
        
 #   for n in range (N):
        
        phase_n = -np.angle(eigvals_edge ) #- can we apply np.angle on a list??
        
        phase_ns_list_listed_by_kx.append(phase_n) # here this phase is a list, because for 2N x2N matrix(N =20)
                                                   # we firstly have 2N eigenvecs corresponding to 2N eiegnvals - 
                                                   #                                          1 eigenval - 1 phase
                                                   
                                                   
                                                   
                                                   # Q: does each eigvec ~ 1 specific eigstate here (I guess so: coz -
                                                   # - the eigen state here is for full revolution, tho each Ui has -
                                                   # - their own eigenVEC, but they are for the same state -, ***?
                                                   # -just under time evolution(so can I understand as that time evolution
                                                   # is for basis? not for eigen state?? But eigenval also changed)
                                                   # But one thing is for sure: using momentum basis-
                                                   #- ( H in momentum space) or using real basis (H 2N by 2N) are exactly
                                                   # corresponding to this 'same eigenstate' diff basis - diff 'eigVEC'
                                                   
                                                   
                                                   
                                                   # hence it's a list of phase [for fixed kx] len(2N)
                                                   # and all lists of phase(len(2N)) are listed by kx again
                                                   
                                                   
                                    # ---- JUST TO CONFIRM: WE DO NEED TO FIND EIGENSTATE BY EIGENVEC OF U_full------
                                                   
                                                   
    
#        print(i, 'eigen vectors shape:', eigvecs_edge.shape)
#        print(i, 'eigen vectors shape:', U_full_ascending_edge(N, k_x[i]).shape)
        # here we checked by printing: U_full_ascending_edge(N, k_x[i]) is 20x20 matrix -- why? wait, yes it should be
        # U_full_ascending_edge(N, k_x[i]) should be 2N by 2N matrix. Here we have N = Ly = 10


    return  phase_ns_list_listed_by_kx

# we input kx in this 'find_epsilon_n_vs_kx’ function, note that kx is a list!!
epsilon_n_list = find_epsilon_n_vs_kx(k_x_list, Ny_plot_yOBC_edge) # here 10 is for N, which is Ly, which is 10 for 2N x2N  = 20x 20 matrix 


#%% print the list of eigvecs for 2N by 2N U_full, find the density amp of each y


# this return the first and second kx (kx[0] and kx[1])'s 20 eigen vec!!
#(roughly checked out)

#%% checking about the edge state density

#**************** BEGINNING OF LOOPS FOR FINDING RHO_EDGE_N, EPSILON_N AND INT PLOTTING ***************

Rho_all_y_listed_by_kx = []

Rho_edge_listed_by_sta_n_listed_by_kx = []

for i in range(len(k_x_list)):  # Here number of i range corresponding to k_x_list
    # take absolute square of eigenvector amplitudes
    Phi_specific_kx =  Phi_n_vec_list[i] # 20 x 20 matrix for eigvecs, corresponding tos specific kx
    
    # here Phi_specific_kx has 20 eigen vecs
    # we only consider the functions within one specific kx now. (within one round labelled by i)
    
    Rho_y_listed_by_y = []
    
    
    for j in range(Ny_plot_yOBC_edge): #j labels y (or range N)
        # here for specific alpha(sublattices)
        
        Phi_vec_n1 = Phi_specific_kx [:, 2*j]  # a vec(array(3, )) with 20 components
        Phi_vec_n2 = Phi_specific_kx [:, 2*j + 1] #the next vec(array(3, )) with 20 components
        
        comp_wise_modulus_sqr_1 = np.abs(Phi_vec_n1)
      
        comp_sum_for_1 = np.sum(comp_wise_modulus_sqr_1)
       
        
        Rho_specific_y = np.sum(np.abs(Phi_vec_n1)) + np.sum(np.abs(Phi_vec_n2))  # Q: what does this actually do? 
            # this corresponding to modulus square for each Phi_n pair in sublattices, which is Rho_n
  
        Rho_y_listed_by_y. append(Rho_specific_y)
        
#        for k in range(len(vec_y_alpha1)): #(actually this is just j)
#actually the above is useless -- let's try this n label below, for each n corresponding to one Phi_n, state_n
    
    Phi_eigvecs_n_list = []
    
    
    Phi_eigvecs_n_Rho_edge_list = []
    
    Phi_eigvecs_n_list_Rho_comp_y_list =[]
    
    for n in range(2*Ny_plot_yOBC_edge): # I guess here is to spam both vec_n1 and vec_n2? ******Why here is 40???
        
        # {{in this loop}}, we have the Phi_list(comp of eigenvec list) for each eigenvec n; corresponding to epsilon_n
        Phi_eigvecs_n = Phi_specific_kx [:, n]
        
        Phi_eigvecs_n_list. append(Phi_eigvecs_n)
        
        
        #here try to find component y', alpha1 and y' alpha2 (y and y' are the same thing)
        Rho_comp_Phi_eigvec_specific_n_list = []
        
#        N = 20 #(only for record: here we have Ly from 1 to 10, and this number of real y basis is N(showed previously))
        
#        Ly = N
        
        for y in range (Ny_plot_yOBC_edge): # {{in this loop}}, we collect component of each eigenvec n (last loop)--
                                            # -- these components are labelled by: 1_alpha1/2, 2_alpha1/2,3,..., Ny_alpha1/2
            Rho_component_y = (np. abs(Phi_eigvecs_n[2*y]))**2 + (np.abs(Phi_eigvecs_n[2*y + 1]))**2 #Define Rho
            
                       
            # here the Rho correponding to each eigvec n' component's modulus sqaure. But np.abs gives the modulus
            
            Rho_comp_Phi_eigvec_specific_n_list.append(Rho_component_y)
            
          # for check 1
              # or go to the var explorer, see Phi_eigvecs_n_list_Rho_each_y_list, this list is for the dos for each y and - next line 
              #each n, should be len(n[len(Ly)])
              
 #       try change the def of 'edge'
            
        Rho_edge_Phi_eigvec_n = Rho_comp_Phi_eigvec_specific_n_list[0]#+ Rho_comp_Phi_eigvec_specific_n_list[Ny_plot_yOBC_edge -1] #+ Rho_comp_Phi_eigvec_specific_n_list[2] + Rho_comp_Phi_eigvec_specific_n_list[Ly -3] + Rho_comp_Phi_eigvec_specific_n_list[Ly -2] + Rho_comp_Phi_eigvec_specific_n_list[Ly -1]
        #total density of edege state #  From those components(1,...y_some, Ny) of each eigenvec_n(also 2Ny num of),
                                      #  identify that of edges; or that of [any position's] -- we adjust here
        #total density of edege state
    # for check 2 
        
        Phi_eigvecs_n_Rho_edge_list.append(Rho_edge_Phi_eigvec_n ) # ****here append all these position identified Rho to all eigenvec n
        
        Phi_eigvecs_n_list_Rho_comp_y_list.append(Rho_comp_Phi_eigvec_specific_n_list) # When finish appending (out of this loop of vec_n)
                                                                                       # This will return a list of all Rho_y=1  -|
                                                                                       #                                 Rho_y=2 -| --Phi_n = 1|
                                                                              #                                             (eigenvec n leballed as for all quasienergies)
    # here I want to check two things: 1. for each Phi n, among desities labelled by y, is there any edge showed
    #                                     higher probability?
    #                                  Conclusion: Not really??? So what does it mean by localised? (Maybe can ask professor)
    
    #                                  2. doubel check if the Rho_edge is really the sum of the two of Rho_comp, 0 and 10(as index 9)?
    #                                  Conclusion: checks out 
    
            
        
    Phi_specific_kx_20_eigvecs = []
    
    for k in range(2*Ny_plot_yOBC_edge): # here to check 20 eigvecs
        
        Phi_specific_kx_20_eigvecs.append(Phi_specific_kx[k])
        

                                   
   # print('i =', i, Rho_y_listed_by_y, 'rho list for all 10 y, should be len = 10') 
    
    Rho_edge_specific_kx = Rho_y_listed_by_y[0] + Rho_y_listed_by_y[Ny_plot_yOBC_edge -1] #+ Rho_y_listed_by_y[Ly -2] + Rho_y_listed_by_y[Ly -3] + Rho_y_listed_by_y[Ly -4] + Rho_y_listed_by_y[Ly -10]# here 19 is Ly -1
    # Q: why increase this can actually increas e the number of Rho_edge_listed_by_sta_n_listed_by_kx
#    for y_basis in range (Ny_plot_yOBC_edge):
#        Rho_specific_y_specific_kx = Rho_y_listed_by_y[0]
#        Rho_sum = + Rho_specific_y_specific_kx
#    print (Rho_sum) # this is labelled by specific kx and specific 
#    Rho_edge_listed_by_kx. append(Rho_edge_specific_kx)
    
    Rho_all_y_listed_by_kx.append(Rho_y_listed_by_y)
    #Rho_y_listed_by_kx.append(Rho_y_from_1_to_Ly_list)
    # here this Rho_y_listed_by_kx is for all kx, packed Rho_y, should be a list of 10 list(for each of 10 kx)
    # and for each of 10 kx, there is 10 Rho for 1....Ly = 10
    
#print( Rho_all_y_listed_by_kx, 'density for specific kx') # this is the eigen vec(with alpha 1, alpha2) 's correponding density 
# this is all state density, not only edge state density! 
#print(Rho_edge_listed_by_kx, 'edge density for specific kx')

# Rho_edge_listed_by_kx, tried kx = kx_list[1], kx = kx_list[2], the edge Rho(desity) = Rho_y=1 + Rho_y = Ly
# select the edge density: Rho_y_y=1 and Pho_y_y = Ly
# Phi_eigvecs_n_Rho_edge_list now appended into kx 

    Rho_edge_listed_by_sta_n_listed_by_kx.append(Phi_eigvecs_n_Rho_edge_list)
    

#%% here in our plot of phase band structure with intensity (labelled by the density of edge)

# which need to be corrected 
 

kx_temperoally_in_p_loop =[]
rho_edge_ns_specific_kx = [] 
epsilon_ns_specific_kx = []


# define the norm of c in color map to make it unified for each kx!

#all_rho = np.concatenate(Rho_edge_listed_by_sta_n_listed_by_kx)
#vmin = all_rho.min()
#vmax = all_rho.max()
#norm = PowerNorm(gamma=0.4, vmin=0, vmax=1)
 

for p in range(len(k_x_list)):  # Here number of i range corresponding to k_x_list
    
    kx_in_p_loop = k_x_list[p] # here we use p as the index for kx
    
    rho_edge_ns_specific_kx = Rho_edge_listed_by_sta_n_listed_by_kx[p] 
    
    Phi_ns_specific_kx =  Phi_n_vec_list[p] # useless??
    epsilon_ns_specific_kx = epsilon_n_list[p]
    
   # append(kx_temperoally_in_p_loop)
    
  #  kx_temperoally_in_p_loop =[]
    
#    for state in range(len(rho_edge_ns_specific_kx)): # state is n in the example
        
#        epsilon_n = epsilon_ns_specific_kx[state]
        
#        rho_edge_n = rho_edge_ns_specific_kx[state]
        
#        print (epsilon_n, 'epsilon_n;', rho_edge_n, 'rho_edge_n, for state', state)
        


#%% real basis x, y PBC, check spectrum
# Plan A: [H(2x2)] to x,y real basis [PBC(2NyNx 2NyNx)]; check spectrum: compare with Nx&Ny generated kx&ky list as vars-all epsilon

# Plan B: [PBC in y H(2Ny x 2Ny)] to x,y real basis [PBC(2NyNx x 2NyNx)]; 
#         check spectrum: compare with Nx generated kx list, same Ny 
#                         - this 'check Plan B spectrum' already compared with 'check Plan A spectrum'


# spectrum need to be find by H(2NxNy x 2NxNy) append 6 piecewisely(will be changed to time discretisation calculate U(T) - epsilon
# AIM: to make the specrtum as the same

#   
#%%
def H2NyNx_PBC_version2(N_y, N_x, t, a_0):
    
    """H_ele = np.zeros((N_y, N_x, N_y, N_y, N_x, N_x, 2, 2), dtype=complex) # here I just want to write a STORAGE: 
                                                   # which is that we pull out values from this matrix later
    
   # k_storage = np.zeros((N_y, N_x, 2, ))
    for iy in range(N_y): #N_x is the dim, ix = 0,,,N_x-1 if call N_x as loop range
        k_y_inside = 2*np.pi*(iy+1)/N_y
        for ix in range(N_x): # HERE fixed k_y& iy , varying ix
            k_x_inside = 2*np.pi*(ix+1)/N_x 
            Hk_22 = H_matrix_2by2(k_x_inside, k_y_inside, t, a_0)
            
            for y_index in range( N_y):
                for y_1_index in range( N_y):
                    for x_index in range( N_x):
                        for x_1_index in range( N_x):      
                            
                            H_ele[iy, ix, y_index, y_1_index, x_index, x_1_index] = (1/N_x) * (1/N_y) * np.exp(1j*  k_y_inside * (y_index - y_1_index ))*np.exp(1j*  k_x_inside * ( x_index - x_1_index)) * Hk_22
#            k[iy, ix] = np.array[k_y_inside, k_x_inside] # this is an (2,) array record ky and kx 
    H_ele_sum_over_k = np.sum (H_ele, axis = (0,1))     
    """

    H_dif_sum = np.zeros((N_y, N_x, 2, 2), dtype = complex)
    for iy in range(N_y): #N_x is the dim, ix = 0,,,N_x-1 if call N_x as loop range
        k_y_inside = 2*np.pi*(iy+1)/N_y
        for ix in range(N_x): # HERE fixed k_y& iy , varying ix
            k_x_inside = 2*np.pi*(ix+1)/N_x 
            Hk_22 = H_matrix_2by2(k_x_inside, k_y_inside, t, a_0)
            
            #coest NxNy
            
            for y_dif in range( N_y):
                for x_dif in range( N_x):
                    H_dif_sum[y_dif, x_dif] += (1/N_x) * (1/N_y) * np.exp(1j*  k_y_inside * y_dif)*np.exp(1j*  k_x_inside * x_dif) * Hk_22
                    #NxNy
                    
# changing based on guess: call storage of results of one function computed 100 times is much faster than
# compute this function literaly 100 times 
    H = np.zeros((2*N_y*N_x, 2*N_y*N_x), dtype =complex)
    for i in range(N_x* N_y):
        for j in range(N_x* N_y):
            x = (i )% N_x 
            x_1 = (j )% N_x 
            delta_x = (x - x_1) % N_x #this gives [the smallest integer which near the (x-x1) and can be divided by N_x]'s difference with (x-x_1)
            y = (i) // N_x            #(x-x_1) -[smalles dividable integer of N-X]
            y_1 = (j) // N_x  # y'
            delta_y = (y - y_1) % N_y
            #here this loops i and j are fine since there isn't sum over them
            
            # 4 is indexed
            H[2*i:2*i+2, 2*j:2*j+2] = H_dif_sum[delta_y, delta_x]   #H_ele_sum_over_k [y, y_1, x, x_1] 
    return H
    #cost 4 *   (NxNy )^2

# total cost of this function: 5 * (NxNy) ^2
    
H_test_version2 = H2NyNx_PBC_version2(10, 15, T/12, a_0)


#%%
def H_6_2NyNx_PBC_version2(N_y, N_x, a_0): #here this function should have same spectrum(EIGEN_yVALS)as the [K BASIS H2x2]
                    #6 refers as list of 6 piecewise H 
    
    H2NyNx_PBC_6_matrices = []
    
    for t_x in range (6): #x refers index for t
    
        t_piecewise = T/6 *  t_x + T/(6 * 2)   # here this x is equivalent to i in the 2x2 H list func
        
        H2NyNx_PBC_piecewise =  H2NyNx_PBC_version2(N_y, N_x, t_piecewise, a_0)
        
        H2NyNx_PBC_6_matrices.append(H2NyNx_PBC_piecewise) # generate 6 2N_y real H
        
    return H2NyNx_PBC_6_matrices
        
def Ufull_descending_PBC_realxy_version2(N_y, N_x, a_0): # this is for whole period T
# var index: 1        
    U_list_1 = []
#reference :
#    H_list_6_2N_y_edge = H_6_2N_y_OBC(N, k_x)
    H_6_2NyNx_list = H_6_2NyNx_PBC_version2(N_y, N_x, a_0)
    
    for i in range (6):

        U_i = expm(-1j* (H_6_2NyNx_list[i])* (T/6))
        
        U_list_1.append(U_i) # append should be in the loop scope, 
                            # one grid inside the for loop title if you want to append the thing produced once loop (per i)
    
    U_reverse_1 = np.flip(U_list_1)
    
    U_full_1 = reduce(np.matmul,  U_reverse_1)
        
    return U_full_1


def find_spectrum_NxNy_version2(N_y, N_x, a_0):
    eigvals, eigvecs = np.linalg.eig(Ufull_descending_PBC_realxy_version2(N_y, N_x, a_0))
    epsilon = - np.angle(eigvals)
    return  epsilon
#%% try to plot a spectrum but not only label the y = 1 or Ny position var of Rho, but for all the Rho 
# content of plotting: Rho against epsilon_n:

# firstly go back to y OBC x PBC, numerical check - including for Rho (1) + Rho(Ny), are they not inside certain epsilon(eigenvec/eigenval?)
# check aim 2: If Rho(1) +... Rho(Ny) = 1(I don't think so) -- but can plot them to see if they show something (against epsilon_n) 
#      Q: Do they show that when time vortex there's no states 
# plot rho(y = 0) or jearby sites for xy OBC's epsilon 

#%%
#--------------------------------------------------------------
#%% test code REWRITE H2NyNx_PBC


#test_H = H2NyNx_PBC(10, 10, 1/12, a_0)
#%%
#H_ky = H_matrix_2by2(k_x, k_y, t, a_0)
#M[2*i:2*i+2, 2*j:2*j+2] += (1/N_y) * (np.exp(1j*  k_y * ((i+1) - (j+1)))) * H_ky
#print(test_H,'H2NyNxtest')
#%%
# ------------After defining the 2NxNy spectrum, subsequent processes----------------
#H_NxNy (PBC) >> H_6 piecewise >> U >> epsilon spectrum
 


#%% COMPARE WITH THE ONLY Y real basis one
def find_spectrum_realy (N_y, k_x): # for calculating: spectrum_test_realy = Ufull_descending_PBC(N_y, k_x_TEST)
    Ufull_descending_PBC_kx_list =[]
    epsilon_n_list_of_kx_list = []
#    U_T = Ufull_descending_PBC(N_y, k_x)
    for i in range (len(k_x)): # This loop's each round run one specific k_x
        
        Ufull_descending_PBC_kx_list.append(Ufull_descending_PBC(N_y, k_x[i]))

        eigvals, eigvecs = np.linalg.eig(Ufull_descending_PBC_kx_list[i]) # will this really return anything? 
        # YES! YOU CAN EXTRACT THE NEW APPENDED LIST IN SAME SCOPE
        
        eigvals_fixed_kx_list = eigvals # just rename
#        print(len(eigvals_fixed_kx_list), eigvals_fixed_kx_list,'eigvals of fixed k_x of Ufull')
#        print(eigvals_kx_2Ny_list, 'eigevals')
        epsilon_n_fixed_kx = [] # phi as phase bands labelled by n
        
        for j in range (len(eigvals_fixed_kx_list)): # j labels phase n, epsilon_n(double check at home)
            
            epsilon_kx_j = - np.angle(eigvals_fixed_kx_list[j])
            
            epsilon_n_fixed_kx. append (epsilon_kx_j) # this is list of phase for specific k_x, length 2N
            
#        print(phase_kx_2Ny_list, 'phases')
       
        epsilon_n_list_of_kx_list. append(epsilon_n_fixed_kx) # this will be a list of different k_x, for each k_x has j
                                                     # -- layer one: len(k_x) -- layer two(each k_x), j_phase
        spectrum_flat = [x for sublist in epsilon_n_list_of_kx_list for x in sublist]

    return  spectrum_flat                             
                                                     
#%%
# after this we can also check with kx and ky, both k basis Ham collecting spectrums for all kx and ky but
# didn't; so if future steps having problem need to go back here to double check 

#%%
N_y_TEST = 10
N_x_TEST = 15
#k_x_TEST = [1 * 2* np.pi/4, 2  *2*  np.pi/4, 3 *2*  np.pi/4, 2*np.pi]
k_x_TEST = generate_k_x_list(0, 2*np.pi, N_x_TEST)
#spectrum_test_realxy = find_spectrum_NxNy(N_y_TEST, N_x_TEST, a_0)

#spectrum_test_realxy_rewrite = find_spectrum_NxNy_rewrite(N_y_TEST, N_x_TEST, a_0)

spectrum_test_realy = find_spectrum_realy(N_y_TEST, k_x_TEST)

#%%
spectrum_test_realxy_version2 = find_spectrum_NxNy_version2(N_y_TEST, N_x_TEST, a_0)



#%% write function of turn off x and y PBC terms, gain OBCs
# aim: double check with that when y (real basis) is OBC (turned off) but x is PBC


#%% same way write only y OBC ##########???????? -- I changed this to xOBC in the function, I don't think it's y OBC
def H2NyNx_xOBC_only(N_y, N_x, t, a_0): # here we have the x OBC only case, from xyPBC (the func essentially for efficiency)
                                        # x y both PBC function **** IMPORTANT **** : [ H2NyNx_PBC_version2] (function)
    H_PBC_generated_x_OBC_ONLY = H2NyNx_PBC_version2(N_y, N_x, t, a_0)

#    rows, cols =H_PBC_generated .shape #rows and cols are just vars - like rows = 15
    
    I, J = np.indices(H_PBC_generated_x_OBC_ONLY.shape)


    cond_i_1 = (I % (2*N_x) == 0) | ((I - 1) % (2*N_x) == 0)        # ij for top right
    cond_j_1 = ((J + 1) % (2*N_x) == 0) |  ((J + 2) % (2*N_x) == 0)
    
    cond_i_2 = ((I + 1) % (2*N_x) == 0) |  ((I + 2) % (2*N_x) == 0) #ij for bottom left
    cond_j_2 = (J % (2*N_x) == 0) | ((J - 1) % (2*N_x) == 0)

    mask = (cond_i_1 & cond_j_1 ) | (cond_i_2 & cond_j_2 )    # both i-rule and j-rule must hold
    
    H_PBC_generated_x_OBC_ONLY[mask] = 0    
       # or whatever values you want
#    print("\nMask (True means 'will be changed'):\n", mask)

    return  H_PBC_generated_x_OBC_ONLY

test_H2NyNx_xOBC_only =  H2NyNx_xOBC_only(3, 3, np.pi/4, 1)

# checked result: this should be X OBC!!! only the small (2x2 size) corners' entrices within each 'same-y' varying x blocks got turned off to 0 
#%% 
#---------*****H_2NxNy_OBC_bothx&y, CENTRE OF WORK for inducing time vortex***** ------------


#This H doesn't produce anything we can check from before
#alternative way to apply a primary check: 
#                          write H_2Nx_OBC(ky) with fixed ky and append all ky's spectrum
#                          check with the H_2NxNy_yPBC_only function above

def H2NyNx_xyOBC(N_y, N_x, t, a_0): 
    
    H_NyNx_generated_both_OBC = H2NyNx_xOBC_only(N_y, N_x, t, a_0)

#    rows, cols =H_PBC_generated .shape #rows and cols are just vars - like rows = 15
    
    I, J = np.indices(H_NyNx_generated_both_OBC.shape)


    
    big_block_top_right_i = (I >= 0) & (I < 2*N_x) 
    big_block_top_right_j = (J >= 2* N_x * (N_y -1)) & (J < 2* N_x * N_y) #Nx(Ny-1)) to j = NxNy
    big_block_bottom_left_i = (I >= 2* N_x * (N_y -1)) & (I < 2* N_x * N_y)
    big_block_bottom_left_j = (J >= 0) &(J < 2*N_x)
    
    mask = ((big_block_top_right_i) & ( big_block_top_right_j )) | (big_block_bottom_left_i) & ( big_block_bottom_left_j )   # both i-rule and j-rule must hold
    
    H_NyNx_generated_both_OBC [mask] = 0    
       # or whatever values you want
#    print("\nMask (True means 'will be changed'):\n", mask)

    return H_NyNx_generated_both_OBC

# Print test: HamNxNy at specific time t = T/12 (first piecewise) x = 3 y = 4 (should be 24x24 dim)

HamNxNy_test = H2NyNx_xyOBC(4, 3, T/12, 1) # test both xy OBC ham Ny = 4 Nx = 3

# 3.18. -- this func got tested out, as in yes, we have xy OBC here, and should start from this point 
# this cell is correct as in it is xy both OBC

#( HamNxNy_test, 'test see what both OBC hamiltonian we wrote is like')
#%% go through the way did before for H_6 U_full, no time vortex as fast run
def U2NyNx_xy_OBC(N_y, N_x, a_0): 
    def H_6_2NyNx_OBC_append(N_y, N_x, a_0): #here this function should have same spectrum(EIGEN_yVALS)as the [K BASIS H2x2]
                        #6 refers as list of 6 piecewise H 
        
        H2NyNx_OBC_6_matrices = []
        
        for t_x in range (6): #x refers index for t
        
            t_piecewise = T/6 *  t_x + T/(6 * 2)   # here this x is equivalent to i in the 2x2 H list func
            
            H2NyNx_OBC_piecewise =  H2NyNx_xyOBC(N_y, N_x, t_piecewise, a_0)
            
            H2NyNx_OBC_6_matrices.append(H2NyNx_OBC_piecewise) # generate 6 2N_y real H
            
        return H2NyNx_OBC_6_matrices

    U_list_new = []
    #reference :
    #    H_list_6_2N_y_edge = H_6_2N_y_OBC(N, k_x)
    H_6_2NyNx_new_OBC_list =  H_6_2NyNx_OBC_append(N_y, N_x, a_0)
        
    for i in range (6):

        U_i = expm(-1j* (H_6_2NyNx_new_OBC_list[i])* (T/6))
            
        U_list_new.append(U_i) # append should be in the loop scope, 
                                # one grid inside the for loop title if you want to append the thing produced once loop (per i)
        
    U_reverse_new = np.flip(U_list_new)
        
    U_full_new = reduce(np.matmul,  U_reverse_new)
            
    
    return  U_full_new
    


#%% -------- introdce time vortex : J >> H_2x2 >> write new H 2x2 in this function 
# ---------***** INTRODUCED TIME VORTEX which is original script*****----------
def H2NyNx_xyOBC_time_vortex_ori(N_y, N_x, t, a_0):  # here this t is the t got sampled, 
                                                 # within this function we will get t_phase_delay to rewrite t
    
# t --- is not the t we stuff in our HxH matrix(which calls J_tilte_a ) anymore
# where is the centra point of the vortex? Plan A: try at original( read to justify)
#                                          Plan B: at middle of all unit cells Nx, Ny plane

    
#    H_generated_both_OBC_no_vor = H2NyNx_xOBC_only(N_y, N_x, t, a_0) # this create an Ham object as input goes into func
    
    # begin to write the vortex Ham
    H_generated_both_OBC_vor = np.zeros((2*N_y*N_x, 2*N_y*N_x), dtype =complex)
    
    for i in range(N_x* N_y):
        for j in range(N_x* N_y):
            x = (i + 1 )% N_x  #double check later (?) -- x & x' are from 0, to Nx -1
            x_1 = (j + 1 )% N_x 
          #  delta_x = (x - x_1) % N_x #this gives [the smallest integer which near the (x-x1) and can be divided by N_x]'s difference with (x-x_1)
            y = ((i) // N_x) + 1 #double check later (?) -- y & y' are from 0, to Ny -1
            y_1 =( (j) // N_x ) + 1 # y'

           # print (i, j, 'test ij in time vortex matrix')

            arv_y = (y + y_1)/2 # here is a question: should we do this way from 0 for all x, x', y, y'? Ans: no. We do from 1
            arv_x = (x + x_1)/2 # 
#-----------------------------------------------------------    
            vortex_theta = np.arctan2(arv_y, arv_x)
            
    # how should we apply theta belongs to [0, 2pi) if we want 
            #theta = np.arctan2(y, x)          # (-pi, pi]
            #theta_0_2pi = np.mod(theta, 2*np.pi)  # [0, 2pi)
            # continuous (no jump) if points are ordered around the circle:
            #theta_cont = np.unwrap(theta)
#------------------------------------------------------------
            t_phase_delayed = t - vortex_theta * T/(2 * np.pi)
            
#            H_generated_both_OBC_vor[2*i:2*i+2, 2*j:2*j+2] = H_generated_both_OBC_no_vor[2*i:2*i+2, 2*j:2*j+2]
            # lhs extract 2x2 little blocks in zero martix (preparted for writing time vortexed Ham), at each loop*(NxNy)
            # rhs rewrite as H - worng: H_generated_both_OBC_no_vor is written outside the loop(labels xx' and yy' - time phase)
            # each loop should assign 2x2 extracted blocks from lhs a different time phase [H_generated_both_OBC_no_vor] with corresponding part of matrix
            # Problem: here this func[H2NyNx_xOBC_only] is a func of whole matrix, is it fine if it's large we call it again and again? 
            
            H_generated_both_OBC_vor[2*i:2*i+2, 2*j:2*j+2] = H2NyNx_xOBC_only(N_y, N_x, t_phase_delayed, a_0)[2*i:2*i+2, 2*j:2*j+2]
            
            
            
    return H_generated_both_OBC_vor

# dk what we are doing here? ####################### double check even this cell is needed later 

#%% This cell, include xy OBC condition, only changed func as [H2NyNx_xyOBC] at 3.18. - details see notes inside the cell 
def H2NyNx_xyOBC_time_vor_core(N_y, N_x, t, a_0, x_0, y_0):  # x_0 and y_0 are the position of time vor core adjustable 
    
#   to be mid point time vor: (1 + N_y)/2 , (1+ N_x)/2 

    H_generated_both_OBC_vor_core = np.zeros((2*N_y*N_x, 2*N_y*N_x), dtype =complex)
    
    
    for i in range(N_x* N_y):
        for j in range(N_x* N_y):
            x = (i + 1 )% N_x  #double check later (?) -- x & x' are from 0, to Nx -1
            x_1 = (j + 1 )% N_x 
          #  delta_x = (x - x_1) % N_x #this gives [the smallest integer which near the (x-x1) and can be divided by N_x]'s difference with (x-x_1)
            y = ((i) // N_x) + 1 #double check later (?) -- y & y' are from 0, to Ny -1
            y_1 =( (j) // N_x ) + 1 # y'

           # print (i, j, 'test ij in time vortex matrix')

            arv_y = ((y + y_1)/2) - y_0 # here is a question: should we do this way from 0 for all x, x', y, y'? Ans: no. We do from 1
            arv_x = ((x + x_1)/2) - x_0 # 
#-----------------------------------------------------------    
            vortex_theta = np.arctan2(arv_y, arv_x)
            
    # how should we apply theta belongs to [0, 2pi) if we want 
            #theta = np.arctan2(y, x)          # (-pi, pi]
            #theta_0_2pi = np.mod(theta, 2*np.pi)  # [0, 2pi)
            # continuous (no jump) if points are ordered around the circle:
            #theta_cont = np.unwrap(theta)
#------------------------------------------------------------
            t_phase_delayed = t - vortex_theta * T/(2 * np.pi)
          

            H_generated_both_OBC_vor_core[2*i:2*i+2, 2*j:2*j+2] = H2NyNx_xyOBC(N_y, N_x, t_phase_delayed, a_0)[2*i:2*i+2, 2*j:2*j+2]
            
            ############################################?????????? >>>>  Should we use H2NyNx_xOBC_only 
            # Now I will change it to: xy OBC - changed 3.18. as H2NyNx_xyOBC -- from H2NyNx_xOBC_only
            
    return H_generated_both_OBC_vor_core


## seems all cells/func include this functions are wrong!!

#%%This cell, include xy OBC condition, only changed func as [H2NyNx_xyOBC] at 3.18. - details see notes inside the cell 
def H2NyNx_xyOBC_No_Vor(N_y, N_x, t, a_0): # H for OBC but no timevortex, exact written in same way with time vor, but no delay 
    
    H_generated_both_OBC_vor_core = np.zeros((2*N_y*N_x, 2*N_y*N_x), dtype =complex)
    
    
    for i in range(N_x* N_y):
        for j in range(N_x* N_y):
            x = (i + 1 )% N_x  #double check later (?) -- x & x' are from 0, to Nx -1
            x_1 = (j + 1 )% N_x 
          #  delta_x = (x - x_1) % N_x #this gives [the smallest integer which near the (x-x1) and can be divided by N_x]'s difference with (x-x_1)
            y = ((i) // N_x) + 1 #double check later (?) -- y & y' are from 0, to Ny -1
            y_1 =( (j) // N_x ) + 1 # y'

           # print (i, j, 'test ij in time vortex matrix')

            arv_y = ((y + y_1)/2)  # here is a question: should we do this way from 0 for all x, x', y, y'? Ans: no. We do from 1
            arv_x = ((x + x_1)/2)  # 
#-----------------------------------------------------------    
            vortex_theta = np.arctan2(arv_y, arv_x)
#------------------------------------------------------------
 #           t_phase_delayed = t - vortex_theta * T/(2 * np.pi)
          
#-----------------------------------MOST IMPORTANT DIFFERENCE HERE!!!----------------------------
            H_generated_both_OBC_vor_core[2*i:2*i+2, 2*j:2*j+2] = H2NyNx_xyOBC(N_y, N_x, t, a_0)[2*i:2*i+2, 2*j:2*j+2]
            
            
            
    return H_generated_both_OBC_vor_core

# check if this cell is the same as H2NyNx_xyOBC itself?  but firstly you can try plot : by this 

#%%

   
#try one H see if this work -- still define time at T/12 and Ny = 4 Nx =3

H_TEST_time_vortex =  H2NyNx_xyOBC_time_vor_core(4, 3, T/12, 1, 2.5 , 2)
# a vague test: test that we did add some shit here: 
#HamNxNy_test = H2NyNx_xyOBC(4, 3, T/12, 1)   

def Ham_same(a, b, rtol=1e-8, atol=1e-10):
    a = np.asarray(a)
    b = np.asarray(b)
    return a.shape == b.shape and np.allclose(a, b, rtol=rtol, atol=atol, equal_nan=True)

# Here they shouldn't be the same 
#%% # here insead of U_full operator for time vor core at origin, U is time vor at any point(x_0,y_0) input
    # (also can be 0 :)
    # most important stuff here for time vortex, full time evolution opt(also include time discretisation)
def U_full_time_vor_core(N_x, N_y, num_time_stages, a_0, x_0, y_0): # which means along lhs to rhs, from (1) earlier/smaller t to (6) later/larger t
                                                    # here x_0 and y_0 denotes the x y positions of time vortex core
    U_list_matrices = []
    t_desicretisation= T /num_time_stages
#    time_stage = num_time_stages
    
    for i in range (num_time_stages):
        
    
        t_i = T - t_desicretisation * (i + 1)
    
        H_i = H2NyNx_xyOBC_time_vor_core (N_x, N_y, t_i, a_0, x_0, y_0)
#expm(-1j* H_i_test* delta_t_test)
        
        U_i = expm(-1j* H_i* (t_desicretisation)) #!!!!! # here using exponential form construct each of U_1 to U_6
        
        U_list_matrices.append(U_i)

        
    U_full = reduce(np.matmul, U_list_matrices) # here adjoint the full time evolution 
        # here U_list_N_matrices = [e^{-iH(T-dt)}]
    return U_full   


def U_full_time_vor_midpoint(N_x, N_y, num_time_stages, a_0):
    x_mid = (1+ N_x)/2
    y_mid = (1+ N_y)/2
    U_full_for_mid_point = U_full_time_vor_core(N_x, N_y, num_time_stages, a_0, x_mid, y_mid)
    
    return U_full_for_mid_point
    
#Test_U_with_time_vor = U_full_time_vor_midpoint(3, 4, 50, 1) # make time steps to be 12 

#%%
def U_full_time_no_vor(N_x, N_y, num_time_stages, a_0): # Recall xyOBC but no vortex
    U_list_matrices = []
    t_descretisation= T /num_time_stages
#    time_stage = num_time_stages
    
    for i in range (num_time_stages):
        
    
        t_i = T - t_descretisation * (i + 1)
    
        H_i =  H2NyNx_xyOBC_No_Vor(N_x, N_y, t_i, a_0)
#expm(-1j* H_i_test* delta_t_test)
        
        U_i = expm(-1j* H_i* (t_descretisation)) #!!!!! # here using exponential form construct each of U_1 to U_6
        
        U_list_matrices.append(U_i)
        

        
    U_full = reduce(np.matmul, U_list_matrices) # here adjoint the full time evolution 
        # here U_list_N_matrices = [e^{-iH(T-dt)}]
    return U_full

#%% here below is an example code  **********Fig 4b reproduce from here!!!    
## IMPORTANT METHOD: when reading codes, please start from: what is input; what is [OUTPUT]
#    return 1
def build_all_sites(N_x, N_y, n_alpha=2):
    """
    All physical sites, stored as rows [x, y, alpha].
    Uses 1-based indexing to match your formulas.
    """
    x = np.arange(1, N_x + 1)  # x will be from 1 to N_x
    y = np.arange(1, N_y + 1)
    alpha = np.arange(1, n_alpha + 1)

    X, Y, A = np.meshgrid(x, y, alpha, indexing='xy')       #Q: What is this? 
    sites = np.column_stack([X.ravel(), Y.ravel(), A.ravel()])#Q: What is this? combine  
    return sites

#all test for above func: 
test_sites33 = build_all_sites(3, 3, n_alpha=2) #build an array of sites just as in [list_of_sites_around_mid_4](but that was list)


#%%
# this part generate the **"[[real basis label]]"** of sites 
# output: column(or row?) of the [real basis label] objects 

def sites_to_comp_indices(sites, N_x, n_alpha=2):  #this function is how we convert each real position labelled sites to components of the |phi_n>
                                                   #which is eigenvec of U(full time evolution)
    """
    Convert rows [x,y,alpha] to basis/component indices.
    """
    #sites = np.asarray(sites, dtype=int) #previous func already build this as ndarray? 

    x = sites[:, 0]
    y = sites[:, 1]
    alpha = sites[:, 2]

    comp_index = (y - 1) * (n_alpha * N_x) + (x - 1) * n_alpha + (alpha - 1)    
    # here it is no problem that no input Ny, because y is the slowest changed index, we only need to know Nx always for that if y>1, Nx will be counted in previous comp_index 
    return comp_index                         #**ALSO, you can consider it as: as the outest layer of UC, y can extented to inf but doesn't affect the structure of how sites got labelled as real positions by |x, y, alpha> 


# But as in that we have written our 'column stack' of all sites in the same way our |phi_n> built
# so input:sites >> func in this cell >> output just from [0,1, ..., 17] ( for Nx x Ny 3x3)

# test of above function: is comp_index an array? a list? corresponding to the 'sites' input? 
#test_sites_as_m_th_comp_of_phi_n = sites_to_comp_indices(test_sites33, 3)
#print(type(test_sites_as_m_th_comp_of_phi_n))
#test result: this is an array as in : X.ravel()|Y.ravel()|A.ravel() -- which(in coding label, start from 0)component of phi (eigenvec of U?)
#                                           ..  | ...     |...       -- component of 
# in comparison
#test_sites_to_comp_not_use_previous_func = sites_to_comp_indices(np.array([[3, 1, 1], [2, 2, 1], [2, 2, 2], [3, 2, 1], [3, 2, 2], [2, 3, 2]]), 4, n_alpha=2) # here we can try to input list_of_sites_around_mid_4 = [[3, 1, 1], [2, 2, 1], [2, 2, 2], [3, 2, 1], [3, 2, 2], [2, 3, 2]]

#print(test_sites_to_comp_not_use_previous_func, 'it should not be 0, 1,..., 16 anymore')
#(type(test_sites_to_comp_not_use_previous_func))

#%% ###need to be changed !!

def Rho_all_sites_epsilon_vor(N_x, N_y, num_time_stages, a_0 , sigma, list_of_sites_needed=None, n_alpha=2):
 #   sort_by_epsilon=True,          #########PROBLEM1#########           **This sigma offset !! Should be exactly include the last 'pi' mode, as in: 
                                                         # for 6-sites summation, we have no peak at for example 2.95. But we don't know about other sites (especially near bound)
 
#    sum_over_filtered_states=True, -- dk why we are using these? 

    U_4 = U_full_time_vor_midpoint(N_x, N_y, num_time_stages, a_0) # real xy OBC matrix with vortex at midpoint
    #############PROBLEM 2################!!!
    # here this number of time stages dirctly influence how many T discretisation you slicing in U full
    # just as T sample in your Fig 4a graph!! To use 6 is really not good 
    
    
    eigvals_U_4, eigvecs_U_4 = np.linalg.eig(U_4) # get eigenvalue
    epsilon_n = -np.angle(eigvals_U_4) # get epsilon; 
#    eigvec_copy_backup = eigvecs_U_4[] no backup but eigenvec and epsilons corresponding to each other as in [filtered_eigvec] we used [:, mask]

    
    if list_of_sites_needed is None: # here None means the sites we extract component array and involved in Rho correspondence is defult as same as the whole plaquette which generate H
        sites_used = build_all_sites(N_x, N_y, n_alpha=n_alpha)
    else:                            # in case if we only what a few sites, like around mid plaquette...so on 
        sites_used = np.asarray(list_of_sites_needed, dtype=int)

    comp_indices = sites_to_comp_indices(sites_used, N_x, n_alpha=n_alpha)
    # sites_used are labelled by quatum nnumber of |x, y, alpha >
    
#    if eps_min <= eps_max:
    eps_min = np.pi - sigma
    eps_max = np.pi + sigma
     
    mask = (epsilon_n >= eps_min) & (epsilon_n <= eps_max)
    
    #mask = np.abs(np.abs(epsilon_n) - np.pi) <= sigma *** This is an important alternative ***
    
#    else:
#        mask = (epsilon_n >= eps_min) | (epsilon_n <= eps_max) # Ordinary interval; if crossing branch cut, use wrapped interval

    filtered_epsilons = epsilon_n[mask]
#    filtered_eigvals = eigvals_U_4[mask]
    filtered_eigvecs = eigvecs_U_4[:, mask]   # shape (dim, N_filtered)
    
    filtered_eigvecs_reorder = filtered_eigvecs
#    if sort_by_epsilon and filtered_epsilons.size > 0: 
    if filtered_epsilons.size > 0:  # there are qausienergies (epsilon_n) in the near-pi range we needed 
        order = np.argsort(filtered_epsilons)  # this [order] is a list of extracted index of filtered near-pi quasiEs
        # how we sort our selected near-pi quasienergies?  #a = np.array([30, 10, 20])
        #idx = np.argsort(a)
        # we will get idx = 1, 2, 0 from SMALLEST to LARGEST
        # NOTE: np.sort return sorted values themselve; np.sort return index of sorted values (from smallest to largest)

        filtered_epsilons = filtered_epsilons[order]                                         
        #filtered_eigvals = filtered_eigvals[order]
        filtered_eigvecs_reorder = filtered_eigvecs[:, order] # filtered eigenvecs as in from eigenvec[:, mask], combined columns
        # here we only reordered eigenvectors, but didn't break the structure of eigvec comp themselves
        # after this line the filtered eigvecs will be reordered in that from SMALLEST epsilon(-pi side(?)) to largest epsilon
        #[filtered_eigvecs_reorder] is an object of that combined eigvecs as **column** of this matrix 

    # rho_by_eps[ comp_index , n -- as in epsilon_n : here we want near -pi] = |phi_n(site_{comp_index})|^2
    rho_all_filtered_all_comp = np.abs(filtered_eigvecs_reorder[comp_indices, :])**2 
    # if we didn't change order of comp correspondence for input sites (this is from how we generate the sites, this is just as filtered_eigvecs_reorder itself)
    # here we will return a matrix of Rho, rows are ordered in same way the comp_index arranged to sites are ordered
    # as in if our input sites are not like generated now, but for example mid_plaqutte sites, we will have
    # filtered _eigvecs_reorder, REORDERED AGAIN that 1st row corresponding to first site, second correspodning to second...but here we don't need this 
    # (comp_indices is just the comp_index returned by previous func; it's an ndarray )
    # but filtered_eigvecs_reorder[comp_indices, :] is a matrix?? How we can use it in abs?
    
    
    
# here a reference of our previous wrote function needed 
# Rho_n_site_i = (np.abs(eigen_vec_n[comp_index])) ** 2 -- here this is sepcific eigvec_n
    #sum_over_filtered_states:
    rho_sum = rho_all_filtered_all_comp.sum(axis=1) # sum across rows, as in different selected (filtered eigvecs), 
                                       # corresponding to that of summing over small interval of epsilon_n ~ pi   # shape (N_sites,)
    #** THIS WILL RETURN an array of Rho in order of all comp_indices 
#    else:
#    rho_sum = None

    return {
        "sites_used": sites_used,                 # (N_sites, 3) STACK OF COLUMN |x|y|alpha> -- useful for next step plot in real space as plaquette 
        "comp_indices": comp_indices,            # (N_sites,)
       # "rho_by_state": rho_by_state,            # (N_sites, N_filtered)
        "rho_sum": rho_sum,                      # (N_sites,)
      #  "filtered_epsilons": filtered_epsilons,  # (N_filtered,)
      #  "filtered_eigvals": filtered_eigvals,    # (N_filtered,)
      #  "filtered_eigvecs": filtered_eigvecs,    # (dim, N_filtered)
    }
    

    #Q1: how would this produce value Q2: speed? if Nx Ny is 8x8?? 
    
    #this will return eigenvectors as in eigenves's mth column, 
    #corresponding to mth eigenval generated at the same time by np.linalg.eig
    # here is a test function to check if we get the corresponding eigenvec for the eigenval (epsilon)

#%%
# TEST AIM : NO_VOR check: try to produce a 3x3 Rho, without vortex, to check if it can give anything on the plaquette plot (with color)
# now for that we need to get num_time_stages = 6 for no vor; 
#  1. how should we choose sigma?  now for 3x3, choose 0.05
# --  trace back (actually from afterwards cell: generate epsilon/s) in function [Rho_sites_all_epsilons_noVor]

 
def Rho_all_sites_epsilon_no_vor(N_x, N_y, num_time_stages, a_0 , sigma, list_of_sites_needed=None, n_alpha=2):
  #   sort_by_epsilon=True,
 #    sum_over_filtered_states=True, -- dk why we are using these? 

     U_4 = U_full_time_no_vor(N_x, N_y, num_time_stages, a_0) # real xy OBC matrix with vortex at midpoint

     eigvals_U_4, eigvecs_U_4 = np.linalg.eig(U_4) # get eigenvalue
     epsilon_n = -np.angle(eigvals_U_4) # get epsilon; 
 #    eigvec_copy_backup = eigvecs_U_4[] no backup but eigenvec and epsilons corresponding to each other as in [filtered_eigvec] we used [:, mask]

     
     if list_of_sites_needed is None: # here None means the sites we extract component array and involved in Rho correspondence is defult as same as the whole plaquette which generate H
         sites_used = build_all_sites(N_x, N_y, n_alpha=n_alpha)
     else:                            # in case if we only what a few sites, like around mid plaquette...so on 
         sites_used = np.asarray(list_of_sites_needed, dtype=int)

     comp_indices = sites_to_comp_indices(sites_used, N_x, n_alpha=n_alpha)
     # sites_used are labelled by quatum nnumber of |x, y, alpha >
     
 #    if eps_min <= eps_max:
     eps_min = np.pi - sigma
     eps_max = np.pi + sigma
      
     mask = (epsilon_n >= eps_min) & (epsilon_n <= eps_max)
     
     #mask = np.abs(np.abs(epsilon_n) - np.pi) <= sigma ** This is an important alternative 
 #    else:
 #        mask = (epsilon_n >= eps_min) | (epsilon_n <= eps_max) # Ordinary interval; if crossing branch cut, use wrapped interval

     filtered_epsilons = epsilon_n[mask]
 #    filtered_eigvals = eigvals_U_4[mask]
     filtered_eigvecs = eigvecs_U_4[:, mask]   # shape (dim, N_filtered)
     
     filtered_eigvecs_reorder = filtered_eigvecs
 #    if sort_by_epsilon and filtered_epsilons.size > 0: 
     if filtered_epsilons.size > 0:  # there are qausienergies (epsilon_n) in the near-pi range we needed 
         order = np.argsort(filtered_epsilons)  # this [order] is a list of extracted index of filtered near-pi quasiEs
         # how we sort our selected near-pi quasienergies?  #a = np.array([30, 10, 20])
         #idx = np.argsort(a)
         # we will get idx = 1, 2, 0 from SMALLEST to LARGEST
         # NOTE: np.sort return sorted values themselve; np.sort return index of sorted values (from smallest to largest)

         filtered_epsilons = filtered_epsilons[order]                                         
         #filtered_eigvals = filtered_eigvals[order]
         filtered_eigvecs_reorder = filtered_eigvecs[:, order] # filtered eigenvecs as in from eigenvec[:, mask], combined columns
         # here we only reordered eigenvectors, but didn't break the structure of eigvec comp themselves
         # after this line the filtered eigvecs will be reordered in that from SMALLEST epsilon(-pi side(?)) to largest epsilon
         #[filtered_eigvecs_reorder] is an object of that combined eigvecs as **column** of this matrix 

     # rho_by_eps[ comp_index , n -- as in epsilon_n : here we want near -pi] = |phi_n(site_{comp_index})|^2
     rho_all_filtered_all_comp = np.abs(filtered_eigvecs_reorder[comp_indices, :])**2 
     # if we didn't change order of comp correspondence for input sites (this is from how we generate the sites, this is just as filtered_eigvecs_reorder itself)
     # here we will return a matrix of Rho, rows are ordered in same way the comp_index arranged to sites are ordered
     # as in if our input sites are not like generated now, but for example mid_plaqutte sites, we will have
     # filtered _eigvecs_reorder, REORDERED AGAIN that 1st row corresponding to first site, second correspodning to second...but here we don't need this 
     # (comp_indices is just the comp_index returned by previous func; it's an ndarray )
     # but filtered_eigvecs_reorder[comp_indices, :] is a matrix?? How we can use it in abs?
     
     
     
 # here a reference of our previous wrote function needed 
 # Rho_n_site_i = (np.abs(eigen_vec_n[comp_index])) ** 2 -- here this is sepcific eigvec_n
     #sum_over_filtered_states:
     rho_sum = rho_all_filtered_all_comp.sum(axis=1) # sum across rows, as in different selected (filtered eigvecs), 
                                        # corresponding to that of summing over small interval of epsilon_n ~ pi   # shape (N_sites,)
     #** THIS WILL RETURN an array of Rho in order of all comp_indices 
 #    else:
 #    rho_sum = None

     return {
         "sites_used": sites_used,                 # (N_sites, 3) STACK OF COLUMN |x|y|alpha> -- useful for next step plot in real space as plaquette 
         "comp_indices": comp_indices,            # (N_sites,)
        # "rho_by_state": rho_by_state,            # (N_sites, N_filtered)
         "rho_sum": rho_sum,                      # (N_sites,)
       #  "filtered_epsilons": filtered_epsilons,  # (N_filtered,)
       #  "filtered_eigvals": filtered_eigvals,    # (N_filtered,)
       #  "filtered_eigvecs": filtered_eigvecs,    # (dim, N_filtered)
     }
 
#test_produceRho3x3 = Rho_all_sites_epsilon_no_vor(3, 3, 6, 1, 0.5)

#%%
#****** UNFINISHED !!!! Must finish today
# here this cell is for plotting the site_used as in real space points on 8x8 plaquette
#plot grids
def site_to_plot_xy(sites, a=1.0):

    sites = np.asarray(sites) # make a clear copy? 
    x = sites[:, 0]
    y = sites[:, 1]
    alpha = sites[:, 2]

    """
    Convert [x, y, alpha] -> real plotting coordinates (Xplot, Yplot).

    You should adapt this to your exact honeycomb convention.
    """
    sites = np.asarray(sites) #? 
    x = sites[:, 0]
    y = sites[:, 1]
    alpha = sites[:, 2]
   
    # Example placeholder mapping only (here we assumed each side of hex is 1 - Q: is this same as a0?)
    #there is a coor_x shift term for y doesn't equal 1(y_coor non zero) x_shift = 0.5 *np.sqrt(3) * (y-1)
    coor_x = (x - 1) * np.sqrt(3) + 0.5 * (alpha - 1) * np.sqrt(3) +  0.5 *np.sqrt(3) * (y-1)
    coor_y = (y - 1) * 1.5 + (alpha -1 ) * 0.5

    return coor_x, coor_y


#%%
#recall function construction: Rho_all_sites_epsilon_vor(N_x, N_y, num_time_stages, a_0 , sigma, list_of_sites_needed=None, n_alpha=2)
#Plaq_dict_UC8X8_vor_mid_offset_12 = Rho_all_sites_epsilon_vor(8, 8, 6, 1, 0.12)
  
#print('offset 0.12', Plaq_dict_UC8X8_vor_mid_offset_12)


#Plaq_dict_UC8X8_vor_mid_offset_18 = Rho_all_sites_epsilon_vor(8, 8, 6, 1, 0.18)

#print('offset 0.18', Plaq_dict_UC8X8_vor_mid_offset_18)


#Plaq_dict_UC8X8_vor_mid_offset_3 = Rho_all_sites_epsilon_vor(8, 8, 6, 1, 0.3)  
#print('offset 0.3', Plaq_dict_UC8X8_vor_mid_offset_3)
#%%
 #VORTEX, PLAQUETTE 
Plaq_dict_UC8X8_vor_mid = Rho_all_sites_epsilon_vor(6, 6, 35, 1, 0.14) #offset need to be setted correctly 
#good news! It runs about 5 min #the offsetadjustion here is a bit tricky we used 0.25 for 3.15 - 2.9
#%%
#WITH VORTEX, PLAQUETTE PLOT
Rho_near_pi8x8_vor = Plaq_dict_UC8X8_vor_mid["rho_sum"]
Xplot, Yplot = site_to_plot_xy(Plaq_dict_UC8X8_vor_mid["sites_used"]) ## here for test, we just use the site plotted as by previous 'build site' function
#    colors = result["rho_sum"] results is?

print('offset 0.12', Plaq_dict_UC8X8_vor_mid)

plt.figure(figsize=(6, 6))
#plt.scatter(Xplot, Yplot, s=80)
plt.scatter(Xplot, Yplot, c=Rho_near_pi8x8_vor, s=180, cmap="inferno_r", linewidths=0.6) # , dgecolors="black"
plt.xlabel("X")
plt.ylabel("Y")
plt.colorbar(label="LDOS near pi")
plt.axis("equal")
plt.grid(True)
plt.show()

#print('print dict of 8x8 with cor, Rho_comp_corresponding, offset 0.3', Plaq_dict_UC8X8_vor_mid)
## here for test, we just use the site plotted as by previous 'build site' function
#    colors = result["rho_sum"] results is?

# YOU CAN ALSO TRY SUMMED OVER THE SITES FOR EACH UNIT CELLS ! BECAUSE IT MIGHT NOT NECESSARILY TO BE SAME FOR 'SITE LOCALISATION' AND 'CELL LOCALISATION'


#%% NO VORTEX, PLAQUETTE 

# create object: dict created by [Rho_all_sites_epsilon_no_vor] 
Plaq_dict_UC8X8 = Rho_all_sites_epsilon_no_vor(6, 6, 35, 1, 0.085) #here if offset 0.12 will get nothing 
# Here what I cahnged is literaterally 
# other reasons check later
#-- later we will return here about running this; but now we use the produced one for speed 

#%% NO VORTEX, PLAQUETTE PLOT
#Plaq_dict_UC8X8_DONE = test_produceRho8x8

Rho_near_pi8x8 = Plaq_dict_UC8X8["rho_sum"]
Xplot, Yplot = site_to_plot_xy(Plaq_dict_UC8X8["sites_used"]) ## here for test, we just use the site plotted as by previous 'build site' function
#    colors = result["rho_sum"] results is?


plt.figure(figsize=(6, 6))
#plt.scatter(Xplot, Yplot, s=80)
plt.scatter(Xplot, Yplot, c=Rho_near_pi8x8, s=100, cmap="inferno_r", linewidths=0.6) # , dgecolors="black"
plt.xlabel("X")
plt.ylabel("Y")
plt.colorbar(label="LDOS near pi")
plt.axis("equal")
plt.grid(True)
plt.show()

# problem: why we didn't get x edge but only get y edge? 

#%% a test for peace of mind: 

# np.sort and np.argsort:

array_being_sorted = np.array([4, 2, 3, 7])
Index_small_to_large = np.argsort(array_being_sorted)
print(Index_small_to_large, 'prediction: 1, 2, 0, 3')  #checks out; also, this output is a list

print(np.sort(array_being_sorted)) # returns the values themselves as from smallest to largest 


# test can we use a list of index to reorder an array?
array_reordered = array_being_sorted[Index_small_to_large] # if put this index list munually, need to be [[]]
print(array_reordered)

# check why this works?
# vector object setting: eig1 = [1, 2, 3, 4] - 1st column as in 0th index; eig2 = [5, 6, 7, 8] - 2nd column, idx:1 ; eig3 = [9, 10, 11, 12] 3rd c, idx 2
column_combined_vec = np.array([[1, 5, 9], [2, 6, 10], [3, 7, 11], [4, 8, 12]]) 

test_order_index = [0, 2, 1]
reorder_combined_vec = column_combined_vec[:, test_order_index ] # here I expec to be eig 1, eig3, eig2 ---
#IMPORTANT : as long as the matrix you are reordering is an array, this works 
#IMPORTANT ** : for matrices, A[:, list] = reorder column of A, but don't touch inside of columns ; A[list, :] reorder rows of A, but don't touch inside of rows  
print(reorder_combined_vec) # ???

# test matrix in np.abs:
test_matrix  = np.array([[1+1j, 3, 5], [2 + 1j, 4, 6]])
abs_matrix = np.abs(test_matrix)**2
print(abs_matrix) #works

# test: sum by axis:
sum_over_row = test_matrix .sum(axis=1) # axis=1: sum down the rows for each row
sum_over_axis = test_matrix .sum(axis=0) # axis=0: sum down the rows for each column

#test if 
A_if_not = 4
print( A_if_not)
if A_if_not> 0:
    A_if_not = 5
print (A_if_not)
#%% MY 'INTENSITY' plot
# -------which sites got summed over --------

list_of_sites_around_mid_4 = [[3, 1, 1], [2, 2, 1], [2, 2, 2], [3, 2, 1], [3, 2, 2], [2, 3, 2]] #each site labelled as [x,y,alpha]
list_of_sites_around_mid_8 = [[5, 3, 1], [4, 4, 1], [4, 4, 2], [5, 4, 1], [5, 4, 2], [4, 5, 2]]                                                                                               # as in alpha is A or B>> 1 or 2
                                                                                               # others from 1 to Nx/y
#print(list_of_sites_around_core[1][1])



#%%
# for sublattices, A is 1, B is 2 (grid index) -- this is showed in eqn (26)>>(28)

def Rho_sites_all_epsilons(N_x, N_y, num_time_stages, a_0, list_of_sites_needed): # this will return a list in order of epsilon_n
    U_4 = U_full_time_vor_midpoint(N_x, N_y, num_time_stages, a_0)   # var: sites summed over
                                                                     # var: epsilon( from NxNy num_time_sites, a_0>>U)    
    eigvals_U_4, eigvecs_U_4 = np.linalg.eig(U_4)
    
    epsilon_n_list = - np.angle(eigvals_U_4) 
    
    
    get_eigen_vecs = []
    for m in range(len(epsilon_n_list)):
        eigen_vec_extract = eigvecs_U_4[:, m]
        get_eigen_vecs. append (eigen_vec_extract)
    eigen_vec_list = get_eigen_vecs
    
    # decompose list_of_sites_needed
    
    Rho_certain_sites_all_n = []
    for eigen_vec_n in eigen_vec_list:
        # for a certain x, y, alpha 
        Rho_n_all_sites = 0
        for site_i in list_of_sites_needed: # here in this loop we are summing over all the site
            
            x, y, alpha = site_i # For alpha's label in sites, we used 1,2 as input instead of 0,1 
            comp_index = (y - 1)* 2 * N_x + (x - 1) * 2 + alpha - 1 # -1 is to fit coding regime 
        #ref: Rho_component_y = (np. abs(Phi_eigvecs_n[2*y]))**2 + (np.abs(Phi_eigvecs_n[2*y + 1]))**2
            Rho_n_site_i = (np.abs(eigen_vec_n[comp_index])) ** 2 
            
            Rho_n_all_sites += Rho_n_site_i 
        
        Rho_n = Rho_n_all_sites
        
        Rho_certain_sites_all_n.append(Rho_n) # this is Rho_sites for all epsilon
        
    return Rho_certain_sites_all_n, epsilon_n_list #should be in order of epsilon_n_list from np.linalg.eig
                                                   #why? ANs: get_eigen_vecs is in the order of epsilo_n_list

#%%
def Rho_sites_all_epsilons_noVor(N_x, N_y, num_time_stages, a_0, list_of_sites_needed): # this will return a list in order of epsilon_n
    U_5 = U_full_time_no_vor(N_x, N_y, num_time_stages, a_0)   # var: sites summed over
                                                                     # var: epsilon( from NxNy num_time_sites, a_0>>U)    
    eigvals_U_5, eigvecs_U_5 = np.linalg.eig(U_5)
    
    epsilon_n_list = - np.angle(eigvals_U_5) #ARRAY
    
    
    get_eigen_vecs = []
    for m in range(len(epsilon_n_list)):
        eigen_vec_extract = eigvecs_U_5[:, m]
        get_eigen_vecs. append (eigen_vec_extract)
    eigen_vec_list = get_eigen_vecs
    
    # decompose list_of_sites_needed
    
    Rho_certain_sites_all_n = []
    for eigen_vec_n in eigen_vec_list:
        # for a certain x, y, alpha 
        Rho_n_all_sites = 0
        for site_i in list_of_sites_needed:
            
            x, y, alpha = site_i
            
            comp_index = (y - 1)* 2 * N_x + (x - 1) * 2 + alpha - 1 # -1 is to fit coding regime 
        #ref: Rho_component_y = (np. abs(Phi_eigvecs_n[2*y]))**2 + (np.abs(Phi_eigvecs_n[2*y + 1]))**2
            Rho_n_site_i = (np.abs(eigen_vec_n[comp_index])) ** 2 
            
            Rho_n_all_sites += Rho_n_site_i 
        
        Rho_n = Rho_n_all_sites
        
        Rho_certain_sites_all_n.append(Rho_n) # this is Rho_sites for all epsilon
        
    return Rho_certain_sites_all_n, epsilon_n_list #should be in order of epsilon_n_list from np.linalg.eig

#%% check func indentation1 
#    def test_eigenpair_alignment(rtol=1e-7, atol=1e-9):
#        check_eigenpair =[]
#        for m in range (2*N_x*N_y):
#            lam = eigvals_U_4[m]
#            v = eigvecs_U_4[:, m] # here eigen vac of np.linalg.eig is an array which columns are eigenvecs 

#            residual = U_4 @ v - lam * v
#            residual_norm = np.linalg.norm(residual) # should be close to zero, *** this can be look later to find inaccuracy of matrix? 

#            passed = np.allclose(U_4 @ v, lam * v, rtol=rtol, atol=atol)
#            
#            check_eigenpair.append((passed, residual_norm))
#        return check_eigenpair # this should return 'true', if not, go back to check     epsilon_n_list & eigen_vec_list

#%% T sample for time discretisation adjusting portal 

t_sampling = 70


#%% plot try scatter # easier to read --WITH TIME VOR
Rho, Epsilon = Rho_sites_all_epsilons(8, 8, t_sampling, 1, list_of_sites_around_mid_8) # it seems the more time sampled, the larger the bulk structure(that range of Epsilon)'s Rho
#print (Rho, type(Rho), 'Rho') # Rho is list
#print (Epsilon, type(Epsilon), 'Epsilon') # Epsilon is array 
idx = np.argsort(Epsilon) # sorted all against the small/large of Epsilon

Eps_sorted = Epsilon[idx]
Rho_sorted = np.array(Rho)[idx]

#%%
#%%
plt.figure()
plt.plot(Eps_sorted, Rho_sorted, 'o-', markersize=3)#linestyle='-', linewidth=1)   # or add marker='.'
plt.xlabel(r'$\epsilon_n T$')
plt.ylabel(r'$\rho$')
plt.title(r'$\rho_n for site around vortex core$ vs $\epsilon_n$ ')
plt.grid(True)
plt.show()
#%% 



# BELOW ARE THE PLOTS FOR FROM UC IS 8X8 TO UC 10X10 TO UC 12X12 



#%% read out from excel file and PLOT -- UC:[ 8x8 ]; TS: 50(?); vor core: Ny+1;Nx+1

df = pd.read_excel(r"D:\OneDrive - Imperial College London\Desktop\master project\if I am lucky\Next-step work\Rho_Eps_more_or_lessdetail.xlsx", header=None, skiprows=1)
#this file with less detail of rho and eps, is used to be substracted by [UC:8x8; no vor; TS: 50(?); vor core: Ny+1;Nx+1] so we assume this file is under same condition with t vor
Rho_readout = df.iloc[:, 0]
Epsilon_readout = df.iloc[:, 1]

peak_extraction = np.asarray(Rho_readout)
i_max = np.argmax(peak_extraction)
x_peak, y_peak = Epsilon_readout[i_max], Rho_readout[i_max]
print('found by old way', x_peak, 'pi mode specific quasienergy', y_peak, 'pi mode specific quasienergy') #**PEAK OF PI MODE LOCATION(extracted by highest peak)**

print(np.array(Epsilon_readout)[-2], 'pi mode specific quasienergy', np.array(Rho_readout)[-2], 'DLOS height')  #**PEAK OF PI MODE LOCATION(extracted by identify near pi epsilon)**


plt.figure()
plt.scatter(x_peak, y_peak, marker='x', color='red')
plt.plot(Epsilon_readout, Rho_readout, 'o-', markersize=3)#linestyle='-', linewidth=1)   # or add marker='.'
plt.xlabel(r'$\epsilon_n T$ ')
plt.ylabel(r'$\rho$')
plt.title(r'$\rho_n$ for site around vortex core vs $\epsilon_n$ ')
plt.grid(True)
plt.show()


#%% convert txt file to excel -- with vor; UC: [10x10]; TS: 60; vor core: Ny;Nx (which is a bit off)

import re

file_path = Path(r"D:\OneDrive - Imperial College London\Desktop\master project\if I am lucky\Next-step work\Mimi_op_before_inter.txt")
text = file_path.read_text(encoding="utf-8")

# match: [eps block] Eps_sorted [rho block] Rho_sorted
m = re.search(
    r'(\[[^\[\]]*?\])\s*Eps_sorted\s*(\[[^\[\]]*?\])\s*Rho_sorted',
    text,
    re.S
)

if not m:
    raise ValueError("Could not find the Eps_sorted / Rho_sorted blocks in the txt file.")

eps_block = m.group(1)
rho_block = m.group(2)

Eps_readouttxt_10 = np.fromstring(eps_block.strip("[]"), sep=' ')
Rho_readouttxt_10 = np.fromstring(rho_block.strip("[]"), sep=' ')

print("len(Eps_sorted) =", len(Eps_readouttxt_10))
print("len(Rho_sorted) =", len(Rho_readouttxt_10))

print("first 5 eps:", Eps_readouttxt_10[:5])
print("first 5 rho:", Rho_readouttxt_10[:5])

if len(Eps_readouttxt_10) != len(Rho_readouttxt_10):
    raise ValueError("Eps_sorted and Rho_sorted have different lengths.")

df = pd.DataFrame({
    "Eps_sorted": Eps_readouttxt_10,
    "Rho_sorted": Rho_readouttxt_10
})

out_path = r"D:\OneDrive - Imperial College London\Desktop\master project\if I am lucky\Next-step work\T_sampleunknow_10x10_new_Ny_add_1_or_not.xlsx"
df.to_excel(out_path, index=False)

print("Saved to:", out_path)

#%% read out above(excel)file and PLOT -- with vor; UC: [10x10]; TS: 60; vor core: Ny;Nx (which is a bit off)

df = pd.read_excel(r"D:\OneDrive - Imperial College London\Desktop\master project\if I am lucky\Next-step work\T_sampleunknow_10x10_new_Ny_add_1_or_not.xlsx", header=None, skiprows=1)

Epsilon_readout_10uc = df.iloc[:, 0]
Rho_readout_10uc = df.iloc[:, 1]

#E_cut_edge = Epsilon_readout_12uc[1:287]
#R_cut_edge = Rho_readout_12uc[1:287]
E_arry_10 = np.array(Epsilon_readout_10uc)
R_arry_10 = np.array(Rho_readout_10uc)

print(E_arry_10[-2], 'pi mode specific quasienergy', R_arry_10[-2], 'DLOS height')  #**PEAK OF PI MODE LOCATION** 
#print('is last mode after pi', E_arry_10[-1], 'pi mode specific quasienergy', R_arry_10[-1], 'DLOS height') -- still smaller than pi;
# Q: why always the second nearest value before pi get Rho peak, as more likely to be 'pi mode'?  -- leave to future 

plt.figure()
plt.scatter(E_arry_10[-2], R_arry_10[-2], marker='x', color='red') #**PEAK OF (+)PI MODE LOCATION** 
plt.scatter(E_arry_10[1], R_arry_10[1], marker='x', color='red')   #**PEAK OF (-)PI MODE LOCATION**
plt.plot(Epsilon_readout_10uc, Rho_readout_10uc, 'o-', markersize=3)#linestyle='-', linewidth=1)   # or add marker='.'
plt.xlabel(r'$\epsilon_n T$ ')
plt.ylabel(r'$\rho$')
plt.title(r'$\rho_n$ for site around vortex core vs $\epsilon_n$ ')
plt.grid(True)
plt.show()

#%%
# read txt peak (2.299454260001524, 0.33087090065525415)

# --- Eps ---
text_Eps_readout = Path(r"D:\OneDrive - Imperial College London\Desktop\master project\if I am lucky\Next-step work\T_sample50_12x12.txt").read_text(encoding="utf-8")
text_Eps_readout = text_Eps_readout.replace("Eps_sorted", "")
text_Eps_readout = text_Eps_readout.replace("[", " ").replace("]", " ")
Eps_readouttxt = np.fromstring(text_Eps_readout, sep=' ')

# --- Rho ---
text_Rho_readout = Path(r"D:\OneDrive - Imperial College London\Desktop\master project\if I am lucky\Next-step work\T_sample50_12x12_Rho.txt").read_text(encoding="utf-8")
text_Rho_readout = text_Rho_readout.replace("Rho_sorted", "")
text_Rho_readout = text_Rho_readout.replace("[", " ").replace("]", " ")
Rho_readouttxt = np.fromstring(text_Rho_readout, sep=' ')

# check lengths
print(len(Eps_readouttxt))
print(len(Rho_readouttxt))

# make one dataframe with two columns
df = pd.DataFrame({
    "Eps_sorted": Eps_readouttxt,
    "Rho_sorted": Rho_readouttxt
})

# save once
df.to_excel(r"T_sample50_12x12_new.xlsx", index=False)

#%% read out above(excel)file and PLOT -- with vor; UC: [12x12]; TS: 50; vor core: Ny+1;Nx+1

df = pd.read_excel(r"D:\OneDrive - Imperial College London\Desktop\master project\if I am lucky\Next-step work\T_sample50_12x12_new.xlsx", header=None, skiprows=1)

Epsilon_readout_12uc = df.iloc[:, 0]
Rho_readout_12uc = df.iloc[:, 1]

#E_cut_edge = Epsilon_readout_12uc[1:287]
#R_cut_edge = Rho_readout_12uc[1:287]
E_arry = np.array(Epsilon_readout_12uc)
R_arry = np.array(Rho_readout_12uc)

print(E_arry[-2], 'pi mode specific quasienergy', R_arry[-2], 'DLOS height')  #**PEAK OF PI MODE LOCATION**

plt.figure()
plt.scatter(E_arry[-2], R_arry[-2], marker='x', color='red') #**PEAK OF (+)PI MODE LOCATION** 
plt.scatter(E_arry[1], R_arry[1], marker='x', color='red')   #**PEAK OF (-)PI MODE LOCATION**
plt.plot(Epsilon_readout_12uc, Rho_readout_12uc, 'o-', markersize=3)#linestyle='-', linewidth=1)   # or add marker='.'
plt.xlabel(r'$\epsilon_n T$ ')
plt.ylabel(r'$\rho$')
plt.title(r'$\rho_n$ for site around vortex core vs $\epsilon_n$ ')
plt.grid(True)
plt.show()


#%%
# read txt peak (2.299454260001524, 0.33087090065525415)

# --- Eps ---
text_Eps_readout_12_58 = Path(r"D:\OneDrive - Imperial College London\Desktop\master project\if I am lucky\Next-step work\T_sample58_12x12_Eps.txt").read_text(encoding="utf-8")
text_Eps_readout_12_58 = text_Eps_readout_12_58.replace("Eps_sorted", "")
text_Eps_readout_12_58 = text_Eps_readout_12_58.replace("[", " ").replace("]", " ")
Eps_readouttxt_12_58 = np.fromstring(text_Eps_readout_12_58, sep=' ')

# --- Rho ---
text_Rho_readout_12_58 = Path(r"D:\OneDrive - Imperial College London\Desktop\master project\if I am lucky\Next-step work\T_sample58_12x12_Rho.txt").read_text(encoding="utf-8")
text_Rho_readout_12_58 = text_Rho_readout_12_58.replace("Rho_sorted", "")
text_Rho_readout_12_58 = text_Rho_readout_12_58.replace("[", " ").replace("]", " ")
Rho_readouttxt_12_58 = np.fromstring(text_Rho_readout_12_58, sep=' ')

# check lengths
print(len(Eps_readouttxt_12_58))
print(len(Rho_readouttxt_12_58))

# make one dataframe with two columns
df = pd.DataFrame({
    "Eps_sorted": Eps_readouttxt_12_58,
    "Rho_sorted": Rho_readouttxt_12_58
})

# save once
df.to_excel(r"T_sample58_12x12.xlsx", index=False)

#%% read out above(excel)file and PLOT -- with vor; UC: [12x12]; TS: 58; vor core: Ny+1;Nx+1

df = pd.read_excel(r"D:\OneDrive - Imperial College London\Desktop\master project\if I am lucky\Next-step work\T_sample58_12x12.xlsx", header=None, skiprows=1)

Epsilon_readout_12uc_58 = df.iloc[:, 0]
Rho_readout_12uc_58 = df.iloc[:, 1]

#E_cut_edge = Epsilon_readout_12uc[1:287]
#R_cut_edge = Rho_readout_12uc[1:287]
E_arry_1258 = np.array(Epsilon_readout_12uc_58)
R_arry_1258 = np.array(Rho_readout_12uc_58)

print(E_arry_1258[-2], 'pi mode specific quasienergy', R_arry_1258[-2], 'DLOS height')  #**PEAK OF PI MODE LOCATION**

plt.figure()
plt.scatter(E_arry_1258[-2], R_arry_1258[-2], marker='x', color='red') #**PEAK OF (+)PI MODE LOCATION** 
plt.scatter(E_arry_1258[1], R_arry_1258[1], marker='x', color='red')   #**PEAK OF (-)PI MODE LOCATION**
plt.plot(Epsilon_readout_12uc_58, Rho_readout_12uc_58, 'o-', markersize=3)#linestyle='-', linewidth=1)   # or add marker='.'
plt.xlabel(r'$\epsilon_n T$ ')
plt.ylabel(r'$\rho$')
plt.title(r'$\rho_n$ for site around vortex core vs $\epsilon_n$ ')
plt.grid(True)
plt.show()



#%% WITH OUT TIME VOR
Rho_0, Epsilon_0 = Rho_sites_all_epsilons_noVor(8, 8, 6, 1, list_of_sites_around_mid_8) 
# we run this now only for checking where does epsilon locate ( quasiE near pi ) and set the simga in Rho func 

#%%  #[[2]]
idx = np.argsort(Epsilon_0)

Eps_sorted_0 = Epsilon_0[idx]
Rho_sorted_0 = np.array(Rho_0)[idx]

plt.figure()
plt.plot(Eps_sorted_0, Rho_sorted_0, 'o-', markersize=3)#linestyle='-', linewidth=1)   # or add marker='.'
plt.xlabel(r'$\epsilon_n T$ ')
plt.ylabel(r'$\rho$')
plt.title(r'$\rho_n$  for sites around vortex core vs $\epsilon_n$ ')
plt.grid(True)
plt.show()

#%% plot difference between Rho_0 & Rho [[2]] T dis 50

Rho_difference = Rho_readout - Rho_sorted_0 # as in Rho with vortex minus Rho without vortex #readout file has T discretisation 50 

peak_extraction = np.asarray(Rho_difference)
i_max = np.argmax(peak_extraction)
x_peak, y_peak = Epsilon_readout[i_max], Rho_difference[i_max]

print('peak', (x_peak, y_peak))


plt.figure()

plt.plot(Epsilon_readout, Rho_difference, 'o-', markersize=3, color = 'tab:blue')#linestyle='-', linewidth=1)   # or add marker='.'
plt.plot(x_peak, y_peak, marker='x', markersize=10, mew=2,
         color='tab:red', label='peak')
plt.xlabel(r'$\epsilon_n T$ ')
plt.ylabel(r'$\rho$') 
plt.title(r'$\rho_n$ difference for sites around vortex core  vs $\epsilon_n$ ')
plt.grid(True)
plt.show()
#%%[[3]] T dis 70 --------------final plot -------------------
Rho_difference_adjusting = Rho_sorted - Rho_sorted_0 # as in Rho with vortex minus Rho without vortex 

peak_extraction = np.asarray(Rho_difference_adjusting)
i_max = np.argmax(peak_extraction)
x_peak, y_peak = Epsilon_readout[i_max], Rho_difference_adjusting[i_max]

print('peak', (x_peak, y_peak))


plt.figure()

plt.plot(Eps_sorted, Rho_difference_adjusting, 'o-', markersize=3, color = 'tab:blue')#linestyle='-', linewidth=1)   # or add marker='.'
plt.plot(x_peak, y_peak, marker='x', markersize=10, mew=2,
         color='tab:red', label='peak')
plt.xlabel(r'$\epsilon_n T$ ')
plt.ylabel(r'$\rho$') 
plt.title(rf'$\rho_n$ difference for sites around vortex core vs $\epsilon_n$, $T discretisation = {t_sampling}$') #rf"Phase bands with edge state density diagnosis, y OBC, $L_y = {Ny_plot_yOBC_edge}$
plt.grid(True)
plt.show()

#%% test np.argsort
x_tail = [6,2,3,4]
y_tail = np.array([7,3,4,5])
dog_tail = np.argsort(y_tail) # returns [in terms of index of elements] smallest to largest; in terms of array  
print(dog_tail, type(dog_tail), 'dog tail')
test_tail = y_tail [np.array([0,  3,  2,  1])]
print(test_tail) #[3] and [1] switched 

#%% ploy try bin- like (still need to figure out what tthey did here)
Rho, Epsilon = Rho_sites_all_epsilons(4, 4, 50, 1, list_of_sites_around_mid_4)
#Epsilon = Epsilon
Rho = np.array(Rho)

bins = 20
edges = np.linspace(-np.pi, np.pi, bins+1)
centers = 0.5*(edges[:-1] + edges[1:])

which = np.digitize(Epsilon, edges) - 1
rho_binned = np.array([Rho[which == i].mean() if np.any(which == i) else np.nan
                       for i in range(bins)])
# -------------------------------------------------------------

plt.figure()
plt.plot(centers, rho_binned, linestyle='-', linewidth=1)
plt.xlabel(r'$\epsilon$')
plt.ylabel(r'$\langle \rho \rangle_{\mathrm{bin}}$')
plt.title('Binned (DOS-like) curve')
plt.grid(True)
plt.show()


#%% test loop append:
hog_list = []
x_test = [1,2,3]
for x_dog in x_test:
     x_hog = x_dog + 3
     hog_list.append(x_hog)  
#print(hog_list) #----- append is in order of the list 
    
            
            

    









#%% generate spectrum(?) epsilon_n --- calcualte xx density -- xx refer to density at a certain position range

#%% CHECK IF XPBC Y OBC MAKE SENSE ---- BELOW JUST FOR CHECKING --------
def H2NyNx_yOBC_for_check(N_y, N_x, t, a_0): 
    
    H_PBC_generated_y_OBC = H2NyNx_PBC_version2(N_y, N_x, t,  a_0)

#    rows, cols =H_PBC_generated .shape #rows and cols are just vars - like rows = 15
    
    I, J = np.indices(H_PBC_generated_y_OBC.shape)


    
    big_block_top_right_i = (I >= 0) & (I < 2*N_x) 
    big_block_top_right_j = (J >= 2* N_x * (N_y -1)) & (J < 2* N_x * N_y) #Nx(Ny-1)) to j = NxNy
    big_block_bottom_left_i = (I >= 2* N_x * (N_y -1)) & (I < 2* N_x * N_y)
    big_block_bottom_left_j = (J >= 0) &(J < 2*N_x)
    
    mask = ((big_block_top_right_i) & ( big_block_top_right_j )) | (big_block_bottom_left_i) & ( big_block_bottom_left_j )   # both i-rule and j-rule must hold
    
    H_PBC_generated_y_OBC [mask] = 0    
       # or whatever values you want
    print("\nMask (True means 'will be changed'):\n", mask)

    return H_PBC_generated_y_OBC

# for checking spectrum, here is the H_6_2NyNx_OBC
def H_6_2Ny_yOBC_check(N_y, N_x, a_0): #here this function should have same spectrum(EIGEN_yVALS)as the [K BASIS H2x2]
                    #6 refers as list of 6 piecewise H 
    
    H2NyNx_yOBC_6_matrices = []
    
    for t_x in range (6): #x refers index for t
    
        t_piecewise = T/6 *  t_x + T/(6 * 2)   # here this x is equivalent to i in the 2x2 H list func
        
        H2NyNx_OBC_piecewise =  H2NyNx_yOBC_for_check(N_y, N_x, t_piecewise, a_0)
        
        H2NyNx_yOBC_6_matrices.append(H2NyNx_OBC_piecewise) # generate 6 2N_y real H
        
    return H2NyNx_yOBC_6_matrices

def Ufull_descending_yOBC_for_check(N_y, N_x, a_0): # this is for whole period T
# var index: 3
    U_list_1 = []
#reference :
#    H_list_6_2N_y_edge = H_6_2N_y_OBC(N, k_x)
    H_6_2NyNx_yOBC_list = H_6_2Ny_yOBC_check(N_y, N_x, a_0)
    

    for i in range (6):

        U_i = expm(-1j* (H_6_2NyNx_yOBC_list[i])* (T/6))
        
        U_list_1.append(U_i) # append should be in the loop scope, 
                            # one grid inside the for loop title if you want to append the thing produced once loop (per i)
    
    U_reverse_3 = np.flip(U_list_1)
    
    U_full_3 = reduce(np.matmul,  U_reverse_3)
        
    return U_full_3


def find_spectrum_for_yOBC_check(N_y, N_x, a_0):
    eigvals, eigvecs = np.linalg.eig(Ufull_descending_yOBC_for_check(N_y, N_x, a_0))
    epsilon_3 = - np.angle(eigvals) 
    return  epsilon_3


# here the U of y real basis x PBC Ham spectrum needed :
    
#H_6_2Ny_OBC(N_y, k_x) is for call y OBC x PBC ham_6
#U for y OBC x PBC? U_full_ascending_edge(N_y, k_x): # which means along lhs to rhs, from (1) earlier/smaller t to (6) later/larger t
    

def find_spectrum_realy_OBC (N_y, k_x): # For the spectrum of only y is wirtten in real basis
    Ufull_edge_kx_list =[]              # Or saying: H_2Ny (fixed kx then append through)
    epsilon_n_list_of_kx_list = []
#    U_T = Ufull_descending_PBC(N_y, k_x)
    for i in range (len(k_x)): # This loop's each round run one specific k_x
        
        Ufull_edge_kx_list.append(U_full_ascending_edge(N_y, k_x[i]))

        eigvals, eigvecs = np.linalg.eig(Ufull_edge_kx_list[i]) # will this really return anything? 
        # YES! YOU CAN EXTRACT THE NEW APPENDED LIST IN SAME SCOPE
        
        eigvals_fixed_kx_list = eigvals # just rename
#        print(len(eigvals_fixed_kx_list), eigvals_fixed_kx_list,'eigvals of fixed k_x of Ufull')
#        print(eigvals_kx_2Ny_list, 'eigevals')
        epsilon_n_fixed_kx = [] # phi as phase bands labelled by n
        
        for j in range (len(eigvals_fixed_kx_list)): # j labels phase n, epsilon_n(double check at home)
            
            epsilon_kx_j = - np.angle(eigvals_fixed_kx_list[j])
            
            epsilon_n_fixed_kx. append (epsilon_kx_j) # this is list of phase for specific k_x, length 2N
            
#        print(phase_kx_2Ny_list, 'phases')
       
        epsilon_n_list_of_kx_list. append(epsilon_n_fixed_kx) # this will be a list of different k_x, for each k_x has j
                                                     # -- layer one: len(k_x) -- layer two(each k_x), j_phase
        spectrum_flat = [x for sublist in epsilon_n_list_of_kx_list for x in sublist]

    return  spectrum_flat  

N_y_TEST = 15
N_x_TEST = 20
#k_x_TEST = [1 * 2* np.pi/4, 2  *2*  np.pi/4, 3 *2*  np.pi/4, 2*np.pi]
k_x_TEST = generate_k_x_list(0, 2*np.pi, N_x_TEST)
#spectrum_test_realxy = find_spectrum_NxNy(N_y_TEST, N_x_TEST, a_0)

#spectrum_test_realxy_rewrite = find_spectrum_NxNy_rewrite(N_y_TEST, N_x_TEST, a_0)

spectrum_test_realy_yOBC =  find_spectrum_realy_OBC(N_y_TEST, k_x_TEST)
spectrum_yOBC_NyNxCheck = find_spectrum_for_yOBC_check(N_y_TEST, N_x_TEST, a_0) #dim: 2NxNy

print('BEGIN: spectrum_test_realxy:',spectrum_yOBC_NyNxCheck, 'spectrum_test_realy', spectrum_test_realy_yOBC, 'END')
#%%


def same_values_ndp_yes_no(a, b, ndp):
    a = np.asarray(a).ravel()
    b = np.asarray(b).ravel()

    a_rounded = np.round(a, ndp)
    b_rounded = np.round(b, ndp)

    return "yes" if np.allclose(np.sort(a_rounded), np.sort(b_rounded)) else "no"

if_sepctrum_are_the_same_y0BC = same_values_ndp_yes_no(spectrum_test_realy_yOBC, spectrum_yOBC_NyNxCheck,5)
print('if_sepctrum_are_the_same', if_sepctrum_are_the_same)


                           
                          
#%%
# TEST CELL FOR LEARNING MASK

#example of rewrite some entrices in a matrix:
# 1) Suppose A already exists
A = np.arange(1, 26).reshape(5, 5)   # a 5x5 matrix with values 1..25
#print("Original A:\n", A)

# 2) Build index grids I, J (same shape as A)
I, J = np.indices(A.shape) #I[i, j] = i and J[i, j] = j
# I contains row indices, J contains column indices
#print( 'begin', I, J, 'np.indices result')
# 3) Define a rule on (i, j): here, pick entries where (i + j) is divisible by 3
mask = (I + J) % 3 == 0
# -- here for any matrix we can apply this '==' check
B = np.array([[5, 0],
              [3, 5]])
B == 5

#print("\nMask (True means 'will be changed'):\n", mask)

# 4) Assign new values ONLY where mask is True
# Example value rule: set A[i, j] = 100 + 10*i + j
A[mask] = 100 + 10*I[mask] + J[mask]

#print("\nUpdated A:\n", A)
#%%

#Plan : [H2NyNx_xyOBC(N_y, N_x, t, a_0)] >> ***(II)*** Task:  introduce time vortex 
#         Task: time discretisation   -- meeting 12.05
#         H_not 6 but discretised time H_(delta t)_list(?)>> U(T)
#         
#       ***(III)*** generate spectrum of this U(T), plot their own Rho_[of what?] along them,
#                   -- meeting 1.09
#       Q as break points: how can Rho be labelled in positions? 
#       Q: after we figured this out,
#       we can try to generate sum of rho in real positions around 0 for a few points 
#       check with graph on article fig_4_a (**final aim for numerical progress make this weekend)



