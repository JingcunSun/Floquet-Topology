
import numpy as np
from functools import reduce
from scipy.linalg import expm
import matplotlib.pyplot as plt

# Define function J_tilde_x(t), J_tilde_y, J_tilde_z, they are all step functions with same step length of t. 
    # -Physical meaning：J_x,y,z (u_{I, J})is the strength of (nearest neighbour) hopping. 
    # -Here when J_tilte_a(t) is a function of time, the 'tilte' also encompass a pair (?) of u_{I, J} choice. 


#J_tilte_a(t) are all defined in between 0 and T, and t is independent variable - here you can ask later about 
#-- how to distinguish t and delta_t, J_a



#def J_tilte_a(t, t_start, t_end, J_a): #delta_t = t_end - t_start
# parameter of I_tilte_a: start point t_a (a for x, y, z), pulse length delta_T, (l, T), amp J_a  
  
# this J_tilte_a function is based on the paragraph above eqn(16)
def J_tilte_a(t, t_a, pulse_length, l, T, pulse_amp): # t_a define the START POINT in each period T for pulse a
    
    t_start_a = l * T + t_a
    
    t_end_a = t_start_a + pulse_length
    
    if t_end_a < T:
    
#    if t_start_a < t < t_end:        
#        return  pulse_amp   
#    else:
#        return 0 -- this if else is not suitable when t is an array

        return np.where((t >= t_start_a) & (t < t_end_a), pulse_amp, 0)
    
    if t_end_a > T:
        
        t_end_a_new = t_end_a - T
        
        return np.where((t >= t_start_a) & (t <= T) | (t >= 0) & (t < t_end_a_new), pulse_amp, 0)
    
#parameter adjustion:

#(here always take T = 1, l = 0, set the numerical platform of time, also unify the y axis as epsilon in fig_1)
T = 1
l = 0
#(as said in the paragraph above (16), t_x = 0, t_y = T/3, t_z = 2*T/3)
t_x = 0
t_y = T/3
t_z = T * (2/3)
    
delta_t = T/2

J_a = ((0.9 * np.pi )/4)/ (T * 0.5)


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

#%% 
#test: matrices multiplication
#A = np.array([[0.0, 1.0], [1.0, 0.0]])
#B = np.array([[1.0, 0.0], [0.0, 1.0]])
#C = np.array([[0.0, 1.0], [0.0, 1.0]])
#matrices_list = [A, B, C]
#result_matrices_multiplication = reduce(np.matmul, matrices_list)
# test successful, so matrix multiplication in codes is: 1. build list 2. reduce func (good that we deal with 2d)

# a template of writing U
#H_i_test = np.array([[0.0, 1.0], [1.0, 2.0]])
#delta_t_test = T/6
#U_i = expm(-1j* H_i_test* delta_t_test) # remember to check the ,. in matrices, also, j isn't im, only '1j'

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

#%%
# try plot: t is from 0 to 1, J_tilte_x, y, z 

t_test = np.linspace(0, 1, 100)

plt.plot(t_test, J_tilte_x(t_test), label='J_tilte_x')
plt.plot(t_test, J_tilte_y(t_test), label='J_tilte_y')
plt.plot(t_test, J_tilte_z(t_test), label='J_tilte_z')

plt.xlabel('t')
plt.ylabel('J_tilte_a')
plt.title('Plot of J_tilta_a(pulse of hopping amp) vs t')
plt.legend()
plt.grid(True)
plt.show()

#plot succeed
#%%

# lattice length a_0:
lattice_length = 1
stage_length = T/6 #(numerical calculations)
# due to equation (28)
#stage_number = 6 #(just N afterwards)
#Pauli matrices:
pauli_x = np.array([[0.0, 1.0], [1.0, 0.0]])
pauli_y = np.array([[0.0, -1j], [1j, 0.0]])
pauli_z = np.array([[1.0, 0.0], [0.0, -1.0]])
    

#delta_t is the actual T/2, pulse length of our case 
def H_matrix_2by2(k_x, k_y, t, a_0):
    
    def J_tilte_x(t):
        return J_tilte_a(t, t_x, delta_t, l, T, J_a)
    def J_tilte_y(t):
        return J_tilte_a(t, t_y, delta_t, l, T, J_a)
    def J_tilte_z(t):
        return J_tilte_a(t, t_z, delta_t, l, T, J_a) 
    
    H = (J_tilte_x(t) + J_tilte_y(t) * np.cos (k_x * a_0) + J_tilte_z(t) * np.cos(k_y * a_0)) * pauli_y + (J_tilte_y(t) * np.sin (k_x * a_0) + J_tilte_z(t) * np.sin(k_y * a_0)) * pauli_x
    
    
    return H

# generate a list of t which takes 6 samples of static H(1 to 6)

#t = T/6 * i + T/12 (T/12 is for if there is any confusion in the boundary points of the six stages)


# run a small test for the future 2N 2N matrix
H_test_for_entrices = H_matrix_2by2(np.pi/2, np.pi/2, 5/12, 1)
a_11 = H_test_for_entrices[0][0]
a_12 = H_test_for_entrices[0][1]
#(1.413716694115407-1.413716694115407j) (1.413716694115407-1.413716694115407j)
#%%
def H_list_6_matrices(k_x, k_y, N):
    
    H_list_6_matrices =[]
    
    for i in range (N):
    
        t_i = T/N * i + T/(N * 2)
    
        H_list_6_matrices.append(H_matrix_2by2(k_x, k_y, t_i, lattice_length)) # add lattice length and separatable t
    
    return H_list_6_matrices


H_6 = H_list_6_matrices(np.pi/3, np.pi/4, 6)
print('check')
print('here is H as a whole', H_6)

#%%

#def U_list_all_stages(k_x, k_y, N):
    
#    H_list_6_matrices =[]
#    U_list_6_matrices =[]
    
#    for i in range (N):
    
#        t_i = (T/N) * i + T/(N * 2)
    
#        H_i = H_matrix_2by2(k_x, k_y, t_i, lattice_length)
#expm(-1j* H_i_test* delta_t_test)
        
#        U_i = expm(-1j* H_i* (T/N)) #!!!!!
        
#        U_list_6_matrices.append(U_i)
        
#    return U_list_6_matrices


#%%

def U_full_ascending(k_x, k_y, N): # which means along lhs to rhs, from (1) earlier/smaller t to (6) later/larger t
    
#    H_list_6_matrices =[]
    U_list_N_matrices =[]
    
    for i in range (N):
    
        t_i = (T/N) * i + T/(N * 2)
    
        H_i = H_matrix_2by2(k_x, k_y, t_i, lattice_length)
#expm(-1j* H_i_test* delta_t_test)
        
        U_i = expm(-1j* H_i* (T/N)) #!!!!!
        
        U_list_N_matrices.append(U_i)
        
#        print(U_list_N_matrices)
        
    U_full = reduce(np.matmul, U_list_N_matrices)
        
    return U_full

#%%
#REMEMBER: same name of variable doesn't affect each other if they are defined in different functions
#def U_full_descending(k_x, k_y, N): # which means along lhs to rhs, from (6) later/larger t to (1) earlier/smaller t
    
#    H_list_6_matrices =[]
#    U_list_6_matrices =[]
    
#    for i in range (N):
    
#        t_i = T - (T/N) * i + T/(N * 2)
    
#        H_i = H_matrix_2by2(k_x, k_y, t_i, lattice_length)
#expm(-1j* H_i_test* delta_t_test)
        
#        U_i = expm(-1j* H_i* (T/N))
        
#        U_list_6_matrices .append(U_i)
        
#    U_full = reduce(np.matmul, U_list_6_matrices)
        
#    return U_full

#%% Ascending:

#U_full_matrix_ascending_sample = U_full_ascending((np.pi)/2, np.pi/2, 6)

#eigen_val_U_asc_T_sample, eigen_vec_U_asc_T_sample= np.linalg.eig(U_full_matrix_ascending_sample)

#phase_T_asc_sample = [np.log(eigen_val_U_asc_T_sample[0]) *1j, np.log(eigen_val_U_asc_T_sample[1]) *1j]

#depends on the choices of time, after this we will fix time 

#%% next step: find the phase defined by U_full_matrix (print out eigen values)
# -- what is the kx, ky? Is that here the time t should just be T? 


#%%test: looking for eigen states/ eigen values of matrices

#A = np.array([[1.0, 2.0], [3.0, 4.0]])
    
#energy_levels, energy_states = np.linalg.eig(A)

#print(energy_levels, '-eigen value.', energy_states, '-eigen vector.')
        
#%%

#Descending:

#U_full_matrix_descending = U_full_descending((np.pi)/2, (np.pi)/2, 6)

#eigen_val_U_des_T , eigen_vec_U_des_T  = np.linalg.eig(U_full_matrix_descending)

#phase_T_des = [np.log(eigen_val_U_des_T[0]) *1j, np.log(eigen_val_U_des_T[1]) *1j]

# a question to find out later: are the eigen values all real?

#%% (firstly try ascending, plot the phase diagram)

#define k_x, fix k_y as in mp.pi/4(just have a try)

#!!! Ask how to write a function which used another function mentioned before?

k_x_list = np.linspace(0,  2 * np.pi, 100)

#k_y_special = np. pi
k_y_test = np.linspace(0, 2 * np.pi, 50)

# the first mission is to create different plots of phase against k_x as a family of plots
#%%
#k_y_test_special = np.pi /2

#Ufull_matrix_asc_list =[]
#eigval_U_asc_T_list = []
#eigvec_U_asc_T_list = []
#eigen_U_asc_T_list = []
#phase1_T_asc_list = []
#phase2_T_asc_list = []

#for i in range (len(k_x_list)):
    
    
#    Ufull_matrix_asc_list.append(U_full_ascending(k_x_list[i], k_y, 6))
#    print(i)
#    print(Ufull_matrix_asc_list)
    
    
#for i in range (len(k_x_list)):
    
#    eigvals, eigvecs = np.linalg.eig(Ufull_matrix_asc_list[i])
#    print(eigvals,'--eigen value;', eigvecs,'--eigen vector', i,'i') # all these eigen values are complex!!
    
    
#    eigval_U_asc_T_list.append(eigvals)
#    eigvec_U_asc_T_list.append(eigvecs)
    
#eigval_U_asc_T_array = np.array(eigval_U_asc_T_list)
# phase corresponding to 1st eigen value (eigvals)
#    phase1_T_asc_list.append (np.log(eigvals[0]) *1j) -- eigvals are complex, again

#    phase1_T_asc_list.append(- np.angle(eigvals[0]))
#    phase2_T_asc_list.append(- np.angle(eigvals[1]))
#    print(eigvals,'--eigen value;') 
#    print(eigvecs,'--eigen vector')
#    print(i,'i;')
#    print(np.angle(eigvals[0]),'phi;')
#    print('-----virtual line----')

#plt.plot(k_x_list, phase1_T_asc_list, label='phase 1(1st eigval), phi_n = phi_0')
#plt.plot(k_x_list, phase2_T_asc_list, label='phase 2(2nd eigval), phi_n = phi_1')

#plt.xlabel('k_x * a_0')
#plt.ylabel('phase(epsilon * T)')
#plt.title('Plot of phase vs k_x')
#plt.legend()
#plt.grid(True)
#plt.show()

#%%
#!!! Here please tidy up a function out! 

plt.figure()

def plot_phase_vs_kx(k_x, k_y):
    Ufull_matrix_asc_list =[]
    eigval_U_asc_T_list = []
#    eigvec_U_asc_T_list = []
    phase1_T_asc_list = []
    phase2_T_asc_list = []
    
    for i in range (len(k_x)):
        
        Ufull_matrix_asc_list.append(U_full_ascending(k_x[i], k_y, 6))
    #    print(i)
#        print(Ufull_matrix_asc_list)
        
        
    for i in range (len(k_x)):
        
        eigvals, eigvecs = np.linalg.eig(Ufull_matrix_asc_list[i])
    #    print(eigvals,'--eigen value;', eigvecs,'--eigen vector', i,'i') # all these eigen values are complex!!

        eigval_U_asc_T_list.append(eigvals)
        phase1_T_asc_list.append(- np.angle(eigvals[0]))
        phase2_T_asc_list.append(- np.angle(eigvals[1]))

# -- Does here to use a list A filled in by last i loop, it needed to have another new i loop to distribute A[i]?

    plt.plot(k_x, phase1_T_asc_list,  'o', label='phase 1(1st eigval), phi_n = phi_0')
    plt.plot(k_x, phase2_T_asc_list,  'o', label='phase 2(2nd eigval), phi_n = phi_1')

    
for i in range (len(k_y_test)):
    
    plot_phase_vs_kx(k_x_list, k_y_test[i])
    
plt.grid()
plt.show()



#%%check equation (29)

particle_hole_symmetry = []

for i in range (len(k_y_test)):
    
    
    for j in range (len(k_x_list)):
    
        eqn_29_check = U_full_ascending(k_x_list[j], k_y_test[i], 6) - np.conj(U_full_ascending(-(k_x_list[j]), -(k_y_test[i]), 6))
    
        particle_hole_symmetry.append(eqn_29_check)
        
#%% 
# plot new matrix H for grouped k_y:

# for real k_y: (set dimension N)



def real_construction_1(N):
    """
    Constructs an N x N complex matrix M where
    M[y-1, y'-1] = sum_{z=1}^{N} exp(i * 2π/N * z * (y - y'))

    Parameters:
    - N: size of the matrix and upper limit of summation

    Returns:
    - M: complex-valued NumPy array of shape (N, N)
    """
    y = np.arange(1, N + 1).reshape(-1, 1)     # shape (N, 1)
    y_prime = np.arange(1, N + 1).reshape(1, -1)  # shape (1, N)
    delta = y - y_prime                        # shape (N, N)

    z = np.arange(1, N + 1).reshape(-1, 1, 1)   # shape (N, 1, 1)
    exponent = 2j * np.pi * z * delta / N      # shape (N, N, N)

    M = np.sum(np.exp(exponent), axis=0)       # shape (N, N)
    return M
    
M = real_construction_1(4)


#%%

def real_construction____(N):
    
    M = np.zeros((N, N), dtype =complex)
    
    for i in range (N):
        
        for j in range (N):
            
            row_list = []
            
            for z in range (1, N+1):
                
                M_i_j = (1/N) * np.sum (np.exp(1j* 2 * (np.pi/N)* z* (i - j)))
                                        
            M[i, j] = M_i_j
                
    return M
                

                                        
#%% test for transform from k_spectrum matrix(diagonalised)
# (1d SSH model with equal neighbouring hopping as 1)
def real_construction_1dSSH(N):
    
    M = np.zeros((N, N), dtype =complex)
    
    for i in range (N):
        
        for j in range (N):
            
            #row_list = []
            
            M_i_j = 0
            
            for z in range (1, N+1):
                
                M_i_j += (1/N) * (np.exp(1j* 2 * (np.pi/N)* z* (i - j))) * 2 * np.cos(2* (np.pi/ N) * z)
                                        
            M[i, j] = M_i_j
                
    return M

M_1dSSH_100 = real_construction_1dSSH(100)
M_1dSSH_5 = real_construction_1dSSH(5)

#%%
#actually try this with sublattices

#for fixed kx, H_yy'
#it will be : for i in range '...' subbing k_x_list[i] into the below function to produce corresponding H_yy'
# here we firstly only try to print one of the six Hamiltonians H1,H2,...,H6

def real_Hyy_fix_k_x(N, k_x, t): # t can be any from [  t_i = T/N * i + T/(N * 2)  ]
    
    M = np.zeros((2*N, 2*N), dtype =complex)
    
    a_0 = 1 # here is the lattices length
    
    
    for i in range (N):
        
        for j in range (N): 
            
            #row_list = []
            
            M_2i_2j = 0
            M_2i_2j_1 = 0
            M_2i_1_2j = 0
            M_2i_1_2j_1 = 0
            
            for z in range (1, N+1):
                
                k_y = 2* (np.pi/ N) * z
                
                H_ky = H_matrix_2by2(k_x, k_y, t, a_0)
              #  H_ky = np.array([[1, 2], [3, 4]])
                
                M_2i_2j += (1/N) * (np.exp(1j* 2 * (np.pi/N)* z* ((i+1) - (j+1)))) * 2 * H_ky[0][0]
                
                M_2i_2j_1 += (1/N) * (np.exp(1j* 2 * (np.pi/N)* z* ((i+1) - (j+1)))) * 2 * H_ky[0][1]
                
                M_2i_1_2j += (1/N) * (np.exp(1j* 2 * (np.pi/N)* z* ((i+1) - (j+1)))) * 2 * H_ky[1][0]
                
                M_2i_1_2j_1 += (1/N) * (np.exp(1j* k_y * ((i+1) - (j+1)))) * 2 * H_ky[1][1]
                
                
# i is for row and j for column
            M[2*i, 2*j] = M_2i_2j
            
            M[2*i, 2*j+ 1 ] =  M_2i_2j_1
            
            M[2*i+ 1, 2*j ] =  M_2i_1_2j
            
            M[2*i + 1, 2*j + 1 ] =  M_2i_1_2j_1
            
            print (2*i, 2*i +1, 2*j, 2*j +1)
                
    return M

test_M_H1 = real_Hyy_fix_k_x(8, np.pi/4, (T/6) * 5+T/12)

# seems only when this J_z isn't turned off, there is the periodic term? (which is corresponding to cos(x))

# 1. Check if it works for calculating k value 
#-- diagonal this real matrix should get the same eigen values with if we don't transform,
# which means if we just subbing the value of k_y on dioganal 2x2 blocks 1-2. -- double check here

# 3. can we add up from H1 to H6 as a whole then turn off corners? or do we only do the part of one stage in 6
# 3-1. should we also do for the U operator? 










