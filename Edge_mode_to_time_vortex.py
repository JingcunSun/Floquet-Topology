

# few general summary of the codes for my CV
# all these codes are working with the Hamiltonian as well as eigen vectors, eigen values of it.
# also include FT (real space, k space)



import numpy as np
from functools import reduce
from scipy.linalg import expm
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from collections import Counter
from collections import Counter
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

    plt.plot(k_x, phase1_T_asc_list,  '.', color='tab:blue', markersize = 3)
    plt.plot(k_x, phase2_T_asc_list,  '.', color='tab:purple', markersize = 3)

    
for i in range (len(k_y_list)):
    
    plot_phase_vs_kx(k_x_list, k_y_list[i])
    
plt.ylim (-np.pi, np.pi)    
    
plt.xlabel(r'$k_x a_0\ \mathrm{(PBC)}$', fontsize=13)
plt.ylabel(r'Phase $\phi_n  = \epsilon_n T$  (T = 1)', fontsize=13)

plt.title(r'Bulk phase band plot vs $k_x$ (PBC in $x,y$)', fontsize=13)


from matplotlib.lines import Line2D # This is purely for plot legend 

legend_elements = [
    Line2D([0], [0], marker='.', color='tab:blue',
           linestyle='None', markersize=6,
           label='Phase band 1 (1st eigenvalue)'),
    Line2D([0], [0], marker='.', color='tab:purple',
           linestyle='None', markersize=6,
           label='Phase band 2 (2nd eigenvalue)')
]

plt.legend(handles=legend_elements, loc='right', bbox_to_anchor=(1, 0.56))
plt.grid()
plt.show()


#%%
# The loop(PLOT) below is: 
#try to not plot 2pi, for that I feel the phase1/phase2 separate phase band plot seems having a bit of strange shift -
# - which the same type of that when I plot only specific k_y value(also, for k_y doesn't equal to 0, 
# - phase1 or phase2's plot doesn't need to obey particle-hole symmetry)

# THIS IS JUST A TEST WHICH DOES NOTHING
#for i in range (len(k_y_test_list_no_including_2pi)):
    
#    plot_phase_vs_kx(k_x_list, k_y_test_list_no_including_2pi[i])
    
plt.xlabel(r'$k_x\ \mathrm{(PBC)}$', fontsize=12)
plt.ylabel(r'Phase $\phi_n$  (T = 1)', fontsize=12)

plt.title(r'Bulk phase band plot vs $k_x$ (PBC in $x,y$), NO y = 2pi**', fontsize=13)

plt.legend()
plt.grid()
plt.close()
# here if change plt.close to plt.show can get the plot

# well, alright - there's almost no difference here.


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


#%% HERE WE WRITE A TEXT FUNCTION FOR PLOT PHASE BAND (equasienergy spectrum) OF real basis PBC in xy H
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
    
    plt.plot(kx_list_for_plot, phi_n_list_of_kx_list, '.', markersize=3)
#%%    
#    print (len(kx_list_for_plot))
#    print(len(phi_n_list_of_kx_list))

#CHIRAL EDGE STATE PLOT
Ny_check = 10 # this is the Ly = N*a0 before Ly defined 
              # if choose same number as len(k_y_list) --> num of y here = initial N_y generate ky list 
              # this replot (plot_check) of phase band will be same as the initial plotfunc [plot_phase_vs_kx]

test_phi_along_kx = plot_check_phase_vs_kx(k_x_list, Ny_check)
#plt.xlim(1.7, 1.8)
plt.ylim(-np.pi, np.pi)
#plt.ylim(-1.75, -0.8)
plt.xlabel('$k_x a_0$')
plt.ylabel('Phase at time T ($\epsilon(k, T)T$)')
plt.grid()
plt.title('Bulk band double check for real basis H2Nx2N, N = 10')
#plt.legend()
plt.tight_layout()
plt.show()
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
        
        print(top_right, x)

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
Ly = 10 # here is the strip in real space y direction of our model -- kind of useless in this section
 
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
    plt.plot(kx_list_for_plot, phase_list_of_list, '.', markersize=3)
            
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
    
plot_edge_vs_kx(k_x_list, Ny_plot_yOBC_edge)
#plt.xlim(1.7, 1.8)
plt.ylim(-np.pi, np.pi)
#plt.ylim(-1.75, -0.8)
plt.xlabel('$k_x a_0$')
plt.ylabel('Phase at time T ($\epsilon(k, T)T$)')
plt.grid()
plt.title(rf'Bulk and Edge Phases vs $k_x a_0$, y direction unit cell number $Ly  = {Ny_plot_yOBC_edge}$')
#plt.legend()
plt.tight_layout()
plt.show()

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

print('Phi n (vec) begin', Phi_n_vec_list, 'vec')
 
# here is a QUESTION!!: why Phi_n_vec_list is for 1 in var explorer is 20x20 (cofusing to read)
print(Phi_n_vec_list[0], Phi_n_vec_list[1]) 

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
        print(comp_wise_modulus_sqr_1, 'COM SQR 1')
        comp_sum_for_1 = np.sum(comp_wise_modulus_sqr_1)
        print(comp_sum_for_1, 'COM SQR 1 SUM')
        
        Rho_specific_y = np.sum(np.abs(Phi_vec_n1)) + np.sum(np.abs(Phi_vec_n2)) 
            # this corresponding to modulus square for each Phi_n pair in sublattices, which is Rho_n
        print(Rho_specific_y, 'rho for each of 10 y')
        Rho_y_listed_by_y. append(Rho_specific_y)
        
#        for k in range(len(vec_y_alpha1)): #(actually this is just j)
#actually the above is useless -- let's try this n label below, for each n corresponding to one Phi_n, state_n
    
    Phi_eigvecs_n_list = []
    
    
    Phi_eigvecs_n_Rho_edge_list = []
    
    Phi_eigvecs_n_list_Rho_comp_y_list =[]
    
    for n in range(2*Ny_plot_yOBC_edge): # I guess here is to spam both vec_n1 and vec_n2? ******Why here is 40???
        
        Phi_eigvecs_n = Phi_specific_kx [:, n]
        
        Phi_eigvecs_n_list. append(Phi_eigvecs_n)
        
        
        #here try to find component y', alpha1 and y' alpha2 (y and y' are the same thing)
        Rho_comp_Phi_eigvec_specific_n_list = []
        
#        N = 20 #(only for record: here we have Ly from 1 to 10, and this number of real y basis is N(showed previously))
        
#        Ly = N
        
        for y in range (Ny_plot_yOBC_edge): # y is in range 10, because y is from 1 to Ly = 10, given previously (N =10)
            
            Rho_component_n = (np. abs(Phi_eigvecs_n[2*y]))**2 + (np.abs(Phi_eigvecs_n[2*y + 1]))**2
            
                       
            # here the Rho correponding to each eigvec n' component's modulus sqaure. But np.abs gives the modulus
            
            Rho_comp_Phi_eigvec_specific_n_list.append(Rho_component_n)
        
        print('n = ', n, 'density of state(Rho) for each y = ', Rho_comp_Phi_eigvec_specific_n_list)
              # for check 1
              # or go to the var explorer, see Phi_eigvecs_n_list_Rho_each_y_list, this list is for the dos for each y and - next line 
              #each n, should be len(n[len(Ly)])
              
 #       try change the def of 'edge'
            
        Rho_edge_Phi_eigvec_n = Rho_comp_Phi_eigvec_specific_n_list[0]+ Rho_comp_Phi_eigvec_specific_n_list[Ny_plot_yOBC_edge -1] #+ Rho_comp_Phi_eigvec_specific_n_list[2] + Rho_comp_Phi_eigvec_specific_n_list[Ly -3] + Rho_comp_Phi_eigvec_specific_n_list[Ly -2] + Rho_comp_Phi_eigvec_specific_n_list[Ly -1]
        #total density of edege state
        #total density of edege state
        print('n = ', n, 'edge density(Rho) for all y =', Rho_edge_Phi_eigvec_n)
        # for check 2 
        
        Phi_eigvecs_n_Rho_edge_list.append(Rho_edge_Phi_eigvec_n )#* 100) # here we can try enlarging edge density rho(useless)
        
        Phi_eigvecs_n_list_Rho_comp_y_list.append(Rho_comp_Phi_eigvec_specific_n_list) # this is just for check, we don't need that much
        
        
    # here I want to check two things: 1. for each Phi n, among desities labelled by y, is there any edge showed
    #                                     higher probability?
    #                                  Conclusion: Not really??? So what does it mean by localised? (Maybe can ask professor)
    
    #                                  2. doubel check if the Rho_edge is really the sum of the two of Rho_comp, 0 and 10(as index 9)?
    #                                  Conclusion: checks out 
    
            
        
    Phi_specific_kx_20_eigvecs = []
    
    for k in range(2*Ny_plot_yOBC_edge): # here to check 20 eigvecs
        
        Phi_specific_kx_20_eigvecs.append(Phi_specific_kx[k])
        
    print ('i=', Phi_specific_kx_20_eigvecs, '20 eigen vec') #here we have 20 eigen vec, each of it has 20 components
    
    print(' check which order is y and alpha')
        
                                   
   # print('i =', i, Rho_y_listed_by_y, 'rho list for all 10 y, should be len = 10') 
    
    Rho_edge_specific_kx = Rho_y_listed_by_y[0] + Rho_y_listed_by_y[Ny_plot_yOBC_edge -1] #+ Rho_y_listed_by_y[Ly -2] + Rho_y_listed_by_y[Ly -3] + Rho_y_listed_by_y[Ly -4] + Rho_y_listed_by_y[Ly -10]# here 19 is Ly -1
    # Q: why increase this can actually increas e the number of Rho_edge_listed_by_sta_n_listed_by_kx
    
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
    
print("len(k_x_list) =", len(k_x_list))
print("len(Rho_edge_listed_by_sta_n_listed_by_kx) =",
      len(Rho_edge_listed_by_sta_n_listed_by_kx))


#%% here in our plot of phase band structure with intensity (labelled by the density of edge)

# which need to be corrected 
    
plt.figure()

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
        
        
    plt.scatter([kx_in_p_loop]*len(rho_edge_ns_specific_kx), epsilon_ns_specific_kx, c=rho_edge_ns_specific_kx, cmap='inferno_r', s=1, vmin = 0, vmax = 1)

    #plt.scatter([kx_in_p_loop]*len(rho_edge_ns_specific_kx), epsilon_ns_specific_kx, c=rho_edge_ns_specific_kx, cmap='magma_r', s=1, vmin=vmin, vmax=vmax)

plt.colorbar()
#plt.xlim(-0.4, 0.17)      # x-axis range if want to see full picture please include 0, coz there is a group of epsi_n at 0
#plt.ylim(-2.2, -1)     # y-axis range
plt.xlabel(r'$k_x a_0$')
plt.ylabel(r'$\epsilon_\alpha(k_x) T$')
#plt.legend()
plt.title(rf"Phase bands with edge state density diagnosis, y OBC, $L_y = {Ny_plot_yOBC_edge}$")
plt.grid()
plt.show()


#%% real basis x, y PBC, check spectrum
# Plan A: [H(2x2)] to x,y real basis [PBC(2NyNx 2NyNx)]; check spectrum: compare with Nx&Ny generated kx&ky list as vars-all epsilon

# Plan B: [PBC in y H(2Ny x 2Ny)] to x,y real basis [PBC(2NyNx x 2NyNx)]; 
#         check spectrum: compare with Nx generated kx list, same Ny 
#                         - this 'check Plan B spectrum' already compared with 'check Plan A spectrum'


# spectrum need to be find by H(2NxNy x 2NxNy) append 6 piecewisely(will be changed to time discretisation calculate U(T) - epsilon
# AIM: to make the specrtum as the same

# try plan A: 

def H2NyNx_PBC(N_y, N_x, t, a_0):
    
    
    H = np.zeros((2*N_y*N_x, 2*N_y*N_x), dtype =complex)
        # clarify an idea: there's a correspondence between k_i and N_i here: k_i = 2* (np.pi/ (N_i*a_0)) * z 
#    dummy_k_y_list = generate_k_y_list(0, 2 * np.pi, N_y) # double check if z need to start somewhere? 0? 2pi? 

#    dummy_k_x_list = generate_k_x_list(0, 2 * np.pi, N_x)
    
    
    #double sum of all kx and ky, total summed terms: N_x * N_y
    kx_test = []
    for z_1 in range (1, N_x+1):
        k_x = 2 * np.pi * z_1/ N_x 
        kx_test.append(k_x)
        
        for z_2 in range (1, N_y+1):
            k_y = 2 * np.pi * z_2/ N_y
# here k_x and k_y is summed over allN_x, N_y, cover rounds and rounds on H entrices 

            for i in range(N_x* N_y):
                for j in range( N_x* N_y):
                    x = (1+i )% N_x+1
                    x_1 = (1+j )% N_x+1
                    y = (1+i) // N_x + 1
                    y_1 = (1+j) // N_x + 1
                    H[2*i:2*i+2, 2*j:2*j+2] += (1/N_x) * np.exp(1j*  k_x * ( x - x_1)) * H_matrix_2by2(k_x, k_y, t, a_0) * (1/N_y) * np.exp(1j*  k_y * (y - y_1 ))
                    
                    
  #                  print(i,j, 'ij')
#    print (kx_test)
    return H

#test_H = H2NyNx_PBC(3, 4, 1/12, a_0)
#%%
#H_ky = H_matrix_2by2(k_x, k_y, t, a_0)
#M[2*i:2*i+2, 2*j:2*j+2] += (1/N_y) * (np.exp(1j*  k_y * ((i+1) - (j+1)))) * H_ky
#print(test_H,'H2NyNxtest')
#%%
def H_6_2NyNx_PBC(N_y, N_x, a_0): #here this function should have same spectrum(EIGEN_yVALS)as the [K BASIS H2x2]
                    #6 refers as list of 6 piecewise H 
    
    H2NyNx_PBC_6_matrices = []
    
    for t_x in range (6): #x refers index for t
    
        t_piecewise = T/6 *  t_x + T/(6 * 2)   # here this x is equivalent to i in the 2x2 H list func
        
        H2NyNx_PBC_piecewise =  H2NyNx_PBC(N_y, N_x, t_piecewise, a_0)
        
        H2NyNx_PBC_6_matrices.append(H2NyNx_PBC_piecewise) # generate 6 2N_y real H
        
    return H2NyNx_PBC_6_matrices
        
def Ufull_descending_PBC_realxy(N_y, N_x): # this is for whole period T
# var index: 1        
    U_list_1 = []
#reference :
#    H_list_6_2N_y_edge = H_6_2N_y_OBC(N, k_x)
    
    for i in range (6):

        U_i = expm(-1j* (H_6_2NyNx_PBC(N_y, N_x, a_0)[i])* (T/6))
        
        U_list_1.append(U_i) # append should be in the loop scope, 
                            # one grid inside the for loop title if you want to append the thing produced once loop (per i)
    
    U_reverse_1 = np.flip(U_list_1)
    
    U_full_1 = reduce(np.matmul,  U_reverse_1)
        
    return U_full_1

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
def find_spectrum_NxNy(N_y, N_x):
    eigvals, eigvecs = np.linalg.eig(Ufull_descending_PBC_realxy(N_y, N_x))
    epsilon = - np.angle(eigvals)
    return  epsilon

#%%
N_y_TEST = 10
N_x_TEST = 15
#k_x_TEST = [1 * 2* np.pi/4, 2  *2*  np.pi/4, 3 *2*  np.pi/4, 2*np.pi]
k_x_TEST = generate_k_x_list(0, 2*np.pi, N_x_TEST)
spectrum_test_realxy = find_spectrum_NxNy(N_y_TEST, N_x_TEST)
spectrum_test_realy = find_spectrum_realy(N_y_TEST, k_x_TEST)
#print('spectrum_test_realxy:', spectrum_test_realxy, 'spectrum_test_realy', spectrum_test_realy)

import numpy as np

def same_values_ndp_yes_no(a, b, ndp):
    a_rounded = np.round(a, ndp)
    b_rounded = np.round(b, ndp)
    return "yes" if np.allclose(sorted(a_rounded), sorted(b_rounded)) else "no"

if_sepctrum_are_the_same = same_values_ndp_yes_no(spectrum_test_realxy, spectrum_test_realy, 5)
print('if_sepctrum_are_the_same', if_sepctrum_are_the_same)

#%% learn how to build dict:
    # LIST of allowed 2-vec s
    # these 2-vec are different from what we are doing in our project: kx, ky list in our porject are list/array of scalars(1 num)
    # H is matrix and eigenvecs are: arrays 
points = [(0.0, 0.0), (1.0, 0.0), (0.5, 0.866)]
# Build dictionary: (kx, ky) → index
A = {} # **********-- what actually build dictionary is {}
for i, (apple, pear) in enumerate(points):
    A[(round(apple, 2), round(pear, 2))] = i

#print(A,'print out dictionary')

# end of the method -- not useful here, a 1d array is more suitable to be used in a 1d list of vectors. We don't need '{}'(dict)
# for a matrix! It is already indexed 
        
#%% TEST:
sum_test1 = 0
for zx in range (1,  4):
    kx = zx
for zy in range (1,  5):
    ky = zy
    
    sum_test1 += kx * ky
    print (kx, ky) # here kx is always 3
    print(sum_test1)
#print(sum_test1)
# so when I tap the same amount in front of both of them: the first runs, finished; then second. and sum is inside the second
#%%
sum_test = 0
for zx in range (1,  4):
    kx = zx
    for zy in range (1,  5):
        ky = zy
    
    sum_test += kx * ky
    print(sum_test)
#print(sum_test)
# for this one, the sum_test runs after zy is fixed(loop finished) as 4
#%% double sum correct template
sum_test2 = 0
for zx in range (1,  4):
    kx = zx
    for zy in range (1,  5):
        ky = zy
    
        sum_test2 += kx * ky
        print(sum_test2)
#print(sum_test2)












