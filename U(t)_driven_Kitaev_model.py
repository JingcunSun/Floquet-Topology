# Code written to compute the total Unitary evolution operator for the periodically driven Kitaev honeycomb model in 2D 
# Code produces plots for the bulk states 
import numpy as np
from scipy.linalg import expm
import matplotlib.pyplot as plt

def kitaev_hamiltonian(kx, ky, Jx, Jy, Jz, a0):
    dx_k = Jx * np.sin(kx * a0) + Jz * np.sin(ky * a0)
    dy_k = Jx + Jy * np.cos(kx * a0) + Jz * np.cos(ky * a0)
    sig_x = np.array([[0, 1], [1, 0]], dtype=complex)
    sig_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    H_k = dx_k * sig_x + dy_k * sig_y
    return H_k



a0 = 1.0
T = 2*np.pi
t_step = T / 6 # For the pulse protocall, this is NOT delta_T
J0 =  0.45   # Trust me this works   
kx = np.linspace(0, 2* np.pi, 100)
ky = np.linspace(0, 2* np.pi, 100)

# 6 pulses per period 
pulse = [
    (J0, 0, J0),    # T/6 Jx and Jz 
    (J0, 0, 0),    # 2*T/6 Jx 
    (J0, J0, 0),    # 3*T/6 Jx and Jy
    (0, J0, 0),    # 4*T/6 Jy 
    (0, J0, J0),    # 5*T/6 Jy and Jz 
    (0, 0, J0),    #  T Jz 
]
eigvals = np.zeros((len(kx), len(ky), 2), dtype=complex)
unitary_ops = []
unitary_ops_neg = []
for i, kx_val in enumerate(kx):
    for j, ky_val in enumerate(ky):
        U_total = np.eye(2, dtype=complex)
        for (Jx, Jy, Jz) in pulse:
            U = expm(-1j * kitaev_hamiltonian(kx_val, ky_val, Jx, Jy, Jz, a0) * t_step)
            U_total = U @ U_total
        eigvals[i, j, :] = np.linalg.eigvals(U_total)
epsilon_T = -np.angle(eigvals)  # as eigenvalues are exp(-i*epsilon*T)
epsilon_T_kyslice = epsilon_T[:, 0, :]  # Eigenvalues for fixed ky, ky = 0 

# Full 2D plot of bands
for i in range(len(ky)):
    epsilon_ky_vals = epsilon_T[:, i, :]
    bands_sorted1 = np.sort(epsilon_ky_vals, axis=1)
    plt.plot(kx, bands_sorted1[:, 0], color='blue')
    plt.plot(kx, bands_sorted1[:, 1], color='orange')
plt.title('Floquet Quasi-energy Bands for Periodically Driven Kitaev Model (All ky)')
plt.xlabel(r'$k_x a_0$')    
plt.ylabel(r'$\varepsilon T$')
plt.xlim([0, 2*np.pi])
#plt.ylim([-np.pi, np.pi])
plt.legend()
plt.savefig('Floquet_bands_kitaev_model_all_ky_test.png')
plt.show()
