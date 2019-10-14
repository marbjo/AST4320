import numpy as np
import matplotlib.pyplot as plt

def W_k(k,R):
    eps = 1e-16
    value = np.where(abs(k) < eps, 2*R, 2*np.sin(k*R)/k)
    return value

def FWHM(W, k):
    #Only for symmetric functions
    half_max = np.max(W)/2
    kmax_idx = np.argmin( np.abs(W - half_max) )
    value = k[kmax_idx]
    width = abs(2*value)
    return value, width, half_max

n = 1000 #iterations, accuracy
k_vec = np.linspace(-20,20,n)
R = 1.

plt.figure()
R_range = 1 #Number of plots
for R in range(1,R_range+1):
    plt.subplot(R_range,1,R)
    W_vec = W_k(k_vec, R)
    plt.plot(k_vec,W_vec, label='R=%s' %R)

    plt.xlabel('k', size= 20)
    plt.ylabel(r'$\widetilde{W}\left(k\right)$', size=20)
    plt.grid(linestyle='--')

    FWHM_value = FWHM(W_vec,k_vec)
    print('Width is: %s' %FWHM_value[1])
    plt.plot([-FWHM_value[0], FWHM_value[0]], [FWHM_value[2], FWHM_value[2] ], 'k', label='FWHM = %.3f' %FWHM_value[1], linewidth=2)
    plt.legend(fontsize=16)

plt.suptitle('Fourier transform of the Top-Hat smoothing function', size=20)
plt.show()
