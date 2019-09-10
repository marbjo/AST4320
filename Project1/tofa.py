import numpy as np
import matplotlib.pyplot as plt

n = 100000
a_start = 10**(-4)
a_end = 1
a = np.linspace(a_start, a_end, n)

T_gas = 2.5*10**(-3)/a**2

T_photon = 2.7/a

plt.loglog(a,T_gas, label=r'$T_{gas}$')
plt.loglog(a,T_photon, label=r'$T_{\lambda}$')
a_dec = 9.082652134*10**(-4)
plt.loglog([a_dec, a_dec], [0, 3000], 'k--', label=r'$a_{decoupling}$')
plt.title('Temperature evolution of gas and photons in time', fontsize= 15)
plt.xlabel('Temperature (K)', fontsize=30)
plt.ylabel('Scale factor a', fontsize= 30)
plt.legend()
plt.show()
