import numpy as np
import matplotlib.pyplot as plt

n = 100000

Omega_m_list = [1.0, 0.3, 0.8]
Omega_L_list = [0.0, 0.7, 0.2]
delta_list = np.zeros((n,3))
f_list = np.zeros((n,3))

H_0 = 1.
a_start = 1e-3
a_end = 1.
da = (a_end - a_start)/n
a = np.linspace(a_start, a_end, n)


def H_func(a, Omega_m, Omega_L):
    value = np.sqrt(H_0**2 *(Omega_m * a**(-3) + Omega_L) )
    return value

def dota_func(a, Omega_m, Omega_L):
    value = np.sqrt(H_0**2 *(Omega_m * a**(-1) + Omega_L*a**2) )
    return value

def dudt_func(H,delta,u):
    value = 3./2 * H**2 * delta - 2*H*u
    return value

for j in range(3):
    Omega_m = Omega_m_list[j]
    Omega_L = Omega_L_list[j]

    H = np.zeros(n); dota = np.zeros(n)
    delta = np.zeros(n); u = np.zeros(n)

    f = np.zeros(n)
    #Initial conditions
    delta[0] = a[0] #delta equal to a in early universe
    u[0] = dota_func( a[0], Omega_m, Omega_L ) #derivatives therefore also equal

    for i in range(n-1):
        H[i] = H_func(a[i], Omega_m, Omega_L)
        dota[i] = dota_func(a[i], Omega_m, Omega_L)

        u[i+1] = u[i] + dudt_func( H[i], delta[i], u[i] )*da/dota[i]
        delta[i+1] = delta[i] + u[i+1]*da/dota[i]

        f[i+1] = np.log(delta[i+1]/delta[i]) / np.log(a[i+1]/a[i])

    delta_list[:,j] = delta
    f_list[:,j] = f
    plt.loglog(a,delta_list[:,j], label=r"$\Omega_m = $%.1f, $\Omega_{\Lambda}=$%.1f" %(Omega_m, Omega_L))

plt.ylabel(r'$\delta$', fontsize=20)
plt.xlabel('a', fontsize=20)
plt.title('Time evolution of perturbed overdensity')
plt.legend()
plt.show()

plt.figure()
plt.ylabel(r'$\delta$', fontsize=20)
plt.xlabel('a', fontsize=20)
plt.title('Time evolution of perturbed overdensity')
plt.loglog(a,delta_list[:,0], label=r"$\Omega_m = $%.1f, $\Omega_{\Lambda}=$%.1f" %(Omega_m_list[0], Omega_L_list[0]))
plt.loglog(a,delta_list[:,1], label=r"$\Omega_m = $%.1f, $\Omega_{\Lambda}=$%.1f" %(Omega_m_list[1], Omega_L_list[1]))
plt.loglog(a,delta_list[:,2], label=r"$\Omega_m = $%.1f, $\Omega_{\Lambda}=$%.1f" %(Omega_m_list[2], Omega_L_list[2]))
plt.legend()
plt.axis([2*10**(-1), 1.1, 2*10**(-1), 1.1] )
plt.show()

z = 1/a - 1
plt.figure()
plt.xlabel('Redshift (z)', fontsize=20)
plt.ylabel('f', fontsize=20)
plt.loglog(z[1:],f_list[1:,0], label=r"$\Omega_m = $%.1f, $\Omega_{\Lambda}=$%.1f" %(Omega_m_list[0], Omega_L_list[0]))
plt.loglog(z[1:],f_list[1:,1], label=r"$\Omega_m = $%.1f, $\Omega_{\Lambda}=$%.1f" %(Omega_m_list[1], Omega_L_list[1]))
plt.loglog(z[1:],f_list[1:,2], label=r"$\Omega_m = $%.1f, $\Omega_{\Lambda}=$%.1f" %(Omega_m_list[2], Omega_L_list[2]))
plt.legend()
plt.show()
