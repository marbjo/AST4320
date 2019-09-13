import numpy as np
import matplotlib.pyplot as plt
from scipy import constants

n = 10000

c = constants.c
G = constants.G
H0 = 2e-18 #per sec
rho_0 = 3*H0**2/(8*np.pi*G)
mp = constants.m_p
k_b = constants.k
mu = 0.62 #temp

def j_l_b(z):
    #before decoupling
    value = (1+z)**(-3./2) * np.sqrt(c**2 *np.pi/(3.*G*rho_0))
    return value

def j_m_b(z):
    #before decoupling
    value = (1+z)**(-3./2) * np.pi**(5./2)/(6.*G**(3./2)*rho_0**(1./2))*(c/np.sqrt(3))**3
    return value

def j_l_a(z):
    #after decoupling
    value = 0.05*(1+z)**(-1./2)*np.sqrt(k_b*np.pi/(mu*mp*G*rho_0))
    return value

def j_m_a(z):
    #after decoupling
    value = (1+z)**(3./2) * 1.25e-4*np.pi**(5./2)/(6.*G**(3./2)*rho_0**(1./2)) * (k_b/(mu*mp))**(3./2)
    return value

z = np.linspace(2000,0,n)
j_length = np.zeros(n) #jeans_length(z)
j_mass = np.zeros(n) #jeans_mass(z)

for i in range(n):
    z_value = z[i]
    if z_value >= 1100:
        j_length[i] = j_l_b(z_value)
        j_mass[i] = j_m_b(z_value)

    if z_value < 1100:
        j_length[i] = j_l_a(z_value)
        j_mass[i] = j_m_a(z_value)

plt.figure()
plt.loglog(z,j_length)
plt.xlabel('Redshift z',fontsize=20)
plt.ylabel(r'$M_J$',fontsize=20)
plt.xlim(3e3,1e-1)

plt.title('Jeans length as a function of redshift')
plt.figure()
plt.loglog(z,j_mass)
plt.xlabel('Redshift z',fontsize=20)
plt.ylabel(r'$\lambda_J$',fontsize=20)
plt.xlim(3e3,1e-1)
plt.title('Jeans mass as a function of redshift')

plt.show()
