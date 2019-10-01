import numpy as np
import random as random
import matplotlib.pyplot as plt

from matplotlib import rc

import scipy.integrate as integrate

class RandomWalk():

    def __init__(self, eps, var_init):
        #Initializing variables
        self.eps = eps
        self.Sc_init = (np.pi/var_init**2)**(1./4)

    def variance(self):
        #Calculating variance
        value = np.pi / self.Sc**4
        return value

    def PDF(self, delta, Sc):
        #Analytic expression for the prob. distribution
        var = np.pi / Sc**4
        value1 = 1./(np.sqrt(2*np.pi*var)) * np.exp(- delta**2 / (2*var))
        return value1

    def PDF_NC(self, delta, Sc):
        #Analytic expression for the prob. distribution where no cross
        var = np.pi / Sc**4

        delta_crit = 1.0

        value2 =  1/(np.sqrt(2*np.pi*var)) * ( np.exp(-delta**2/(2*var)) - np.exp(-(2*delta_crit-delta)**2/(2*var)))
        return value2

    """
    def Walk(self):

        #One random walk

        #Resetting the variables
        self.Sc = self.Sc_init
        self.delta = 0

        while self.Sc > 1:

            #delta_crit_list = []
            #Sc_crit_list = []

            var = self.variance()

            #New S_c reduced by epsilon
            self.Sc = self.Sc - self.eps

            #New variance is difference between old variances
            var_diff = self.variance() - var

            #Updating density
            std_dev = np.sqrt(var_diff) #Standard deviation
            self.delta = self.delta + np.random.normal(scale = std_dev)

            #Only if it doesn't cross
            #if self.delta < 1:
            #    delta_crit_list.append(self.delta)
            #    Sc_crit_list.append(self.Sc)

        return self.delta , self.Sc
        """
    def WalkCrit(self):

        #One random walk, which never passes delta_crit = 1

        #Resetting the variables
        self.Sc = self.Sc_init
        #self.delta = 0

        var = self.variance()
        self.delta = np.random.normal(scale=np.sqrt(var))

        check = 0

        while self.Sc > 1:

            var = self.variance()

            #New S_c reduced by epsilon
            self.Sc = self.Sc - self.eps

            #New variance is difference between old variances
            var_diff = self.variance() - var

            #Updating density
            std_dev = np.sqrt(var_diff) #Standard deviation
            self.delta = self.delta + np.random.normal(scale = std_dev)

            #Only if it doesn't cross
            if self.delta >= 1:
                check = 1

        return self.delta , self.Sc, check

    def MultipleWalk(self,N):
        #Function to do multiple walks

        Sc_arr = np.zeros(N)
        delta_arr = np.zeros(N)
        delta_crit = []
        Sc_crit = []

        for i in range(N):
            d, S_c, check = self.WalkCrit()

            delta_arr[i] = d
            Sc_arr[i] = S_c

            if check==0:
                delta_crit.append(d)
                Sc_crit.append(S_c)

        return delta_arr, Sc_arr, np.array(delta_crit), np.array(Sc_crit)

    def Distribution(self, N):
        #Plotting prob. distribution

        #Performing the walks
        delta_arr , Sc_arr , delta_crit_arr, Sc_crit_arr = self.MultipleWalk(N)

        #Plotting vector and analytic values
        delta_vec = np.linspace(np.min(delta_arr), np.max(delta_arr), len(delta_arr))
        analytic = self.PDF(delta_vec, Sc_arr)

        plt.hist(delta_arr, bins=50, histtype='bar', ec='black', color = 'green', label='Random Walk', density = True)
        plt.plot(delta_vec, analytic, color = "black", linewidth = 1, label = "Analytic PDF")

        plt.xlabel(r'$\delta$', fontsize=20)
        plt.ylabel(r'P($\delta | $ M)', fontsize=20)
        plt.title('Gaussian random walk versus analytical expression', fontsize=20)
        plt.legend()
        plt.show()

    def NC_Distribution(self, N):
        #Plotting prob. distribution

        #Performing the walks
        delta_arr , Sc_arr , delta_crit_arr, Sc_crit_arr = self.MultipleWalk(N)

        #Plotting vector and analytic values
        delta_crit_vec = np.linspace(np.min(delta_crit_arr), np.max(delta_crit_arr), len(delta_crit_arr))
        analytic2 = self.PDF_NC(delta_crit_vec, Sc_crit_arr)

        analytic2norm = analytic2 / np.trapz(analytic2,delta_crit_vec)

        plt.hist(delta_crit_arr, bins=50, histtype='bar', ec='black', color = 'green', label='Random Walk', density = True)
        plt.plot(delta_crit_vec, analytic2norm, color = "black", linewidth = 1, label = "Analytic PDF") #ONLY FITS WITH RANDOM SCALING OF ~2.5 ?????


        plt.xlabel(r'$\delta$', fontsize=20)
        plt.ylabel(r'$P_{nc}(\delta | $ M)', fontsize=20)
        plt.title('Gaussian random walk versus analytical expression', fontsize=20)
        plt.legend()
        plt.show()


if __name__ == '__main__':
    initial_var = 0.5e-4 #Half of 10^-4
    eps = 1e-2
    N = int(1e4)
    RW = RandomWalk(eps,initial_var)

    """
    #Test if analytic function is normalized???
    def PDF_test(delta):
        #Analytic expression for the prob. distribution
        var = np.pi / 1**4
        value1 = 1./(np.sqrt(2*np.pi*var)) * np.exp(- delta**2 / (2*var))
        return value1

    I = integrate.quad(PDF_test, -6,1)
    print(I)
    """
    output = RW.Distribution(N)
    #output = RW.NC_Distribution(N)
