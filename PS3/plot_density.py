'''
Density plotter for Q3.3
'''

import numpy as np
import matplotlib.pyplot as plt

def carrier_conc_intrinsic(T, mn, mp, eg):
    n = 2.0 * ( (kb * T/(2.0 * np.pi * hbar**(2.0)))**(3.0/2.0) ) * ( (mn * mp)**(3.0/4.0) ) * np.exp( -eg / (2.0*kb*T) )
    return n

def main():
    global kb
    global hbar

    kb = 8.617e-5   # eV K^-1
    hbar = 6.582e-16   # eV s^-1

    me = 0.51e6  # eV c^-2
    ev = 1.6e-18   # convert eV to J
    # Si data
    mn = 1.08 * me
    mp = 0.81 * me
    eg = 1.12

    nlist = []
    Tlist = np.linspace(100, 600, 10000)


    for T in Tlist:
        n_t = carrier_conc_intrinsic(T, mn, mp, eg)
        nlist.append(n_t)

    nlist = np.array(nlist)
    Tlist = np.array(Tlist)

    # List manipulation
    # Tlist = 1000 / Tlist

    plt.plot(Tlist, nlist)
    plt.xlabel(r'1000/T ($K^{-1}$)')
    plt.ylabel(r'Carrier concentration (cm$^3$)')
    plt.yscale('log')
    plt.show()

main()
