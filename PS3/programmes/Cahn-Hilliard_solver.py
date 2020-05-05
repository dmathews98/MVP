'''
Animation and solver for Cahn-Hilliard equation
'''

from __future__ import division
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy import signal

class CahnHilliard(object):

    def __init__(self, N, nsweeps, dx, dt, phi_0):
        self.N = N
        self.nsweeps = nsweeps
        self.M = 0.1
        self.kappa = 0.1
        self.dx = dx
        self.dt = dt
        
        self.phi_0 = phi_0

        self.noise_dampening = 100.0
        self.check = 10

        self.a = 0.1
        self.b = 0.1

        self.phi_lattice = np.empty((self.N, self.N))
        self.mu_lattice = np.empty((self.N, self.N))

        self.initialise_phi_lattice(phi_0)

    '''
    Initialises Phi lattice with noise and init condition
    '''
    def initialise_phi_lattice(self, phi_0):
        for i in range(self.N):
            for j in range(self.N):
                self.phi_lattice[i][j] = np.random.randint(-10, 11)/self.noise_dampening + phi_0

########################### - CONVULTION METHODS - #############################
    '''
    Calculates and returns laplacian of supplied lattice
    '''
    def laplacian_convolve(self, grid):
        # 2d laplacian kernel
        kernel = [[0., 1., 0.],
                  [1., -4., 1.],
                  [0., 1., 0.]]
        return signal.convolve2d(grid,kernel,boundary='wrap',mode='same')

    '''
    Calculates and returns mu based on current phi
    '''
    def mu_convolve(self):
        return -self.a*self.phi_lattice + self.b*(self.phi_lattice**3.0) - (self.kappa/(self.dx**2.0))*self.laplacian_convolve(self.phi_lattice)

    '''
    Updates phi through calculating mu
    '''
    def phi_convolve(self):
        self.phi_lattice = self.phi_lattice + (self.M * self.dt/(self.dx**2.0)) * self.laplacian_convolve(self.mu_convolve())

############################### - LOOP METHODS - ###############################
    '''
    Modulus method for image use per periodic boundary conditions
    '''
    def pbc(self, i):
        # Boolean so mod is only done on extraneous points needing image checked
        if (i > self.N-1) or (i < 0):
            image_i = i % self.N

            return int(image_i)

        else:
            return i

    '''
    Updates lattice point on mu lattice
    '''
    def update_mu_lattice_point(self, i, j):
        self.mu_lattice[i][j] = -self.a * self.phi_lattice[i][j] + self.b * (self.phi_lattice[i][j]**3.0) \
                                - (self.kappa/(self.dx**2.0)) * ( self.phi_lattice[self.pbc(i+1)][j] + self.phi_lattice[self.pbc(i-1)][j] \
                                + self.phi_lattice[i][self.pbc(j+1)] + self.phi_lattice[i][self.pbc(j-1)] - (4*self.phi_lattice[i][j]) )

    '''
    Updates lattice point on phi lattice
    '''
    def update_phi_lattice_point(self, i, j):
        self.phi_lattice[i][j] = self.phi_lattice[i][j] + (self.M * self.dt/(self.dx**2.0)) * ( self.mu_lattice[self.pbc(i+1)][j] + self.mu_lattice[self.pbc(i-1)][j] \
                                + self.mu_lattice[i][self.pbc(j+1)] + self.mu_lattice[i][self.pbc(j-1)] - (4*self.mu_lattice[i][j]) )

    '''
    Iterates over lattice by given method
    '''
    def lattice_iterator(self, method):
        for i in range(self.N):
            for j in range(self.N):
                method(i, j)

################################################################################
#                             DATA ANALYSIS METHODS                            #
################################################################################
    '''
    Calculates and returns free energy of phi lattice
    '''
    def calc_free_energy(self):
        grad = np.gradient(self.phi_lattice)
        F = -(self.a/2.0)*(self.phi_lattice**2.0) + (self.b/4.0)*(self.phi_lattice**4.0) + (self.kappa/2.0)*(grad[0]**2.0 + grad[1]**2.0)
        F_sum = 0
        for i in range(self.N):
            for j in range(self.N):
                F_sum += F[i][j]
        return F_sum * (self.dx**2.0)

    '''
    Runs over given number of sweeps and gathers free energy data which is output to a CH-f-data.txt
    '''
    def run_data_gather(self):
        data = []
        for k in range(self.nsweeps):
            self.phi_convolve()
            if (k > 100) and (k % self.check == 0):  # free energy check every set number of sweeps
                data.append([k, self.calc_free_energy()])

        title = 'CH-FE-phi'+str(self.phi_0)+'-data'+str(self.N)+'.txt'
        np.savetxt(title, data)

################################################################################
#                               ANIMATION METHODS                              #
################################################################################
    '''
    Runs simulaton for desired number of sweeps
    '''
    def run_dynamics(self):
        # - Loop - #
        # for s in range(self.check):   # only updates animation every self.check sweeps
        #   self.lattice_iterator(self.update_mu_lattice_point)
        #   self.lattice_iterator(self.update_phi_lattice_point)

        # - Convolution - #
        for s in range(self.check):   # only updates animation every self.check sweeps
            self.phi_convolve()

    def animate(self, i):
        self.run_dynamics()
        self.mat.set_data(self.phi_lattice)
        return self.mat,

    def run(self):
        fig, ax = plt.subplots()
        self.mat = ax.imshow(self.phi_lattice, interpolation='gaussian', cmap='magma')

        anim = FuncAnimation(fig, self.animate, interval = 10, blit = True)

        plt.show()

if __name__ == "__main__":
    # I/O
    
    if len(sys.argv) != 6:
        print('Usage: N nsweeps dx dt phi_0')
        quit()
    else:
        [file, N, nsweeps, dx, dt, phi_0] = sys.argv

    sim = CahnHilliard(int(N), int(nsweeps), float(dx), float(dt), float(phi_0))
    sim.run()

    # --- Data Analysis --- #

    # N = 100
    # nsweeps = 100000
    # dx = 1.0
    # dt = 1.0
    # phi_0 = 0.0
    #
    # sim = CahnHilliard(int(N), int(nsweeps), float(dx), float(dt), float(phi_0))
    # sim.run_data_gather()
    #
    # phi_0 = 0.5
    #
    # sim = CahnHilliard(int(N), int(nsweeps), float(dx), float(dt), float(phi_0))
    # sim.run_data_gather()
