'''
Animation and solver for Jacobi algorithm of Poissons Eqn
'''

from __future__ import division
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy import signal

class Gauss_Seidel(object):

    def __init__(self, N, lim, phi_0, omega, type):
        self.type = type   # used for monopole/wire
        self.N = N
        self.c_lim = lim   # limit of convergence
        self.dx = 1.0

        self.phi_0 = phi_0

        self.centre_pos = [int(self.N/2.0), int(self.N/2.0), int(self.N/2.0)]   # used for access to central point

        self.noise_dampening = 100.0
        self.check = 10   # for animation, every 10 sweeps shown

        self.omega = omega   # over-relaxation parameter

        # initialise arrays
        self.lattice = np.empty((self.N, self.N, self.N))
        self.new_lattice = np.empty((self.N, self.N, self.N))
        self.rho_lattice = np.zeros((self.N, self.N, self.N))

        self.initialise_lattice(self.lattice, phi_0)
        self.initialise_rho_lattice()

    '''
    Initialises lattice with noise and init condition
    '''
    def initialise_lattice(self, lattice, phi_0):
        for i in range(self.N):
            for j in range(self.N):
                for k in range(self.N):
                    if (i == 0) or (i == self.N-1) or (j == 0) or (j == self.N-1) or (k == 0) or (k == self.N-1):   # boundaries at 0
                        lattice[i][j][k] = 0
                    else:
                        lattice[i][j][k] = np.random.randint(-10, 11)/self.noise_dampening + phi_0   # noise plus standard value to initialise elsewhere

        self.old_lattice = np.copy(self.lattice)

    '''
    Initialise rho lattice for monopole or wire problem
    '''
    def initialise_rho_lattice(self):
        if self.type == 'monopole':
            self.rho_lattice[self.centre_pos[0],self.centre_pos[1],self.centre_pos[2]] = 1
        elif self.type == 'wire':
            for k in range(self.N):
                self.rho_lattice[int(self.N/2.0)][int(self.N/2.0)][k] = 1

    '''
    Initialises lattices with noise and init condition in 2d for omega data gathering
    '''
    def initialise_lattice_2d(self):
        self.lattice = np.empty((self.N, self.N))
        self.new_lattice = np.empty((self.N, self.N))
        self.rho_lattice = np.zeros((self.N, self.N))
        for i in range(self.N):
            for j in range(self.N):
                if (i == 0) or (i == self.N-1) or (j == 0) or (j == self.N-1):   # boundaries at 0
                    self.lattice[i][j] = 0
                else:
                    self.lattice[i][j] = np.random.randint(-10, 11)/self.noise_dampening + self.phi_0   # noise plus standard value to initialise elsewhere

        self.old_lattice = np.copy(self.lattice)

        self.rho_lattice[int(self.N/2.0)][int(self.N/2.0)] = 1
############################### - LOOP METHODS - ###############################
    '''
    Method to determine when convergence reached by standard magntiude of vector difference
    Returns True/False for convergence
    '''
    def check_convergence(self):
        diff = self.lattice - self.old_lattice
        mag_diff = np.linalg.norm(diff)
        if mag_diff <= self.c_lim:
            return True
        else: return False

    '''
    Gauss-Seidel algorithm modified for SOR by loop method, for standard Gauss-Seidel use omega = 1
    '''
    def SOR_algorithm(self, i, j, k):
        if (i == 0) or (j == 0) or (k == 0):   # boundary is 0
            self.lattice[i][j][k] = 0
        else:
            self.lattice[i][j][k] = (1.0 - self.omega) * self.lattice[i][j][k] + \
                                    self.omega * (1.0/6.0) * (self.lattice[i+1][j][k] + self.lattice[i-1][j][k] +\
                                                              self.lattice[i][j+1][k] + self.lattice[i][j-1][k] +\
                                                              self.lattice[i][j][k+1] + self.lattice[i][j][k-1] +\
                                                              ((self.dx**2.0) * self.rho_lattice[i][j][k]))

    '''
    Gauss-Seidel algorithm modified for SOR by loop method in 2D for omega v convergence time plot
    '''
    def SOR_algorithm_2d(self, i, j):
        if (i == 0) or (j == 0):   # boundary is 0
            self.lattice[i][j] = 0
        else:
            self.lattice[i][j] = (1.0 - self.omega) * self.lattice[i][j] + \
                                    self.omega * (1.0/4.0) * (self.lattice[i+1][j] + self.lattice[i-1][j] +\
                                                              self.lattice[i][j+1] + self.lattice[i][j-1] +\
                                                              ((self.dx**2.0) * self.rho_lattice[i][j]))

    '''
    Iterates over lattice by given method, for SOR-GS required only to N-1
    '''
    def lattice_iterator(self, method):
        for i in range(self.N-1):
            for j in range(self.N-1):
                for k in range(self.N-1):
                    method(i, j, k)

    '''
    Iterates over lattice by given method, for SOR-GS required only to N-1
    '''
    def lattice_iterator_2d(self, method):
        for i in range(self.N-1):
            for j in range(self.N-1):
                method(i, j)

################################################################################
#                             DATA ANALYSIS METHODS                            #
################################################################################
    '''
    Distance from centre calculator
    '''
    def dist(self, i, j ,k):
        a = abs(int(self.N/2.0) - i)
        b = abs(int(self.N/2.0) - j)
        c = abs(int(self.N/2.0) - k)
        return np.sqrt(a**2.0 + b**2.0 + c**2.0)

    '''
    Magnitude of 3D vector
    '''
    def mag(self, x, y, z):
        return np.sqrt(x**2.0 + y**2.0 + z**2.0)

    '''
    Data list generator and output method for monopole charge
    '''
    def output_data_monopole(self):
        contour_data = []
        vector_data = []
        dist_data = []

        Ef = np.gradient(self.lattice)   # grad of E = potential

        for i in range(self.N):
            for j in range(self.N):
                for k in range(self.N):
                    if k == int(self.N/2.0):
                        contour_data.append([i, j, self.lattice[i, j, k]])
                        vector_data.append([i, j, -Ef[0][i][j][k], -Ef[1][i][j][k], -Ef[2][i][j][k]])
                        dist_data.append([i, j, self.dist(i,j,k), self.mag(Ef[0][i][j][k], Ef[1][i][j][k], Ef[2][i][j][k])])

        contour_data = np.array(contour_data)
        vector_data = np.array(vector_data)
        dist_data = np.array(dist_data)

        np.savetxt('monopole_contour_data'+str(self.N)+'.txt', contour_data)
        np.savetxt('monopole_vector_data'+str(self.N)+'.txt', vector_data)
        np.savetxt('monopole_dist_data'+str(self.N)+'.txt', dist_data)

    '''
    Data list generator and output method for wire charge
    '''
    def output_data_wire(self):
        contour_data = []
        vector_data = []
        dist_data = []

        # Finding curl values of 3D magnetic field
        Af = np.gradient(self.lattice)
        Bx = Af[1] - Af[2]
        By = Af[2] - Af[0]
        Bz = Af[0] - Af[1]

        for i in range(self.N):
            for j in range(self.N):
                for k in range(self.N):
                    if k == int(self.N/2.0):
                        contour_data.append([i, j, self.lattice[i,j,k]])
                        vector_data.append([i, j, Bx[i][j][k], By[i][j][k], Bz[i][j][k]])
                        dist_data.append([self.dist(i,j,k), self.lattice[i,j,k], self.mag(Bx[i][j][k], By[i][j][k], Bz[i][j][k])])

        contour_data = np.array(contour_data)
        vector_data = np.array(vector_data)
        dist_data = np.array(dist_data)

        np.savetxt('wire_contour_data'+str(self.N)+'.txt', contour_data)
        np.savetxt('wire_vector_data'+str(self.N)+'.txt', vector_data)
        np.savetxt('wire_dist_data'+str(self.N)+'.txt', dist_data)

    '''
    Runs the simulation until given convergence, prints number of sweeps and can output data to txt files. Graphs final potential
    '''
    def run_convergence(self):
        convergence_res = False   # convergence check
        nsweeps = 0

        while convergence_res == False:
            self.lattice_iterator(self.SOR_algorithm)

            convergence_res = self.check_convergence()
            nsweeps += 1
            self.old_lattice = np.copy(self.lattice)

        print('Convergence of '+str(self.c_lim)+' in '+str(nsweeps)+' sweeps')

        # --- Data Output and Imaging --- #
        # if self.type == 'monopole':
        #     self.output_data_monopole()
        # elif self.type == 'wire':
        #     self.output_data_wire()

        fig, ax = plt.subplots()
        self.mat = ax.imshow(self.lattice[:,:,self.centre_pos[2]], interpolation='gaussian', cmap='magma')
        fig.colorbar(self.mat)
        # title = str(self.type)+str(self.N)+'_potential.pdf'
        # plt.savefig(title)
        plt.show()

    '''
    Runs the simulation for range of omega values and outputs number of sweeps until desired convergence
    '''
    def gather_omega_data(self):
        w_list = np.arange(1.00, 2.00, 0.01)
        sweeps_list = np.empty(len(w_list))

        for i in range(len(w_list)):
            convergence_res = False   # convergence check
            nsweeps = 0
            self.omega = w_list[i]
            # self.initialise_lattice_2d()   # reset lattice 2d
            self.initialise_lattice(self.lattice, self.phi_0)   # reset lattice 3d

            while convergence_res == False:
                # self.lattice_iterator_2d(self.SOR_algorithm_2d)  # 2d
                self.lattice_iterator(self.SOR_algorithm)   # 3d

                convergence_res = self.check_convergence()
                nsweeps += 1
                self.old_lattice = np.copy(self.lattice)

            sweeps_list[i] = nsweeps

        w_list = np.transpose(w_list)
        sweeps_list = np.transpose(sweeps_list)

        data = np.transpose(np.vstack((w_list, sweeps_list)))

        title = 'SOR_data-3d-'+str(self.N)+'.txt'
        np.savetxt(title, data)

################################################################################
#                               ANIMATION METHODS                              #
################################################################################
    '''
    Runs simulaton and updates animation every desired number of sweeps
    '''
    def run_dynamics(self):
        # - Loop - #
        for s in range(self.check):   # only updates animation every self.check sweeps
            self.lattice_iterator(self.SOR_algorithm)

    def animate(self, i):
        self.run_dynamics()
        self.mat.set_data(self.lattice[:,:,self.centre_pos[2]])
        return self.mat,

    def run(self):
        fig, ax = plt.subplots()
        self.mat = ax.imshow(self.lattice[:,:,self.centre_pos[2]], interpolation='gaussian', cmap='magma')
        fig.colorbar(self.mat)

        anim = FuncAnimation(fig, self.animate, interval = 10, blit = True)

        plt.show()


if __name__ == "__main__":
    # I/O

    if len(sys.argv) != 7:
        print('Usage: N lim phi_0 omega(=1 standard) monopole/wire anim/conv')
        quit()
    elif sys.argv[5] not in ['monopole', 'wire']:
        print('Type should be monopole or wire')
        quit()
    elif sys.argv[6] not in ['anim', 'conv', 'gather-omega']:
        print('View should be anim or conv')
        quit()
    else:
        [file, N, lim, phi_0, omega, type, view] = sys.argv

    sim = Gauss_Seidel(int(N), float(lim), float(phi_0), float(omega), str(type))

    if view == 'anim':
        sim.run()
    elif view == 'conv':
        sim.run_convergence()
    elif view == 'gather-omega':
        sim.gather_omega_data()
