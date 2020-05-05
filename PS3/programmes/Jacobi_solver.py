'''
Animation and solver for Poisson's Equation using the Jacobi algorithm by both
loops methods and convolution for a single monopole charge and a single wire
'''

from __future__ import division
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy import signal

class Jacobi(object):

    def __init__(self, N, lim, phi_0, type):
        self.type = type   # used for monopole/wire
        self.N = N
        self.c_lim = lim   # limit of convergence
        self.dx = 1.0

        self.centre_pos = [int(self.N/2.0), int(self.N/2.0), int(self.N/2.0)]   # used for access to central point

        self.noise_dampening = 100.0
        self.check = 10   # for animation, every 10 sweeps shown

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

    '''
    Initialise rho lattice for monopole or wire problem
    '''
    def initialise_rho_lattice(self):
        if self.type == 'monopole':
            self.rho_lattice[self.centre_pos[0],self.centre_pos[1],self.centre_pos[2]] = 1
        elif self.type == 'wire':
            for k in range(self.N):
                self.rho_lattice[int(self.N/2.0)][int(self.N/2.0)][k] = 1


########################### - CONVULTION METHODS - #############################
    '''
    Calculates and returns laplacian of supplied lattice
    '''
    def laplacian_convolve_3d(self, grid):
        # 3D laplacian kernel
        kernel = [[[0., 0., 0.],
                   [0., 1., 0.],
                   [0., 0., 0.]],

                  [[0., 1., 0.],
                   [1., 0., 1.],
                   [0., 1., 0.]],

                  [[0., 0., 0.],
                   [0., 1., 0.],
                   [0., 0., 0.]]]

        return signal.fftconvolve(grid,kernel, mode='same')

    '''
    Updates phi through convolution
    '''
    def phi_convolve_3d(self):
        self.new_lattice = (1.0/6.0)*(self.laplacian_convolve_3d(self.lattice) + (self.dx**2.0)*self.rho_lattice)

############################### - LOOP METHODS - ###############################
    '''
    Method to determine when convergence reached by standard magntiude of vector difference
    Returns True/False for convergence
    '''
    def check_convergence(self):
        diff = self.new_lattice - self.lattice
        mag_diff = np.linalg.norm(diff)
        if mag_diff <= self.c_lim:
            return True
        else: return False

    '''
    Jacobi algorithm by loop method, included for option and completeness but substantially slower
    Achieves convergence in same steps as convolution
    '''
    def jacobi_algorithm(self, i, j, k):
        if (i == 0) or (i == self.N-1) or (j == 0) or (j == self.N-1) or (k == 0) or (k == self.N-1):   # keep boundary at 0
            self.new_lattice[i][j][k] = 0
        else:
            self.new_lattice[i][j][k] = (1.0/6.0) * (self.lattice[i+1][j][k] + self.lattice[i-1][j][k] +\
                                                 self.lattice[i][j+1][k] + self.lattice[i][j-1][k] +\
                                                 self.lattice[i][j][k+1] + self.lattice[i][j][k-1] + (self.dx**2.0)*self.rho_lattice[i][j][k])

    '''
    Iterates over lattice by given method
    '''
    def lattice_iterator_jacobi(self, method):
        for i in range(self.N):
            for j in range(self.N):
                for k in range(self.N):
                    method(i, j, k)

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
                        dist_data.append([self.dist(i,j,k), self.lattice[i,j,k], self.mag(Ef[0][i][j][k], Ef[1][i][j][k], Ef[2][i][j][k])])
                        contour_data.append([i, j, k, self.lattice[i, j, k]])
                        vector_data.append([i, j, k, -Ef[0][i][j][k], -Ef[1][i][j][k], -Ef[2][i][j][k]])

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
                        dist_data.append([self.dist(i,j,k), self.lattice[i,j,k], self.mag(Bx[i][j][k], By[i][j][k], Bz[i][j][k])])
                        contour_data.append([i, j, k, self.lattice[i,j,k]])
                        vector_data.append([i, j, k, Bx[i][j][k], By[i][j][k], Bz[i][j][k]])

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
            # - Loop - #
            # self.lattice_iterator_jacobi(self.jacobi_algorithm)

            # - Convolution - #
            self.phi_convolve_3d()

            convergence_res = self.check_convergence()
            nsweeps += 1
            self.lattice = np.copy(self.new_lattice)

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

################################################################################
#                               ANIMATION METHODS                              #
################################################################################
    '''
    Runs simulaton and updates animation every desired number of sweeps
    '''
    def run_dynamics(self):
        # - Loop - #
        # for s in range(self.check):   # only updates animation every self.check sweeps
        #     self.lattice_iterator_jacobi(self.jacobi_algorithm)
        #     self.lattice = np.copy(self.new_lattice)

        # - Convolution - #
        for s in range(self.check):   # only updates animation every self.check sweeps
            self.phi_convolve_3d()
            self.lattice = np.copy(self.new_lattice)

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

    if len(sys.argv) != 6:
        print('Usage: N c_lim phi_0 monopole/wire anim/conv')
        quit()
    elif sys.argv[4] not in ['monopole', 'wire']:
        print('Type should be monopole or wire')
        quit()
    elif sys.argv[5] not in ['anim', 'conv']:
        print('Method should be anim or conv')
        quit()
    else:
        [file, N, lim, phi_0, type, method] = sys.argv

    sim = Jacobi(int(N), float(lim), float(phi_0), str(type))

    if method == 'anim':
        sim.run()
    elif method == 'conv':
        sim.run_convergence()
