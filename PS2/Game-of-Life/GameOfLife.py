'''
Simulation class for Game of Life model
'''

import sys
import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import time as t

class States(object):

    def __init__(self, N, nsweeps, InitialState, runs=1):
        self.N = N
        self.nsweeps = nsweeps
        self.runs = runs
        self.istate = InitialState

        self.lattice = np.zeros((self.N,self.N))
        self.next_lattice = np.zeros((self.N,self.N))

        self.tot_alive = 0

        self.stopSim = False # Flag to stop the simulation thread

###############################################################################
#                         SIMULATION METHODS                                  #
###############################################################################

    '''
    Used for stopping the simulation when running
    '''
    def stop(self):
        self.stopSim = True

    '''
    Creates lattice based on input size and desired state
    '''
    def create_lattice(self):
        if self.istate == 'random':
            self.lattice = np.random.choice([0,1], (self.N, self.N))

        elif self.istate == 'glider':
            c = int(self.N/2)
            self.lattice[c][c] = 1
            self.lattice[c+1][c] = 1
            self.lattice[c-1][c] = 1
            self.lattice[c+1][c-1] = 1
            self.lattice[c][c-2] = 1

        elif self.istate == 'oscillator':
            c = int(self.N/2)
            self.lattice[c][c-1] = 1
            self.lattice[c][c] = 1
            self.lattice[c][c+1] = 1

        elif self.istate == 'explosion':
            c = int(self.N/2)
            self.lattice[c-1][c] = 1
            self.lattice[c-2][c] = 1
            self.lattice[c-3][c] = 1

            self.lattice[c+1][c] = 1
            self.lattice[c+2][c] = 1
            self.lattice[c+3][c] = 1

            self.lattice[c][c-1] = 1
            self.lattice[c][c-2] = 1
            self.lattice[c][c-3] = 1

            self.lattice[c][c+1] = 1
            self.lattice[c][c+2] = 1
            self.lattice[c][c+3] = 1

        elif self.istate == 'oscillator_cool':
            c = int(self.N/2)
            self.lattice[c-2][c] = 1
            self.lattice[c-3][c] = 1
            self.lattice[c-4][c] = 1

            self.lattice[c+2][c] = 1
            self.lattice[c+3][c] = 1
            self.lattice[c+4][c] = 1

            self.lattice[c][c-2] = 1
            self.lattice[c][c-3] = 1
            self.lattice[c][c-4] = 1

            self.lattice[c][c+2] = 1
            self.lattice[c][c+3] = 1
            self.lattice[c][c+4] = 1

        elif self.istate == 'pulsar':
            self.lattice = np.zeros((17, 17))
            self.lattice[2, 4:7] = 1
            self.lattice[4:7, 7] = 1
            self.lattice += self.lattice.T
            self.lattice += self.lattice[:, ::-1]
            self.lattice += self.lattice[::-1, :]

        self.next_lattice = np.copy(self.lattice)
        self.tot_alive = self.get_total_alive()

        # print(self.lattice)

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
    Returns total alive
    '''
    def get_total_alive(self):
        tot_alive = 0
        for i in range(self.N):
            for j in range(self.N):
                if self.lattice[i][j] == 1:
                    tot_alive += 1
        return tot_alive

    '''
    Counts and returns number of alive nearest neighbours
    '''
    def count_alive_neighbours(self, i, j):
        num_alive = 0
        for n in range(-1, 2):
            for m in range(-1, 2):
                if n == 0 and m == 0:   # dont include self
                    continue
                elif self.lattice[self.pbc(i+n)][self.pbc(j+m)] == 1:
                    num_alive += 1

        return num_alive

    '''
    Checks and updates state based on rules
    '''
    def update_state(self, i, j):
        number_alive = self.count_alive_neighbours(i, j)

        # for alive state: less than 2 or more than 3 alive neighbours means death at next step
        if self.lattice[i][j] == 1:
            if number_alive < 2:
                self.next_lattice[i][j] = 0
                self.tot_alive -= 1
            elif number_alive > 3:
                self.next_lattice[i][j] = 0
                self.tot_alive -= 1

        # for a dead state, 3 alive neighbours means life at next step
        elif self.lattice[i][j] == 0:
            if number_alive == 3:
                self.next_lattice[i][j] = 1
                self.tot_alive += 1

    '''
    Controls the looping over the system and updating the lattice each timestep
    '''
    def iterator(self):
        # loops over all states
        for i in range(self.N):
            for j in range(self.N):
                self.update_state(i,j)
        self.lattice = np.copy(self.next_lattice)
        # print(self.lattice)

################################################################################
#                              DATA ANALYSIS METHODS                           #
################################################################################
    '''
    Runs simulation while getting data
    '''
    def analysis_iterator(self):
        sweeps = np.arange(self.nsweeps)
        histocount = []
        for r in range(self.runs):
            self.create_lattice()
            alivecount = []
            for k in range(self.nsweeps):
                # loops over all states
                for i in range(self.N):
                    for j in range(self.N):
                        self.update_state(i,j)
                self.lattice = np.copy(self.next_lattice)

                alivecount.append(self.tot_alive)

                if (k > 5) and (alivecount[k] == alivecount[k-1]) and (alivecount[k] == alivecount[k-2]) and (alivecount[k] == alivecount[k-3]):
                    histocount.append(k)
                    break

        histocount = np.array(histocount)
        # plt.plot(sweeps, alivecount)

        #histonorm = histocount/self.runs
        np.savetxt('histodata.txt', histocount)
        #plt.figure(1)
        #plt.hist(histocount, bins=30)
        #plt.savefig('histattempt2.pdf')
        # plt.show()

    '''
    Gets the centre of mass of the glider, if at the boundary returns false
    '''
    def get_CoM(self):
        xvals = []
        yvals = []
        for i in range(self.N):
            for j in range(self.N):
                if self.lattice[i][j] == 1:
                    if (i == self.N-1) or (i == 0): #boundary checks
                        return False
                    elif (j == self.N-1) or (j == 0):
                        return False
                    else:
                        xvals.append(j)
                        yvals.append(i)

        com = [np.mean(xvals), np.mean(yvals)]

        return com

    '''
    Runs the simulation to gather the glider data
    '''
    def get_glider_speed_data(self):
        com_data = []
        sweeps_data = []
        for k in range(self.nsweeps):
            CoM = self.get_CoM()
            self.iterator()
            if CoM == False:   # check for pbs, means skip to next
                continue
            else:
                com_data.append(CoM)
                sweeps_data.append(k)


        com_data = np.transpose(np.array(com_data))
        sweeps_data = np.transpose(np.array(sweeps_data))

        data = np.transpose(np.vstack((sweeps_data, com_data)))   # stacks into [[sweep, x, y]]

        np.savetxt('glider_data.txt', data)

################################################################################
#                               ANIMATION METHODS                              #
################################################################################
    '''
    Runs simulaton for desired number of sweeps
    '''
    def run_dynamics(self):
        self.stopSim = False   # for stopping simulation when window closed
        for k in range(self.nsweeps):
            if (self.stopSim): break   # if window closes then simulation stopped
            self.iterator()

    def animate(self, i):
        self.run_dynamics()
        self.mat.set_data(self.lattice)
        return self.mat,

    def run(self):
        fig, ax = plt.subplots()
        self.mat = ax.imshow(self.lattice)

        anim = FuncAnimation(fig, self.animate, interval = 50, blit = True)

        plt.show()

if __name__ == "__main__":
    # I/O

    if len(sys.argv) != 4:
        print('Incorrect number of arguments. Please enter lattice-size number-sweeps Initial-state')
        quit()
    elif sys.argv[3] not in ['random', 'glider', 'oscillator', 'oscillator_cool', 'explosion', 'pulsar']:
        print('Please enter random, glider or oscillator for the initial state')
        quit()
    else:
        [file, N, nsweeps, Initial_State] = sys.argv

    lattice = States(int(N), int(nsweeps), str(Initial_State))
    lattice.create_lattice()
    lattice.run()

    # --- Data analyis --- #

    # Histogram

    # runs = 500
    # nsweeps = 3000
    # N = 50
    # Initial_State = 'random'
    #
    # lattice = States(int(N), int(nsweeps), str(Initial_State), runs)
    # lattice.analysis_iterator()

    # Glider

    # runs = 1
    # nsweeps = 500
    # N = 50
    # Initial_State = 'glider'
    #
    # lattice = States(int(N), int(nsweeps), str(Initial_State), runs)
    # lattice.create_lattice()
    # lattice.get_glider_speed_data()
