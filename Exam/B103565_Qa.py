'''
B103565
MVP Exam
Contact process
'''

import sys
import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import time as t
import random as r

class ContactProcess(object):

    def __init__(self, L, p, nsweeps):
        self.L = L
        self.N = L*L
        self.p = p
        self.nsweeps = nsweeps
        self.lattice = np.empty((L,L))

        self.create_lattice(L)

    def create_lattice(self, L):
        self.lattice = np.random.choice([0, 1], size=[L, L]) # 0 inactive, 1 active
    
    def create_survival_lattice(self):
        i = r.randint(0, self.L-1)
        j = r.randint(0, self.L-1)
        self.lattice = np.zeros((self.L,self.L))
        self.lattice[i][j] = 1

###############################################################################
#                         SIMULATION METHODS                                  #
###############################################################################
    '''
    Used for stopping the simulation when running
    '''
    def stop(self):
        self.stopSim = True

    '''
    Modulus method for image use per periodic boundary conditions
    '''
    def pbc(self, i):
        # Boolean so mod is only done on extraneous points needing image checked
        if (i > self.L-1) or (i < 0):
            image_i = i % self.L

            return int(image_i)

        else:
            return i

    '''
    Checks if chosen active site should become inactive by prob 1-p
    '''
    def test_inactive_site(self):
        check = r.random()
        if check < (1-self.p):
            return True
        else:
            return False

    '''
    Checks if neighbour site should become active by prob p
    '''
    def test_inactive_neighbour(self):
        check = r.random()
        if check < (self.p):
            return True
        else:
            return False

    '''
    Randomly selects nearest neighbour and tests to be upated if inactive by prob p
    '''
    def spread_neighbour(self, i, j):
        vert = r.randint(-1,1)
        horz = r.randint(-1,1)
        if (vert == 0) and (horz == 0):
            self.spread_neighbour(i, j)  # recursively recalls function if 0,0 as same point
        else:
            if self.lattice[self.pbc(i+vert)][self.pbc(j+horz)] == 0:   # if inactive test to become active
                test = self.test_inactive_neighbour()
                if test == True:
                    self.lattice[self.pbc(i+vert)][self.pbc(j+horz)] = 1

    '''
    Checks and updates state based on rules
    '''
    def update_state(self, i, j):
        # active so becomes inactive with prob 1-p
        if self.lattice[i][j] == 1:
            test = self.test_inactive_site()
            if test == True:
                self.lattice[i][j] = 0
            
            self.spread_neighbour(i, j)   # checks random neighbour

    '''
    Controls the looping and selecting random point
    '''
    def iterator(self):
        for n in range(self.N):
            i = r.randint(0,self.L-1)
            j = r.randint(0,self.L-1)

            self.update_state(i,j)

################################################################################
#                              DATA ANALYSIS METHODS                           #
################################################################################

    def get_active_fraction(self):
        num_active = np.sum(self.lattice)
        return (num_active/float(self.N))

    def get_delta(self, vals):
        av_vals = np.mean(vals)
        av_vals_2 = np.mean(vals**2.0)
        return (av_vals_2 - av_vals**2.0)/float(self.N)

    '''
    Resampling via bootstrap method
    '''
    def resampling(self, vals):
        samplevals = []
        for i in range(len(vals)):
            j = np.random.randint(0, len(vals))
            samplevals.append(vals[j])

        return np.array(samplevals)

    '''
    Calculates variance error vis resampling
    '''
    def get_delta_err_via_resample(self, vals, resamples):
        var_list = []
        for k in range(resamples):   # calculates var for each sample point
            samp = self.resampling(vals)
            var = self.get_delta(samp) * self.N
            var_list.append(var)

        var_list = np.array(var_list)   # convert to numpy array

        # Calculating error
        var_av = np.mean(var_list)
        av_var_2 = np.mean(var_list**2.0)

        delta_err = np.sqrt(av_var_2 - var_av**2.0)/float(self.N)

        return delta_err

    def part_b(self):
        data = []
        sweeps = []
        ac_frac = []

        for i in range(self.nsweeps):
            if i > 100:   # wait 100 sweeps for equib
                self.iterator()
                current_frac = self.get_active_fraction()
                sweeps.append(i)
                ac_frac.append(current_frac)

        # data = np.transpose(np.vstack( (np.transpose(sweeps), np.transpose(ac_frac)) ))

        plt.plot(sweeps, ac_frac)
        plt.show()

    def part_c(self):
        p_range = np.arange(0.55, 0.7001, 0.005)
        av_inf = []

        for i in range(len(p_range)):
            self.p = p_range[i]
            self.create_lattice(self.L)
            ac_frac = []

            for n in range(self.nsweeps):
                if n > 100:   # wait 100 sweeps for equib
                    self.iterator()
                    current_frac = self.get_active_fraction()
                    ac_frac.append(current_frac)

            if min(ac_frac) == 0:
                av_inf.append(0)
            else:
                av_inf.append(np.mean(ac_frac))

        data = np.transpose(np.vstack( (np.transpose(p_range), np.transpose(av_inf)) ))

        np.savetxt('part_c_data.txt', data)

        plt.plot(p_range, av_inf)
        plt.show()

    def part_d(self, resamples):
        p_range = np.arange(0.55, 0.7001, 0.005)
        delta_A_list = []
        delta_A_err_list = []

        for i in range(len(p_range)):
            self.p = p_range[i]
            self.create_lattice(self.L)
            A_list = []

            for n in range(self.nsweeps):
                if n > 100:   # wait 100 sweeps for equib
                    self.iterator()
                    A = self.get_active_fraction() * self.N
                    A_list.append(A)

            delta_A_list.append(self.get_delta(np.array(A_list)))
            delta_A_err_list.append(self.get_delta_err_via_resample(np.array(A_list), resamples))

        data = np.transpose(np.vstack( (np.transpose(p_range), np.transpose(delta_A_list), np.transpose(delta_A_err_list)) ))

        np.savetxt('part_d_data.txt', data)
     
    def part_e(self, reruns):
        sweeps = np.arange(300)
        data = np.zeros(300)

        for n in range(reruns):
            self.create_survival_lattice()
            data[0] += self.get_active_fraction() * self.N
            for i in range(300):
                self.iterator()
                data[i] += self.get_active_fraction() * self.N

        data = data/(float(reruns)*self.N)

        plt.plot(sweeps, data)
        plt.show()

    def part_f(self, reruns):
        p_range = np.array([0.6, 0.625, 0.65])
        sweeps = np.arange(300)
        data = np.zeros((3,300))
        
        for p in range(len(p_range)):
            self.p = p_range[p]
            for n in range(reruns):
                self.create_survival_lattice()
                data[p][0] += self.get_active_fraction() * self.N
                for i in range(300):
                    self.iterator()
                    data[p][i] += self.get_active_fraction() * self.N

        data = data/(float(reruns)*self.N)

        data = np.transpose(np.vstack((sweeps, data)))

        np.savetxt('part_f_data.txt', data)


################################################################################
#                               ANIMATION METHODS                              #
################################################################################

    def run_dynamics(self):
        self.stopSim = False   # for stopping simulation when window closed
        for k in range(10):
            if (self.stopSim): break   # if window closes then simulation stopped
            self.iterator()

    def animate(self, i):
        self.iterator()
        self.mat.set_data(self.lattice)
        return self.mat,

    def run(self):
        fig, ax = plt.subplots()
        self.mat = ax.imshow(self.lattice, cmap='seismic')
        fig.colorbar(self.mat)

        anim = FuncAnimation(fig, self.animate, interval = 100)

        plt.show()

if __name__ == "__main__":
    # I/O

    if len(sys.argv) != 5:
        print('Incorrect number of arguments. Please enter L p nsweeps part')
        quit()
    else:
        [file, L, p, nsweeps, part] = sys.argv

    simulation = ContactProcess(int(L), float(p), int(nsweeps))

    if part == 'anim':
        simulation.run()
    elif part == 'b':
        simulation.part_b()
    elif part == 'c':
        simulation.part_c()
    elif part == 'd':
        resamples = 10
        simulation.part_d(resamples)
    elif part == 'e':
        reruns = 20
        simulation.part_e(reruns)
    elif part == 'f':
        reruns = 20
        simulation.part_f(reruns)