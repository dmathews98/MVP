'''
Class to model and simulate
'''

import sys
import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import random as r

class SIRS(object):

    def __init__(self, N, p1, p2, p3, nsweeps, immunity):
        self.N = N
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.nsweeps = nsweeps
        self.wuhan = np.zeros((N,N))

        self.f_i = immunity

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
    Creates lattice based on input size
    '''
    def create_lattice(self):
        # for i in range(self.N):
        #     for j in range(self.N):
        #         test = np.random.random()  # two equal probability outcomes
        #         if test > 0.66:
        #             self.wuhan[i][j] = 1   # recovered
        #         elif test < 0.33:
        #             self.wuhan[i][j] = -1   # infected
        #         else:
        #             self.wuhan[i][j] = 0   # susceptible
        self.lattice = np.random.choice([-1, 0, 1, 2], size=[self.N, self.N], p=[(1-self.f_i)/3.0, (1-self.f_i)/3.0, (1-self.f_i)/3.0, self.f_i])
        # sets with equal probabilities for 3 states and the set probability for immune state

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
    Counts and returns True if there is a nearby infected state
    '''
    def check_neighbours(self, i, j):
        # for n in range(-1, 2):
        #     for m in range(-1, 2):
        #         if n == 0 and m == 0:   # dont include self
        #             continue
        #         elif n == -1 or n == 1 and m != 0:
        #             continue
        #         elif m == -1 or m == 1 and n != 0:
        #             continue
        #         elif self.wuhan[self.pbc(i+n)][self.pbc(j+m)] == -1:
        #             return True
        #         else:
        #             return False
        if self.lattice[self.pbc(i+1)][j] == -1:   # down
            return True
        elif self.lattice[self.pbc(i-1)][j] == -1:  # up
            return True
        elif self.lattice[i][self.pbc(j+1)] == -1:   # right
            return True
        elif self.lattice[i][self.pbc(j-1)] == -1:   # left
            return True
        else:
            return False

    '''
    Checks and updates state based on rules
    '''
    def update_state(self, i, j):
        # for susceptible state: if near infected has prob p1 to become infected
        if self.lattice[i][j] == 0:
            near_i = self.check_neighbours(i, j)
            if near_i == True:
                test = r.random()
                if test <= self.p1:
                    self.lattice[i][j] = -1

        # for an infected state, p2 prob to become recovered
        elif self.lattice[i][j] == -1:
            test = r.random()
            if test <= self.p2:
                self.lattice[i][j] = 1

        # for a recovered state, p3 prob to become susceptible
        elif self.lattice[i][j] == 1:
            test = r.random()
            if test <= self.p3:
                self.lattice[i][j] = 0

    '''
    Controls the looping and selecting random point
    '''
    def iterator(self):
        for n in range(self.N*self.N):
            i = r.randint(0,self.N-1)
            j = r.randint(0,self.N-1)

            self.update_state(i,j)

################################################################################
#                              DATA ANALYSIS METHODS                           #
################################################################################

    '''
    Returns average number of infected
    '''
    def count_infected(self):
        num_infected = 0
        for i in range(self.N):
            for j in range(self.N):
                if self.lattice[i][j] == -1:
                    num_infected += 1.0

        return num_infected

    '''
    Resampling via bootstrap method
    '''
    def resampling(self, vals):
        samplevals = []
        for i in range(len(vals)):
            j = np.random.randint(0, len(vals))
            samplevals.append(vals[j])

        return samplevals

    '''
    Calculates variance error vis resampling
    '''
    def get_var_err_via_resample(self, vals, resamples):
        var = []
        for k in range(resamples):   # calculates var for each sample point
            samp = self.resampling(vals)
            samp = np.array(samp)
            av = np.mean(samp)
            variance = np.var(samp)
            var.append(variance)

        var = np.array(var)   # convert to numpy array for faster calculation

        # Calculating error

        var_av = np.mean(var)
        av_var2 = float(np.sum(var**2.0))/len(var)

        var_err = np.sqrt(av_var2 - var_av**2.0)/(self.N**2.0)

        return var_err

    '''
    Gathers the data for the no immune states plots
    '''
    def run_no_immune_data_gather(self, resamples):
        p1_range = np.arange(0, 1.01, 0.025)   # prob ranges for varying p1 and p3
        # p3_range = np.arange(0, 1.01, 0.025)
        p3_range = np.array([0.5])   #for cut in plane as need 10k sweeps
        p_data = []   # list initialised for ease of data output to file in gnuplot form

        av_I_data = np.array([])
        var_I_data = np.array([])

        resample_data = []

        # av_I_data = np.empty((len(p1_range), len(p3_range)))   # initialises arrays for average Infected states
        # var_I_data = np.empty((len(p1_range), len(p3_range)))   # and variance

        for p1_i in range(len(p1_range)):   # over p1 values
            self.p1 = p1_range[p1_i]
            for p3_j in range(len(p3_range)):   # over p3 values
                self.p3 = p3_range[p3_j]

                p_data.append([p1_range[p1_i], p3_range[p3_j]])
                self.create_lattice()   # new lattice needed each time
                I_count = []   # will stored infected count at each sweep
                for k in range(self.nsweeps):
                    # ignore first 100 sweeps to ensure reached steady state
                    if k > 100:
                        I = self.count_infected()   # count infected
                        # if num infected reaches 0 then its an absorbing state so set list to 0 for average infected to be 0
                        if I == 0:
                            I_count = [0]
                            break   # break out of loop as have reached absorbing state
                        I_count.append(I)

                    self.iterator()   # update state for next sweep

                I_count = np.array(I_count)
                I_av = np.mean(I_count)/(self.N**2.0)   # find average infected per site for current p1,p3
                I_var = np.var(I_count)/(self.N**2.0)  # variance per site for current p1,p3

                av_I_data = np.append(av_I_data, I_av)
                var_I_data = np.append(var_I_data, I_var)

                # when p3 = 0.5 also carry out resampling for variance error in the slice
                if p3_range[p3_j] == 0.5:
                    I_var_err = self.get_var_err_via_resample(I_count, resamples)
                    resample_data.append([p1_range[p1_i], p3_range[p3_j], I_av, I_var, I_var_err])

                # av_I_data[p1_i, p3_j] = I_av   # set p1,p3 average infected per site for contour
                # var_I_data[p1_i, p3_j] = I_var   # set p1,p3 varaiance per site for contour

        p_data = np.array(p_data)
        p_data = np.transpose(p_data)

        av_I_data = np.transpose(av_I_data)
        var_I_data = np.transpose(var_I_data)

        data = np.transpose(np.vstack((p_data, av_I_data, var_I_data)))
        resample_data = np.array(resample_data)

        #np.savetxt('SIRS_no_immune_cont_data.txt', data)
        np.savetxt('SIRS_no_immune_slice_data.txt', resample_data)

        # plt.pcolormesh(p1_range, p3_range, av_I_data)
        # plt.savefig('i_sirs.png')
        # plt.pcolormesh(p1_range, p3_range, var_I_data)
        # plt.savefig('var_sirs.png')

    '''
    Gathers the data for the immune plots
    '''
    def run_immune_data_gather(self, runs):
        Immune_range = np.arange(0, 0.51, 0.01)
        data = []

        # loop over immune fraction range
        for i in range(len(Immune_range)):
            self.f_i = Immune_range[i]
            mean_I_vals_at_fi = []   # used for getting error on average infected over number of runs

            for r in range(runs):   # number of runs for errors
                self.create_lattice()   # new lattice needed each time
                I_count = []   # will stored infected count at each sweep

                for k in range(self.nsweeps):
                    # ignore first 100 sweeps to ensure reached steady state
                    if k > 100:
                        I = self.count_infected()   # count infected

                        # if num infected reaches 0 then its an absorbing state so set list to 0 for average infected to be 0
                        if I == 0:
                            I_count = [0]
                            break   # break out of loop as have reached absorbing state

                        I_count.append(I)

                    self.iterator()   # update state for next sweep

                I_count = np.array(I_count)
                I_av_run = np.mean(I_count)  # average for run

                mean_I_vals_at_fi.append(I_av_run)

            I_av_of_runs = np.mean(mean_I_vals_at_fi)/(self.N**2.0)
            I_err_of_runs = np.std(mean_I_vals_at_fi)/(self.N**(2.0))   # standard error on mean is standard deviation/root N, /N^2 for per site -> N^5/2

            data.append([self.f_i, I_av_of_runs, I_err_of_runs])  # append [f_i, I_av, I_err]

        data = np.array(data)

        np.savetxt('SIRS_immune_data.txt', data)


################################################################################
#                               ANIMATION METHODS                              #
################################################################################

    def run_dynamics(self):
        self.stopSim = False   # for stopping simulation when window closed
        for k in range(self.nsweeps):
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

    if len(sys.argv) != 7:
        print('Incorrect number of arguments. Please enter lattice-size p1 p2 p3 number-sweeps immunity')
        quit()
    else:
        [file, N, p1, p2, p3, nsweeps, immunity] = sys.argv

    pandemic = SIRS(int(N), float(p1), float(p2), float(p3), int(nsweeps), float(immunity))
    pandemic.create_lattice()
    pandemic.run()

    # --- Data Gathering --- #

    # No immuninity

    # N = 50
    # p1 = 0
    # p2 = 0.5
    # p3 = 0
    # nsweeps = 10000
    # resamples = 30
    #
    # immunity = 0.0
    #
    # pandemic = SIRS(int(N), float(p1), float(p2), float(p3), int(nsweeps), float(immunity))
    # pandemic.run_no_immune_data_gather(resamples)

    # Immunity

    # N = 50
    # p1 = 0.5
    # p2 = 0.5
    # p3 = 0.5
    # nsweeps = 1000
    # runs = 5
    #
    # immunity = 0.0
    #
    # pandemic = SIRS(int(N), float(p1), float(p2), float(p3), int(nsweeps), float(immunity))
    # pandemic.run_immune_data_gather(runs)
