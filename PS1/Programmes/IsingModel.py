'''
Ising and Kawasaki model for lattice spins
'''

import sys
import time as t
import math as m
import numpy as np

'''
Class for spin lattice and manipulation
'''
class Spins(object):

    def __init__(self, N, T, J, nsweeps, pfreq, Method):
        self.J = J   # used if want proper units
        self.k = 1.0   # used if want proper units
        self.N = N
        self.T = T
        self.nsweeps = nsweeps
        self.pfreq = pfreq
        self.Method = Method

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
    Creates random initial lattice based on size input
    '''
    def create_lattice(self):
        latt = np.zeros((self.N,self.N))
        for i in range(self.N):
            for j in range(self.N):
                test = np.random.random()  # two equal probability outcomes
                if test < 0.5:
                    latt[i][j] = -1.0
                elif test >= 0.5:
                    latt[i][j] = 1.0

        self.lattice = latt

    '''
    Modulus method for image use per periodic boundary conditions
    '''
    def image_check(self, i):
        # Boolean so mod is only done on extraneous points needing image checked
        if (i > int(self.N-1)) or (i < 0):
            image_i = i%self.N

            return int(image_i)

        else:
            return i

    '''
    Calculates nearest neighbour energy
    '''
    def nn_energy(self, i, j):
        Enn = 0
        Enn += self.lattice[self.image_check(i+1)][j]   # right
        Enn += self.lattice[self.image_check(i-1)][j]   # left
        Enn += self.lattice[i][self.image_check(j+1)]   # above
        Enn += self.lattice[i][self.image_check(j-1)]   # below
        Enn = 2 * self.J * Enn * self.lattice[i][j]

        return Enn

    '''
    Metropolis test used to determine whether to change state or not
    - T used as input for exponential
    '''
    def Metropolis_test(self, E_diff):
        if E_diff <= 0:
            return True

        elif E_diff > 0:
            p = min(1, m.exp(-1.0*(E_diff)/(self.k*self.T)))
            r = np.random.random()
            if r <= p:
                return True
            elif r > p:
                return False

    '''
    Glauber dynamics method of choosing random single point
    '''
    def Glauber(self):
        # Select random point
        i = np.random.randint(0, self.N)
        j = np.random.randint(0, self.N)

        E_diff = self.nn_energy(i,j)   # Calculate energy diff

        outcome = self.Metropolis_test(E_diff)   # Metropolis test to check if flip is accepted

        if outcome == True:
            self.lattice[i][j] *= -1   # flip accepted

    '''
    Kawasaki dynamics method of choosing two random points and swapping sequentially
    '''
    def Kawasaki(self):
        # Select two random points
        i1 = np.random.randint(0, self.N)
        j1 = np.random.randint(0, self.N)

        i2 = np.random.randint(0, self.N)
        j2 = np.random.randint(0, self.N)

        # If both spins are the same, swapping will not change anything so ignore
        if self.lattice[i1][j1] != self.lattice[i2][j2]:

            # Calculate nearest neighbour energy diff from first flip
            E_diff = self.nn_energy(i1,j1)

            # Calculate nearest neighbour energy diff of second flip
            E_diff += self.nn_energy(i2,j2)

            # If same column nearest neighbour then remove double count
            if (i1 == i2) and ((j1 == self.image_check(j2+1)) or (j1 == self.image_check(j2-1))):
                E_diff += 4 * self.J

            # If same row nearest neighbour then remove double count
            elif (j1 == j2) and ((i1 == self.image_check(i2+1)) or (i1 == self.image_check(i2-1))):
                E_diff += 4 * self.J

            outcome = self.Metropolis_test(E_diff)   # Metropolis test to check if swap is accepted

            # If outcome accepted then flip
            if outcome == True:

                #!NOTE! Should I just multiply each by -1 ???

                prev1 = self.lattice[i1][j1]   # store
                self.lattice[i1][j1] = self.lattice[i2][j2]   # swap states
                self.lattice[i2][j2] = prev1

    '''
    Runs the simulation, used for animation only
    '''
    def Run_dynamics(self, func):
        self.stopSim = False   # for stopping simulation when window closed

        for i in range(self.nsweeps):   # runs over desired sweeps

            if (self.stopSim): break   # if window closes then simulation stopped

            for n in range(self.N*self.N):   #runs over all sites per sweep
                func()   # runs the desired update method for site

###############################################################################
#                         DATA GATHERING AND ANALYSIS                         #
###############################################################################

    '''
    Calculates total and average energy
    '''
    def total_E(self):
        totE = 0
        for i in range(self.N):
            for j in range(self.N):
                # Only check right and above of all so no double count
                eij = self.lattice[self.image_check(i+1)][j]   # right
                eij += self.lattice[i][self.image_check(j+1)]   # above
                totE += -1.0 * self.J * eij * self.lattice[i][j]

        self.totE = totE

    '''
    Calculates total absolute magnetisation
    '''
    def get_tot_M(self):
        self.totM = float(abs(np.sum(self.lattice)))

    '''
    Calculates Heat Cap / N
    '''
    def get_heat_cap_N(self, Evals, Eav):
        # np.var  maybe ???
        sumE2 = np.sum(Evals**2.0)
        av_E2 = float(sumE2)/len(Evals)
        Cv_N = (av_E2 - Eav**2.0)/(self.N**2.0 * self.k * self.T**2.0)

        return Cv_N

    '''
    Calculates magnetic susceptibility
    '''
    def get_mag_sus(self, Mvals, Mav):
        # np.var  maybe ???
        sumM2 = np.sum(Mvals**2.0)
        av_M2 = float(sumM2)/len(Mvals)
        m_sus = (av_M2 - Mav**2.0)/(self.N**2.0 * self.k * self.T)

        return m_sus

    '''
    Get average energy error
    '''
    def get_E_err(self, Evals, Eav):
        sumE2 = np.sum(Evals**2.0)
        av_E2 = float(sumE2)/len(Evals)
        E_err = np.sqrt(((av_E2 - (Eav)**2.0)/(len(Evals))))

        return E_err

    '''
    Get average absolute magnetism error
    '''
    def get_M_err(self, Mvals, Mav):
        sumM2 = np.sum(Mvals**2.0)
        av_M2 = float(sumM2)/len(Mvals)
        M_err = np.sqrt(((av_M2 - (Mav)**2.0)/(len(Mvals))))

        return M_err

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
    Calculates heat capacity error
    '''
    def get_Cv_err(self, Evals, resamples):
        c = []
        for k in range(resamples):   # calculates heat capacity for each sample point
            e_samp = self.resampling(Evals)
            e_samp = np.array(e_samp)
            e_av = np.mean(e_samp)
            cv = self.get_heat_cap_N(e_samp, e_av)
            c.append(cv)

        c = np.array(c)   # convert to numpy array for faster calculation

        # Calculating error

        c_av = np.mean(c)

        sumc2 = np.sum(c**2.0)
        av_c2 = float(sumc2)/len(c)

        c_err = np.sqrt(av_c2 - c_av**2.0)

        return c_err

    '''
    Calculates magnetic susceptibility error
    '''
    def get_chi_err(self, Mvals, resamples):
        chi = []
        for k in range(resamples):   # calculates chi for each sample point
            m_samp = self.resampling(Mvals)
            m_samp = np.array(m_samp)
            m_av = np.mean(m_samp)
            magsus = self.get_mag_sus(m_samp, m_av)
            chi.append(magsus)

        chi = np.array(chi)   # convert to numpy array for faster calculation

        # Calculating error

        chi_av = np.mean(chi)

        sumchi2 = np.sum(chi**2.0)
        av_chi2 = float(sumchi2)/len(chi)

        chi_err = np.sqrt(av_chi2 - chi_av**2.0)

        return chi_err

    '''
    Runs the simulation over a temperature range in reverse and gathers the data to be output for graphing
    - Reverse allows equilibration period to still be <200 sweeps at high temperature and so effectively like cooling a system
    '''
    def gather_data(self, func, temp, resamples):
        Temperature = []
        Average_Energy = []
        AvEnergy_Error = []
        Heat_Capacity = []
        HeatCap_Error = []
        if self.Method == 'g':   # only initialise if Glauber as Kawasaki doesnt use
            Average_AbsMagnetism = []
            AvMag_Error = []
            Mag_Susceptibility = []
            MagSus_Error = []

        for t in range(len(temp)-1, -1, -1):   # Working backwards (high to low) temp so can use random initialisation
            self.T = temp[t]   # set temperature for simulation at working temperature
            Temperature.append(self.T)   # create new array; used to ensure correct temperature for data in file

            Evals = []   # initialise array for energy values
            Mvals = []   # initialise array for magnetisation values

            for i in range(self.nsweeps):   # runs over desired sweeps

                for n in range(self.N*self.N):   #runs over all sites per sweep
                    func()   # runs the desired update method for site

                if i > 200:   # if sweep past 100 sweeps; relaxation period
                    if i % self.pfreq == 0:   # if at desired printing sweep
                        self.total_E()
                        Evals.append(self.totE)   # add energy to array
                        if self.Method == 'g':   # if glauber also add absolute magnetism
                            self.get_tot_M()
                            Mvals.append(abs(self.totM))

            Evals = np.array(Evals)   # numpy faster calculation
            Mvals = np.array(Mvals)

            # Getting data and errors

            Eav = np.mean(Evals)
            Average_Energy.append(Eav)   # average energy output to file
            AvEnergy_Error.append(self.get_E_err(Evals, Eav))

            Heat_Capacity.append(self.get_heat_cap_N(Evals, Eav))   # heat capity per N output to file
            HeatCap_Error.append(self.get_Cv_err(Evals, resamples))

            if self.Method == 'g':   # only get magnetisation things for Glauber
                Mav = np.mean(Mvals)
                Average_AbsMagnetism.append(Mav)   # average absolute magnetism to file
                AvMag_Error.append(self.get_M_err(Mvals, Mav))

                Mag_Susceptibility.append(self.get_mag_sus(Mvals, Mav))   # magnetic susceptibility to file
                MagSus_Error.append(self.get_chi_err(Mvals, resamples))

        # Print to files
        if self.Method == 'k':
            self.write_to_kawasaki_files(Temperature, Average_Energy, AvEnergy_Error, Heat_Capacity, HeatCap_Error)
        elif self.Method == 'g':
            self.write_to_glauber_files(Temperature, Average_Energy, AvEnergy_Error, Heat_Capacity, HeatCap_Error, Average_AbsMagnetism, AvMag_Error, Mag_Susceptibility, MagSus_Error)

    '''
    Outputs kawasaki data to files for the graphing programme
    '''
    def write_to_kawasaki_files(self, temp, Eav, Eerr, Cv, cverr):
        f = open('data-k.txt', 'a')
        for i in range(len(Eav)):
            f.write(str(temp[i]) + ' , '
            + str(Eav[i]) + ' , ' + str(Eerr[i]) + ' , '
            + str(Cv[i]) + ' , ' + str(cverr[i]) + '\n')
        f.close()

    '''
    Outputs glauber data to files for the graphing programme
    '''
    def write_to_glauber_files(self, temp, Eav, Eerr, Cv, cverr, Mav, Merr, Mchi, chi_err):
        f = open('data-g.txt', 'a')
        for i in range(len(Eav)):
            f.write(str(temp[i]) + ' , '
            + str(Eav[i]) + ' , ' + str(Eerr[i]) + ' , '
            + str(Cv[i]) + ' , ' + str(cverr[i]) + ' , '
            + str(Mav[i]) + ' , ' + str(Merr[i]) + ' , '
            + str(Mchi[i]) + ' , ' + str(chi_err[i]) + '\n')
        f.close()

if __name__ == "__main__":
    # I/O
    # if len(sys.argv) != 6:
    #     print("Incorrect number of arguments.\nPlease input: N T Sweeps Print-frequency Method")
    #     quit()
    # elif sys.argv[5] not in ['g', 'k']:
    #     print('Please enter g for Glauber or k for Kawasaki')
    #     quit()
    # else:
    #     [file, N, T, nsweeps, pfreq, Method] = sys.argv
    #
    # J = 1.0
    #
    # spin = Spins(int(N), float(T), float(J), int(nsweeps), int(pfreq), Method)
    # spin.create_lattice()
    # if Method == 'g':
    #     func = spin.Glauber
    # elif Method == 'k':
    #     func = spin.Kawasaki
    # spin.Run_dynamics(func)

    # For data gathering for plotting

    J = 1.0   # Allows actual dimensions to be used if desired
    N = 50
    T = 1.0   # for initisation
    T_range = np.arange(1, 3.1, 0.1)
    pfreq = 10
    nsweeps = 1000
    resamples = 5

    # Glauber
    Method = 'g'
    spin = Spins(int(N), float(T), J, int(nsweeps), int(pfreq), Method)
    spin.create_lattice()
    func = spin.Glauber
    spin.gather_data(func, T_range, resamples)

    # Kawasaki
    Method = 'k'
    spin = Spins(int(N), float(T), J, int(nsweeps), int(pfreq), Method)
    spin.create_lattice()
    func = spin.Kawasaki
    spin.gather_data(func, T_range, resamples)
