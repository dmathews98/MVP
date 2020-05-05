'''
Animation class for PS1
'''

import sys
from threading import Thread
from IsingModel import Spins
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np

class SpinView(object):

    '''
    Initialises class, model lattice and figure
    '''
    def __init__(self, N, T, J, nsweeps, pfreq, Method):
        # Intialise class, lattice and figure
        self.model = Spins(N, T, J, nsweeps, pfreq, Method)
        self.model.create_lattice()
        self.Method = Method

        self.fig, self.ax = plt.subplots()
        self.mat = self.ax.imshow(self.model.lattice, cmap='winter')

        self.ani = None # For storing the animation object

    '''
    Runs thread and animation based on desired method
    '''
    def run(self, func):
        # Check for model method
        if self.Method == 'g':
            func = self.model.Glauber
        elif self.Method == 'k':
            func = self.model.Kawasaki

        thread = Thread(target=self.model.Run_dynamics, args=(func,))

        thread.start()

        self.ani = FuncAnimation(self.fig, self.animate, interval = 100, blit=True)
        plt.show()

        self.model.stop()

    '''
    Used for animation update of data - doesnt redraw artist just updates data
    '''
    def animate(self, i):
        self.mat.set_data(self.model.lattice)

        return self.mat,

if __name__ == "__main__":
    # I/O
    if len(sys.argv) != 6:
        print("Incorrect number of arguments.\nPlease input: N T Sweeps Print-frequency Method")
        quit()
    elif sys.argv[5] not in ['g', 'k']:
        print('Please enter g for Glauber or k for Kawasaki')
        quit()
    else:
        [file, N, T, nsweeps, pfreq, Method] = sys.argv

    J = 1.0   # Used if wanting to adjust units

    spin = SpinView(int(N), float(T), J, int(nsweeps), int(pfreq), Method)
    if Method == 'g':
        func = spin.model.Glauber
    elif Method == 'k':
        func = spin.model.Kawasaki
    spin.run(func)
