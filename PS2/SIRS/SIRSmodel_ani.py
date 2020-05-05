'''
Animation class for PS1
'''

import sys
from threading import Thread
from SIRSmodel import Coronavirus
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np

class SIRSView(object):

    '''
    Initialises class, model lattice and figure
    '''
    def __init__(self, N, p1, p2, p3, nsweeps):
        # Intialise class, lattice and figure
        self.model = Coronavirus(N, p1, p2, p3, nsweeps)
        self.model.create_lattice()

        self.fig, self.ax = plt.subplots()
        self.mat = self.ax.imshow(self.model.wuhan, cmap='seismic')
        self.fig.colorbar(self.mat)

        self.ani = None # For storing the animation object

    '''
    Runs thread and animation based on desired method
    '''
    def run(self):
        thread = Thread(target=self.model.run_dynamics, args=())

        thread.start()

        self.ani = FuncAnimation(self.fig, self.animate, interval = 100)
        plt.show()

        self.model.stop()

    '''
    Used for animation update of data - doesnt redraw artist just updates data
    '''
    def animate(self, i):
        self.mat.set_data(self.model.wuhan)

        return self.mat,

if __name__ == "__main__":
    # I/O
    if len(sys.argv) != 6:
        print('Incorrect number of arguments. Please enter lattice-size p1 p2 p3 number-sweeps')
        quit()
    else:
        [file, N, p1, p2, p3, nsweeps] = sys.argv

    pandemic = SIRSView(int(N), float(p1), float(p2), float(p3), int(nsweeps))
    pandemic.run()
