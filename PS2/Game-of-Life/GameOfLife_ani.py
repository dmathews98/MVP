'''
Animation class for PS1
'''

import sys
from threading import Thread
from GameOfLife import States
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np

class StateView(object):

    '''
    Initialises class, model lattice and figure
    '''
    def __init__(self, N, nsweeps, InitialState):
        # Intialise class, lattice and figure
        self.model = States(N, nsweeps, InitialState)
        self.model.create_lattice()

        self.fig, self.ax = plt.subplots()
        self.implot = self.ax.imshow(self.model.lattice, cmap='winter')

        self.ani = None # For storing the animation object

    '''
    Runs thread and animation based on desired method
    '''
    def run(self):
        thread = Thread(target=self.model.run_dynamics, args=())

        thread.start()

        self.ani = FuncAnimation(self.fig, self.animate, interval = 100, blit=True)
        plt.show()

        self.model.stop()

    '''
    Used for animation update of data - doesnt redraw artist just updates data
    '''
    def animate(self, i):
        self.implot.set_data(self.model.lattice)

        return self.implot,

if __name__ == "__main__":
    # I/O
    if len(sys.argv) != 4:
        print("Incorrect number of arguments.\nPlease input: N Initial-state")
        quit()
    elif sys.argv[3] not in ['random', 'glider']:
        print('Please enter random or glider for Initial-state')
        quit()
    else:
        [file, N, nsweeps, Initial_State] = sys.argv

    gameoflife = StateView(int(N), int(nsweeps), str(Initial_State))
    gameoflife.run()
