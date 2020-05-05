'''
Histogram plotter
'''

import sys
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    # I/O

    if len(sys.argv) != 3:
        print('Usage python plot-hist.py filename bins')
        quit()
    else:
        [init, filename, bins] = sys.argv

    data = np.loadtxt(str(filename))

    plt.hist(data, bins=int(bins))
    plt.show()
