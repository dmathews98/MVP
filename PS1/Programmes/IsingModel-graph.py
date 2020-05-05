'''
Plots the data for the magnetic susceptibility against Temperature
'''
import matplotlib.pyplot as plt
import numpy as np
import sys

def Av_E(method):
    if method == 'g':
        filename = 'data-g.txt'
    elif method == 'k':
        filename = 'data-k.txt'
    yvals = []
    xvals = []
    errors = []
    filein = open(filename, 'r')
    for line in filein.readlines():
        vals = line.split(' , ')
        xvals.append(float(vals[0]))
        yvals.append(float(vals[1]))
        errors.append(float(vals[2]))
    filein.close()

    plt.figure(1)
    plt.errorbar(xvals, yvals, yerr=errors, ecolor='r', markersize=2, capsize=3, elinewidth=0.5, barsabove=True)
    plt.title('Average energy vs Temeperature')
    plt.xlabel('Temperature')
    plt.ylabel('Average Energy')
    if method == 'g':
        plt.savefig('AvE-g.pdf')
    elif method == 'k':
        plt.savefig('AvE-k.pdf')

def Cv_N(method):
    if method == 'g':
        filename = 'data-g.txt'
    elif method == 'k':	
        filename = 'data-k.txt'
    yvals = []
    xvals = []
    errors = []
    filein = open(filename, 'r')
    for line in filein.readlines():
        vals = line.split(' , ')
        xvals.append(float(vals[0]))
        yvals.append(float(vals[3]))
        errors.append(float(vals[4]))
    filein.close()

    plt.figure(2)
    plt.errorbar(xvals, yvals, yerr=errors, ecolor='r', markersize=2, capsize=5, elinewidth=0.5, barsabove=True)
    plt.title('Heat Capacity per site vs Temeperature')
    plt.xlabel('Temperature')
    plt.ylabel('Heat capacity / N')

    if method == 'g':
        plt.savefig('CvN-g.pdf')
    elif method == 'k':
        plt.savefig('CvN-k.pdf')

def Av_M():
    filename = 'data-g.txt'
    yvals = []
    xvals = []
    errors = []
    filein = open(filename, 'r')
    for line in filein.readlines():
        vals = line.split(' , ')
        xvals.append(float(vals[0]))
        yvals.append(float(vals[5]))
        errors.append(float(vals[6]))
    filein.close()

    plt.figure(3)
    plt.errorbar(xvals, yvals, yerr=errors, ecolor='r', markersize=2, capsize=5, elinewidth=0.5, barsabove=True)
    plt.title('Average Absolute Magnetism vs Temeperature')
    plt.xlabel('Temperature')
    plt.ylabel('Average absolute Magnetism')
    plt.savefig('AvM-g.pdf')

def Chi():
    filename = 'data-g.txt'
    yvals = []
    xvals = []
    errors = []
    filein = open(filename, 'r')
    for line in filein.readlines():
        vals = line.split(' , ')
        xvals.append(float(vals[0]))
        yvals.append(float(vals[7]))
        errors.append(float(vals[8]))
    filein.close()

    plt.figure(4)
    plt.errorbar(xvals, yvals, yerr=errors, ecolor='r', markersize=2, capsize=5, elinewidth=0.5, barsabove=True)
    plt.title('Magnetic susceptibility vs Temeperature')
    plt.xlabel('Temperature')
    plt.ylabel('Magnetic susceptibility')
    plt.savefig('Chi-g.pdf')


def main():
    if len(sys.argv) != 2:
        print("Enter correct number of inputs")
        quit()
    elif sys.argv[1] not in ['g', 'k']:
        print("Enter g for Glauber or k for Kawasaki")
        quit()
    else:
        [file, method] = sys.argv

    Av_E(method)
    Cv_N(method)
    if method == 'g':
        Av_M()
        Chi()

main()
