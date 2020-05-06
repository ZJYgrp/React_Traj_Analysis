# This script should be revised to an object-oriented program, which defines TS extraction class, trajectories classfication class, etc.
import numpy as np
import sys
import os
import glob
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import LineCollection
from matplotlib.collections import PatchCollection
from matplotlib.colors import colorConverter
from matplotlib.ticker import FuncFormatter
from scipy import stats
from scipy.stats import norm


# get atom index first
def get_atomindex():
    atom = [int(x) for x in
            input("Input four 4 atomic indexs corresponding to 2 bonds, using ' ' as delimiter.\n").split()]
    print('atom 1 ', atom[0], 'atom 2 ', atom[1], 'atom 3 ', atom[2], 'atom 4 ', atom[3])
    judge = input("Do you think the four indexes are reasonable ?(y/n)")
    if judge != 'y':
        print('I have to quit, sorry!')
        exit(1)
    print
    'run ...'
    return atom


# Creating new directory
def mkdir():
    os.system('rm -r ./trajTS')
    os.system('rm -r ./reorder')
    os.system('mkdir ./trajTS')
    os.system('mkdir ./reorder')
    return 1


# all the following files are used for extracting 1. TS structure, 2. Time-distance-distance(TDD) and 3. reordered traj to be from R to P
def open_file(filename):
    filein = open(filename, 'r')
    lines = filein.readlines()
    n_lines = len(lines)
    n_atoms = int(lines[0].split()[0])
    n_idx = int((n_lines / (n_atoms + 2)))
    fileout_TS_xyz = open('./trajTS/trajTs.xyz', 'a')
    fileout_TS = open('./trajTS/trajTs.txt', 'a')
    fileout_reorder = open('./reorder/' + filename.split('/')[2], 'w')
    return lines, n_lines, n_atoms, n_idx, filein, fileout_TS_xyz, fileout_TS, fileout_reorder


def looking_for_ts(lines, n_lines, n_atoms, n_idx):
    runpoint = lines[1].split()[6]
    if int(runpoint) != 1:
        print('The traj file has been reordered, break!')
        return '1', '1'
    for i in range(1, n_idx):
        if runpoint == lines[1 + i * (n_atoms + 2)].split()[6]: break
    n1 = i
    n2 = n_idx - i
    return n1, n2


def get_distance(lines, n, n_atoms, atom):
    x = [0, 0, 0, 0]
    y = [0, 0, 0, 0]
    z = [0, 0, 0, 0]
    for i in range(0, 4): x[i] = float(lines[(n_atoms + 2) * n + atom[i] + 1].split()[1])
    for i in range(0, 4): y[i] = float(lines[(n_atoms + 2) * n + atom[i] + 1].split()[2])
    for i in range(0, 4): z[i] = float(lines[(n_atoms + 2) * n + atom[i] + 1].split()[3])
    bond1 = round(((x[1] - x[0]) ** 2 + (y[1] - y[0]) ** 2 + (z[1] - z[0]) ** 2) ** .5, 3)
    bond2 = round(((x[3] - x[2]) ** 2 + (y[3] - y[2]) ** 2 + (z[3] - z[2]) ** 2) ** .5, 3)
    return bond1, bond2


def determine_order(n1, n2, lines, n_atoms, atom, fileout_TS):
    bond_ts1, bond_ts2 = get_distance(lines, 0, n_atoms, atom)
    bond_p1, bond_p2 = get_distance(lines, n1 - 1, n_atoms, atom)
    bond_r1, bond_r2 = get_distance(lines, n1 + n2 - 1, n_atoms, atom)
    print('From reactant to product, the bond 1 changes from R:', bond_r1, ' to TS:', bond_ts1, ' then to P:', bond_p1)
    print('From reactant to product, the bond 2 changes from R:', bond_r2, ' to TS:', bond_ts2, ' then to P:', bond_p2)
    print('Presumbly bond 1 forms from R to P')
    if (bond_r1 > bond_p1): return 'y'
    if (bond_r1 < bond_p1): return 'n'


def r_to_p(n1, n2, lines, n_atoms, n_idx, atom, fileout_reorder):
    for i in range(0, n2):
        for j in range(0, n_atoms + 2): fileout_reorder.write(lines[(n_idx - 1 - i) * (n_atoms + 2) + j])
    for i in range(0, n1):
        for j in range(0, n_atoms + 2): fileout_reorder.write(lines[i * (n_atoms + 2) + j])
    print('R to P')
    return 1


def p_to_r(n1, n2, lines, n_atoms, n_idx, atom, fileout_reorder):
    for i in range(0, n1):
        for j in range(0, n_atoms + 2): fileout_reorder.write(lines[(n1 - 1 - i) * (n_atoms + 2) + j])
    for i in range(0, n2):
        for j in range(0, n_atoms + 2): fileout_reorder.write(lines[(i + n1) * (n_atoms + 2) + j])
    print('P to R')
    return 1


def output_TS(lines, n_atoms, atom, fileout_TS, fileout_TS_xyz):
    bond_ts1, bond_ts2 = get_distance(lines, 0, n_atoms, atom)
    fileout_TS.write('1 , ' + str(bond_ts1) + ' , ' + str(bond_ts2) + "\n")
    for i in range(0, n_atoms + 2): fileout_TS_xyz.write(lines[i])
    return 1


def Comb(filename, atom):
    print('Working on ', filename)
    if os.stat(filename).st_size == 0:
        print('The file is empty !')
        return 1
    lines, n_lines, n_atoms, n_idx, filein, fileout_TS_xyz, fileout_TS, fileout_reorder = open_file(filename)
    if len(lines) == 0: return 1
    n1, n2 = looking_for_ts(lines, n_lines, n_atoms, n_idx)
    if n1 == n2 == 1: return 1
    output_TS(lines, n_atoms, atom, fileout_TS, fileout_TS_xyz)
    filein.close()
    fileout_TS.close()
    fileout_reorder.close()
    return 1


# All following functions are used for generating hist diagram
def hist(filename):
    data = getdata(filename)
    ax_hist = initialization_hist()
    making2hist(ax_hist, data, 'm', 'n')
    plt.savefig('./trajTS/TS zone histogram.png', dpi=100)
    mean_m, interval_m = mean_interval(data['m'])
    mean_n, interval_n = mean_interval(data['n'])
    fileout_hist = open('./trajTS/TS_zone', 'w')
    fileout_hist.write(str("%.2f" % round(interval_m[0], 2)) + "\n")
    fileout_hist.write(str("%.2f" % round(interval_m[1], 2)) + "\n")
    fileout_hist.write(str("%.2f" % round(interval_n[0], 2)) + "\n")
    fileout_hist.write(str("%.2f" % round(interval_n[1], 2)) + "\n")
    fileout_hist.write('mean_m ' + str("%.2f" % round(mean_m, 2)) + "\n")
    fileout_hist.write('mean_n ' + str("%.2f" % round(mean_n, 2)) + "\n")
    fileout_hist.close()
    return "%.2f" % round(interval_m[0], 2), "%.2f" % round(interval_m[1], 2), "%.2f" % round(interval_n[0], 2), "%.2f" % round(interval_n[1], 2)


def making2hist(ax, data, a, b):
    n1, bins1, patches1 = ax.hist(data[a], bins=10, color='b', normed=1,
                                  weights=np.zeros_like(data[a]) + 1. / data[a].size, label='Bond 1')
    n2, bins2, patches2 = ax.hist(data[b], bins=10, color='r', normed=1,
                                  weights=np.zeros_like(data[b]) + 1. / data[b].size, label='Bond 2')
    ax.legend()
    return bins1, bins2


def getdata(filename):
    if os.stat(filename).st_size == 0:
        print('The file is empty')
        return 1
    data = np.genfromtxt(filename, delimiter=',',
                         dtype=[('l', 'float32'),
                                ('m', 'float32'),
                                ('n', 'float32')])
    return data


def initialization_hist():
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    ax.set_xlabel('Bond Length ($\AA$)', size=15)
    ax.set_ylabel('Frequency (%)', size=15)
    return ax


def mean_interval(data):
    mean = np.mean(data)
    mu = np.std(data)
    interval = norm.interval(alpha=0.98, loc=mean, scale=mu)
    return mean, interval


# main func
def main():
    atom = get_atomindex()
    mkdir()
    for filename in glob.glob('./ntraj/*.xyz'): Comb(filename, atom)
    print('Done trine!')
    TS = hist('./trajTS/trajTs.txt')
    print('Done histogram :', TS)
    print('Complete !')


if __name__ == "__main__":
    main()
