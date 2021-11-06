import glob
import os

from matplotlib.ticker import MultipleLocator
from scipy.stats import norm
import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy import *
# from mpl_toolkits.mplot3d import Axes3D
# import matplotlib.patches as mpatches

mpl.use('Agg')


# TS sampling
def open_tsfile(file_name):
    filein = open(file_name, 'r')
    lines = filein.readlines()
    n_lines = len(lines)
    n_atoms = int(lines[0].split()[0])
    n_idx = int((n_lines / (n_atoms + 2)))
    fileout_TS = open('./trajTS/trajTs.txt', 'a')
    return lines, n_lines, n_atoms, n_idx, filein, fileout_TS


def TS(file_name, atom):
    print('Working on ', file_name)
    if os.stat(file_name).st_size == 0:
        print('The file is empty !')
        return 1
    lines, n_lines, n_atoms, n_idx, filein, fileout_TS = open_tsfile(file_name)
    if len(lines) == 0: return 1
    fileout_TS = outputs(lines, n_atoms, n_idx, atom, fileout_TS)
    fileout_TS.close()
    filein.close()
    return 1


def outputs(lines, n_atoms, n_idx, atom, fileout_ts):
    for i in range(0, n_idx):
        bond1, bond2, bond3 = get_distance(lines, i, n_atoms, atom)
        fileout_ts.write(str(i) + ' , ' + str(bond1) + ' , ' + str(bond2) + ' , ' + str(bond3) + '\n')
    return fileout_ts


def hist(file_name):
    data = get_data(file_name)
    #  ax_hist=initialization_hist()
    #  making2hist(ax_hist,data,'m','n')
    #  plt.savefig('./trajTS/TS zone histogram.png',dpi=100)
    mean_m, interval_m = mean_interval(data['m'])
    mean_n, interval_n = mean_interval(data['n'])
    mean_o, interval_o = mean_interval(data['o'])
    fileout_hist = open('./trajTS/TS_zone', 'w')
    fileout_hist.write(str('%.2f' % round(interval_m[0], 2)) + '\n')
    fileout_hist.write(str('%.2f' % round(interval_m[1], 2)) + '\n')
    fileout_hist.write(str('%.2f' % round(interval_n[0], 2)) + '\n')
    fileout_hist.write(str('%.2f' % round(interval_n[1], 2)) + '\n')
    fileout_hist.write(str('%.2f' % round(interval_o[1], 2)) + '\n')
    fileout_hist.write(str('%.2f' % round(interval_o[1], 2)) + '\n')
    fileout_hist.write('mean_m ' + str('%.2f' % round(mean_m, 2)) + '\n')
    fileout_hist.write('mean_n ' + str('%.2f' % round(mean_n, 2)) + '\n')
    fileout_hist.write('mean_o ' + str('%.2f' % round(mean_o, 2)) + '\n')
    fileout_hist.close()
    return '%.2f' % round(interval_m[0], 2), '%.2f' % round(interval_m[1], 2), '%.2f' % round(interval_n[0],
                                                                                              2), '%.2f' % round(
        interval_n[1], 2)


def get_data(file_name):
    if os.stat(file_name).st_size == 0:
        print('The file is empty')
        return 1

    f = open(file_name)
    length = len(f.readline().split())
    f.close()
    data = genfromtxt(file_name, skip_header=1, delimiter=',', dtype=None, encoding=None,
                         names=[chr(i) for i in range(ord('l'), ord('l') + length)])
    #  data = np.genfromtxt(filename, skip_header=1 ,delimiter=',', dtype=None, encoding=None, names=['l', 'm', 'n', 'o'])
    #  data = np.genfromtxt(filename, skip_header=1 ,delimiter=',',
    #	    dtype=[('l','S11'),('m','float32'),('n','float32'),('o','float32')])
    return data


def initialization_hist():
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    ax.set_xlabel('Bond Length ($\AA$)', size=15)
    ax.set_ylabel('Probability (%)', size=15)
    return ax


def making2hist(ax, data, a, b):
    n1, bins1, patches1 = ax.hist(data[a], bins=10, color='b', normed=1,
                                  weights=zeros_like(data[a]) + 1. / data[a].size, label='Bond 1')
    n2, bins2, patches2 = ax.hist(data[b], bins=10, color='r', normed=1,
                                  weights=zeros_like(data[b]) + 1. / data[b].size, label='Bond 2')
    ax.legend()
    return bins1, bins2


def mean_interval(data):
    data_mean = mean(data)
    mu = std(data)
    interval = norm.interval(alpha=0.98, loc=data_mean, scale=mu)
    return data_mean, interval


# all the following files are used for extracting  Time-distance-distance(TDD)

def open_file(file_name):
    filein = open(file_name, 'r')
    lines = filein.readlines()
    n_lines = len(lines)
    n_atoms = int(lines[0].split()[0])
    n_idx = int((n_lines / (n_atoms + 2)))
    return lines, n_lines, n_atoms, n_idx, filein


def get_distance(lines, n, n_atoms, atom):
    x = [0, 0, 0, 0, 0, 0]
    y = [0, 0, 0, 0, 0, 0]
    z = [0, 0, 0, 0, 0, 0]
    for i in range(0, 6): x[i] = float(lines[(n_atoms + 2) * n + atom[i] + 1].split()[1])
    for i in range(0, 6): y[i] = float(lines[(n_atoms + 2) * n + atom[i] + 1].split()[2])
    for i in range(0, 6): z[i] = float(lines[(n_atoms + 2) * n + atom[i] + 1].split()[3])
    bond1 = round(((x[1] - x[0]) ** 2 + (y[1] - y[0]) ** 2 + (z[1] - z[0]) ** 2) ** .5, 3)
    bond2 = round(((x[3] - x[2]) ** 2 + (y[3] - y[2]) ** 2 + (z[3] - z[2]) ** 2) ** .5, 3)
    bond3 = round(((x[5] - x[4]) ** 2 + (y[5] - y[4]) ** 2 + (z[5] - z[4]) ** 2) ** .5, 3)
    return bond1, bond2, bond3


def get_lines(ax, data_len):
    if data_len == 4:
        for traj in glob.glob('./TDD_r2pX/traj*.xyz.txt'):
            if os.stat(traj).st_size == 0: break
            data1 = get_data(traj)
            ax.plot(data1['n'], data1['o'], 'g-', lw=2.0, zorder=-1000)
        for traj in glob.glob('./TDD_r2pY/traj*.xyz.txt'):
            if os.stat(traj).st_size == 0: break
            data1 = get_data(traj)
            ax.plot(data1['n'], data1['o'], 'b-', lw=2.0, zorder=1000)
    elif data_len == 3:
        for traj in glob.glob('./TDD_r2pX/traj*.xyz.txt'):
            if os.stat(traj).st_size == 0: break
            data1 = get_data(traj)
            ax.plot(data1['m'], data1['n'], 'g-', lw=2.0, zorder=-1000)
        for traj in glob.glob('./TDD_r2pY/traj*.xyz.txt'):
            if os.stat(traj).st_size == 0: break
            data1 = get_data(traj)
            ax.plot(data1['m'], data1['n'], 'b-', lw=2.0, zorder=1000)

    '''
    for traj in glob.glob('./TDD_r2pA/traj*.xyz.txt'):
        if os.stat(traj).st_size == 0: break
        data1 = get_data(traj)
        if len(data1[0]) == 4:
            # ax.plot(data1['n'],data1['m'],data1['o'],'g-',lw=2.0,zorder=-1000)
            ax.plot(data1['n'], data1['o'], 'g-', lw=2.0, zorder=-1000)
        elif len(data1[0]) == 3:
            ax.plot(data1['m'], data1['n'], 'g-', lw=2.0, zorder=-1000)
    for traj in glob.glob('./TDD_r2pB/traj*.xyz.txt'):
        if os.stat(traj).st_size == 0: break
        data1 = get_data(traj)
        if len(data1[0]) == 4:
            # ax.plot(data1['n'],data1['m'],data1['o'],'b-',lw=2.0,zorder=1000)
            ax.plot(data1['n'], data1['o'], 'b-', lw=2.0, zorder=1000)
        elif len(data1[0]) == 3:
            ax.plot(data1['m'], data1['n'], 'b-', lw=2.0, zorder=1000)
    '''
    return ax


def get_scatter(data, ax):
    if len(data[0]) == 4:
        # ax.scatter(data['n'],data['m'],data['o'], c='r', alpha=0.3,s=50,zorder=1001)
        ax.scatter(data['n'], data['o'], c='r', alpha=0.3, s=50, zorder=1001)
    elif len(data[0]) == 3:
        ax.scatter(data['m'], data['n'], c='r', alpha=0.3, s=50, zorder=1001)
    return data


# Trajectory graph D vs D, share functions with D vs. T
def dvd(filename):
    data = get_data('./trajTS/trajTs.txt')
    ax = initialization_dvd(data)
    get_scatter(data, ax)
    get_lines(ax, len(data[0]))
    plt.savefig('./trajTS/' + filename.split('/')[1] + '_dvd.png', dpi=100)


def initialization_dvd(data):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    if len(data[0]) == 4:
        # ax=fig.gca(projection='3d')

        ax.set_xlim(1.5, 4.0)
        ax.set_ylim(1.5, 4.0)
        # ax.set_zlim(1.3, 5.0)
        ax.set_title('Bond Distances')

        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_major_locator(MultipleLocator(0.5))
        # ax.zaxis.set_major_locator(MultipleLocator(0.5))

        ax.grid(color='#000000', linewidth=1, linestyle=':')
        # ax.xaxis._axinfo['grid'].update({'linewidth':1, 'color' : 'black'})

        # ax.yaxis._axinfo['grid'].update({'linewidth':1, 'color' : 'black'})
        # ax.zaxis._axinfo['grid'].update({'linewidth':1, 'color' : 'black'})
        # ax.zaxis._axinfo['grid']['linestyle'] = ':'
        ax.set_xlabel('Bond 2 (' + r'$\mathrm{\AA}$' + ')')  # r'$\mathrm{\AA}$' is unitalicized Angstrom abbr.
        ax.set_ylabel('Bond 3 (' + r'$\mathrm{\AA}$' + ')')
        # ax.set_zlabel('Bond C (' + r'$\mathrm{\AA}$' + ')')
    #    ax.w_xaxis._axinfo.update({'grid' : {'color': (0, 0, 0, 0.3)}})
    #    ax.w_yaxis._axinfo.update({'grid' : {'color': (0, 0, 0, 0.3)}})
    #    ax.w_zaxis._axinfo.update({'grid' : {'color': (0, 0, 0, 0.3)}})
    # ax.invert_xaxis()
    elif len(data[0]) == 3:
        # ax=fig.add_axes([0, 0, 1, 1])
        ax.set_xlim(1.5, 4.0)
        ax.set_ylim(1.5, 4.0)
        ax.grid(color='#000000', linewidth=1, linestyle=':')
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_major_locator(MultipleLocator(0.5))

        ax.set_xlabel('Bond 1 (' + r'$\mathrm{\AA}$' + ')')  # r'$\mathrm{\AA}$' is unitalicized Angstrom abbr.
        ax.set_ylabel('Bond 2 (' + r'$\mathrm{\AA}$' + ')')
        ax.set_title('Bond Distances')
    else:
        print('I don\'t know how to plot this data!')
        exit(1)

    # ax.scatter( 2.994 ,2.677 , 2.008,c='k', s=50)
    # ax.scatter( 2.676 ,3.066 , 1.985,c='k', s=50)
    # ax.scatter( 2.65 ,2.65 , 1.64,c='g', s=50)
    # p1 = Ellipse((2.95, 3.1), 1.0, 0.9,-45,facecolor = 'y',alpha=0.7)
    # ax.add_patch(p1)
    # art3d.pathpatch_2d_to_3d(p1, z=2.05,zdir='z')
    # p2 = Ellipse((2.9, 3.1),0.3,1.2,45, facecolor = 'g',alpha=0.7)
    # ax.add_patch(p2)
    # art3d.pathpatch_2d_to_3d(p2, z=1.4)
    # p3 = Circle((3.7, 2.0),0.5, facecolor = 'k',alpha=0.2,zorder=-2000)
    # ax.add_patch(p3)
    # art3d.pathpatch_2d_to_3d(p3, z=1.3)
    # p4 = Circle((2.0, 3.7),0.5, facecolor = 'k',alpha=0.2,zorder=-2000)
    # ax.add_patch(p4)
    # art3d.pathpatch_2d_to_3d(p4, z=1.3)
    # p5 = Circle((3.5, 3.0),1.0, facecolor = 'k',alpha=0.3,zorder=-2000)
    # ax.add_patch(p5)
    # art3d.pathpatch_2d_to_3d(p5, z=4.0)
    # ax.set_title('Distance vs Distance')
    # plt.xlabel('Distance ($\AA$)')
    # plt.ylabel('Distance ($\AA$)')
    return ax


# main func
def main():
    if len(os.listdir('./TDD_r2pX')) > 0 or len(os.listdir('./TDD_r2pY')) > 0:
        print('Now, making r2p_DVD plot')
        dvd('./TDD_r2p/*.xyz.txt', './trajTS/trajTs.txt')
        print('Complete!')
    else:
        print('No files found in TDD_r2p! What am I supposed to analyze?')


if __name__ == '__main__':
    main()
