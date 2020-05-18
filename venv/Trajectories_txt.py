import sys
import os
import numpy as np
import glob
from scipy.stats import norm
from itertools import combinations
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import LineCollection
from matplotlib.collections import PatchCollection
from matplotlib.colors import colorConverter
from matplotlib.ticker import FuncFormatter
#from scipy import stats
#from scipy.stats import norm
# Time gap function reads in the TDD_r2p trajectories and
# analyze the timing of bond formation
class Trajectories_txt:
    def __init__(self, file1, file2, N_col=3):
        # trajectory format for ProgDyn output. The four criteria need to take different tests.
        if (not os.path.exists(file1+'trajTs.txt')) and (not os.listdir(file2)):
            print('A trajTs.txt or traj*.xyz.txt is not found. Analysis can not proceed.')
            return
        elif os.path.exists(file1+'trajTs.txt') and (not file2):
            if os.stat(file1 + 'trajTs.txt').st_size == 0:
                print('trajTs.txt is found but is empty. Analysis can not proceed.')
                return
            else:
                # Creating new folders for the following analysis
                # Open new file handles and get parameters, delimiter might be hard to find, takes trial and error
                data = np.genfromtxt(file1+'trajTs.txt', delimiter=',', autostrip=True, dtype=float)
                # using isfinite can help remove all nan,inf,etc. But be cautious about the dtype option.
                data2 = np.where(np.isfinite(data), data, 0)
                # Reshape is necessary to regularize the data format, just in case only one line of data is involved.
                Num_lines = len(open(file1+'trajTs.txt', 'r').readlines())
                data2 = np.reshape(data2, (Num_lines, -1))
                col = [i for i in range(0, len(data2[0])) if np.mean(data2, axis=0)[i] != 0]
                self.TS_data = data2[:, col]
        elif (not os.path.exists(file1+'trajTs.txt')) and os.listdir(file2):
            comb = [subset for subset in combinations(np.arange(N_col), 2)]
            self.lines_col=[ [] for i in range(0,len(comb))]
            for filename in glob.glob(file2+"*.xyz.txt"):
                if os.stat(filename).st_size == 0:
                    print(Filename+' is found but is empty. Jump to next.')
                    break
                else:
                    # Creating new folders for the following analysis
                    # Open new file handles and get parameters, delimiter might be hard to find, takes trial and error
                    data = np.genfromtxt(filename, delimiter=',', autostrip=True)
                    for j in range(0, len(comb)):
                        self.lines_col[j].append(zip(data[:,comb[j][0]], data[:,comb[j][1]]))

    def TS_plot(self, custom_bin_num=0, custom_range=[]):
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_subplot(111)
        ax.set_xlabel('Bond Length ($\AA$)', fontname="Arial", size=12)
        ax.set_ylabel('Frequency (%)', fontname="Arial", size=12)
        color_code=['b', 'r', 'g', 'k', 'm', 'y', 'c']
        # if a customized range has been proposed
        if custom_range:
            hist_range = custom_range
            bin_num = custom_bin_num
            print("Cutomized histogram bins applied,"+str(bin_num))
            print("Cutomized histogram range applied, from "+str(custom_range[0])+" to "+str(custom_range[1]))
        # if a customized range was not provided
        else:
            TS_MIN = np.min(self.TS_data, axis=0)
            TS_MAX = np.max(self.TS_data, axis=0)
            gap=(TS_MAX[0]-TS_MIN[0])/20
            bin_num= int((np.max(TS_MAX)-np.min(TS_MIN))/gap)
            hist_range = [np.min(TS_MIN),np.max(TS_MAX)]
            print("Default number of histogram bins applied," + str(bin_num))
            print("Default range applied, from " + str(hist_range[0]) + " to " + str(hist_range[1]))
        for i in range(0,len(self.TS_data[0])):
            ax.hist(self.TS_data[:,i], bins=bin_num, range = hist_range,
                    color = color_code[i], density=True, label='Bond '+str(i))
        plt.xticks(np.arange(round(hist_range[0],1), round(hist_range[1],1), step=(hist_range[1]-hist_range[0])/5))
        ax.set_xticklabels([str("%.2f" % round(float(label), 2)) for label in
                            np.arange(round(hist_range[0],1), round(hist_range[1],1), step=(hist_range[1]-hist_range[0])/5)])
        ax.legend()
        plt.savefig('./trajTS/TS zone histogram.pdf', dpi=300)
        fileout=open('./trajTS/TS_zone.dat','w')
        mean=np.mean(self.TS_data, axis=0)
        mu=np.std(self.TS_data, axis=0)
        Conf_Int = norm.interval(alpha=0.997, loc=mean, scale=mu)
        for i in range(0,len(Conf_Int)):
            print("Printing TS information for Bond "+str(i+1))
            fileout.write('Bond '+str(i+1)+' info: \n')
            fileout.write('Mean is: ' + str(mean[i]) + '\n')
            fileout.write('Std. Dv. is: ' + str(mu[i]) + '\n')
            fileout.write('TS zone (99.7%) upper bound is: ' + str(Conf_Int[i][0]) + '\n')
            fileout.write('TS zone (99.7%) lower bound is: ' + str(Conf_Int[i][1]) + '\n')
        fileout.close()

    def Traj_plot(self):
        for i in range(0,len(self.lines_col)):
            Pass

