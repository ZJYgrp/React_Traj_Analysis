# This script should be revised to an object-oriented program, which defines TS extraction class, trajectories classfication class, etc.
import numpy as np
import sys
import os
import glob
from matplotlib import pyplot
from scipy import stats
from scipy.stats import norm
from XYZtraj_Trajs import Trajectories 

def mkdir():
        os.system('rm -rf ./trajTS')
        os.system('rm -rf ./reorder')
        os.system('rm -rf ./TDD')
        os.system('rm -rf ./TDD_r2p*')
        os.system('rm -rf ./TDD_r2r')
        os.system('rm -rf ./TDD_p2p*')
        os.system('rm -rf ./TDD_inter')
        os.system('mkdir ./trajTS')
        os.system('mkdir ./reorder')
        os.system('mkdir ./TDD')
        os.system('mkdir ./TDD_r2pA')
        os.system('mkdir ./TDD_r2pB')
        os.system('mkdir ./TDD_r2r')
        os.system('mkdir ./TDD_p2pA')
        os.system('mkdir ./TDD_p2pB')
        os.system('mkdir ./TDD_inter')
        return 1

def read_conf_file():
    f = open('traj.conf')
    lines = [a.strip() for a in f.readlines()]
    f.close()
    reset = lines[0]
    if not lines[1].isdigit():
        print('Invalid mode in conf file, must be a number')
        exit(1)
    else:
        mode = int(lines[1])
        if mode < 1 or mode > 2:
            print('Invalid mode in conf file, too large or too small')
            exit(1)
    atoms = []
    for num in lines[2].split():
        if not num.isdigit():
            print('Invalid atom index')
            exit(1)
        else:
            atoms.append(int(num))
    if (len(atoms) % 2) != 0:
        raise TypeError('Odd number of atomic indices have been received–this is ODD!')
    elif (mode == 1 and len(atoms) != 6) or (mode == 2 and len(atoms) != 4):
        print('Invalid number of atoms. Exiting program')
        exit(1)
    print ('Run ...')
    return reset, mode, atoms


def Get_mode():
    print('Available modes:\nMode 1:\t1 bond always forms, then 1 of 2 other bonds forms to create 2 products\nMode 2:\t1 of 2 possible bonds form to create 2 products\n')
    mode = input('Please choose analysis mode: ')
    if not mode.isdigit() or int(mode) < 1 or int(mode) > 2:
        print('Invalid mode selection. Exiting program')
        exit(1)
    return int(mode)

def Get_atomindex(mode):
    if mode == 1:
        print('Enter indices for bond that always forms, then for the bonds corresponding to products A and B')
    elif mode == 2:
        print('Enter indices for bonds corresponding to product A, then product B')
    atoms = [int(x) for x in
            input('Input atomic indices corresponding to N bonds, using ' ' as delimiter. Make sure the forming bond goes first\n').split()]
    for i in range(0,len(atoms)):
        print('atom ', str(i+1), ' ', str(atoms[i]))
    judge = input('Do you think the indices are reasonable?(y/n): ')
    if judge != 'y':
        print('I have to quit, sorry!')
        exit(1)
    if (len(atoms) % 2) != 0:
        raise TypeError('Odd number of atomic indices have been received–this is ODD!')
    elif (mode == 1 and len(atoms) != 6) or (mode == 2 and len(atoms) != 4):
        print('Invalid number of atoms. Exiting program')
        exit(1)
    print ('Run ...')
    return atoms

def log_results(A, B, revert, inter, total):
  out = open('./trajTS/traj_log', 'w+')
  out.write('Results\nTotal number of trajectories: '+str(total)+'\nTotal forming product: '+str(A+B)+'\nA: '+str(A)+' B: '+str(B)+' Reactant: '+str(revert)+' Intermediate: '+str(inter)+'\nPercent product A: '+str((A*100/(A+B)))+'%\nPercent product B: '+str(B*100/(A+B))+'%')

def plot_gaps(X, Y):

    #bins = np.linspace(min(X) + 10, max(Y) + 10, 80)
    xAB = list(map(abs, X))
    yAB = list(map(abs, Y))

    bins = np.arange(min(xAB), max(yAB) + 5, 5)
    pyplot.hist(xAB, bins, alpha=0.5, label='X')
    pyplot.hist(yAB, bins, alpha=0.5, label='Y')
    pyplot.legend(loc='upper right')
    pyplot.axvline(x=np.median(xAB), color='b', linestyle='dashed', linewidth=2)
    pyplot.axvline(x=np.median(yAB), color='orange', linestyle='dashed', linewidth=2)
    #pyplot.title('Bond Formation Time Gaps')
    #pyplot.xlabel('Time Gap (fs)')
    #pyplot.ylabel('Frequency')
    #pyplot.show()
    pyplot.savefig('plt.png')

def plot_bonds(X, Y):

    xAB = list(map(abs, X))
    yAB = list(map(abs, Y))
    #bins = np.linspace(min(X) + 10, max(Y) + 10, 80)
    bins = np.arange(min(xAB), max(yAB) + 5, 5)

    pyplot.hist(xAB, bins, alpha=0.5, label='X')
    pyplot.hist(yAB, bins, alpha=0.5, label='Y')
    pyplot.legend(loc='upper right')
    pyplot.axvline(x=np.median(xAB), color='b', linestyle='dashed', linewidth=2)
    pyplot.axvline(x=np.median(yAB), color='orange', linestyle='dashed', linewidth=2)
    #pyplot.title('Bond Formation Times')
    #pyplot.xlabel('Formation Time (fs)')
    #pyplot.ylabel('Frequency')
    #pyplot.show()
    pyplot.savefig('plt.png')


# main func
def main():
# Remember to add a choice function regarding the removal of current folders
    if os.path.exists('traj.conf'):
        reset, mode, atom = read_conf_file()
        if reset == 'y': mkdir()
    else:
        judge = input('Do you want to start analyzing from the very beginning? (y/n) Type y to remove all analysis folders and n to keep the current folder (e.g. reorder, etc.) for analysis: ')
        if judge == 'y': mkdir()
        mode = Get_mode()
        atom = Get_atomindex(mode)
# The attribute of Trajectories Class involves
#        self.Get_distance(1)
#        self.TS_finder()
#        self.Rearrangement()
    for filename in glob.glob('./ntraj/*.xyz'):
        Trajectories(filename, atom, mode).TS_finder()
        Trajectories(filename, atom, mode).Rearrangement()
#       self.Classification() Note: Classification has to be called after Rearrangement or when reorder has been performed.
    A, B = 0, 0
    time_file = open('./trajTS/trajTimes.txt', 'w')
    if mode == 1:
        bond1_times, bond2_times, bond3_times = [], [], []
        gapA, gapB = [], []
        totalA, totalB, concertedA, concertedB = 0, 0, 0, 0
        for filename in glob.glob('./reorder/*.xyz'):
            result = Trajectories(filename, atom, mode).Formation_Time()
            if type(result) is list:
                time_file.write(', '.join([os.path.basename(filename)]+list(map(str, result)))+'\n')
                
                if result[0] == 'A':
                    bond1_times.append(result[1])
                    bond2_times.append(result[2])
                    gapA.append(result[2] - result[1])
                    totalA += 1
                    if result[2]-result[1]<60:
                        concertedA += 1
                elif result[0] == 'B':
                    bond1_times.append(result[1])
                    bond3_times.append(result[2])     
                    gapB.append(result[2] - result[1])
                    totalB += 1
                    if result[2]-result[1]<60:
                        concertedB += 1
        print('Trajectory analysis complete!')
        print('Median time for bond 1: ', np.median(list(map(abs, bond1_times))))
        print('Median time for bond 2: ', np.median(list(map(abs, bond2_times))))
        print('Median time for gapA: ', np.median(list(map(abs, gapA))))
        print('Median time for bond 3: ', np.median(list(map(abs, bond3_times))))
        print('Median time for gapB: ', np.median(list(map(abs, gapB))))
        print('Percentage dynamically concerted A: ', concertedA/totalA*100)
        print('Percentage dynamically concerted B: ', concertedB/totalB*100)
        # print(','.join(list(map(str, [sum(bond1_times)/len(bond1_times), sum(bond2_times)/len(bond2_times), sum(gapA)/len(gapA), sum(bond3_times)/len(bond3_times), sum(gapB)/len(gapB)]))))
        plot_gaps(gapA, gapB)
    elif mode == 2:
        bond1_times, bond2_times = [], []
        for filename in glob.glob('./reorder/*.xyz'):
            result = Trajectories(filename, atom, mode).Formation_Time()
            if type(result) is list:
                time_file.write(', '.join([os.path.basename(filename)]+list(map(str, result)))+'\n')
                if result[0] == 'A':
                    bond1_times.append(result[1])
                elif result[0] == 'B':
                    bond2_times.append(result[1])     
                    if (result[1] > 300):
                        print(filename, result)
        print('Trajectory analysis complete!')
        print('Median time for bond 1: ', np.median(list(map(abs, bond1_times))))
        print('Median time for bond 2: ', np.median(list(map(abs, bond2_times))))
        plot_bonds(bond1_times, bond2_times)
    time_file.close()


if __name__ == '__main__':
    main()
                                                 
