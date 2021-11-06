import glob
import os
import numpy as np
from trajectories import Trajectories


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
    print('Run ...')
    return reset, mode, atoms


def get_mode():
    print('Available modes:\n'
          'Mode 1:\t1 bond always forms, then 1 of 2 other bonds forms to create 2 products\n'
          'Mode 2:\t1 of 2 possible bonds form to create 2 products\n')
    mode = input('Please choose analysis mode: ')
    if not mode.isdigit() or int(mode) < 1 or int(mode) > 2:
        print('Invalid mode selection. Exiting program')
        exit(1)
    return int(mode)


def get_atom_index(mode):
    if mode == 1:
        print('Enter indices for bond that always forms, then for the bonds corresponding to products A and B')
    elif mode == 2:
        print('Enter indices for bonds corresponding to product A, then product B')
    atoms = [int(x) for x in
             input('Input atomic indices corresponding to N bonds, using ' ' as delimiter. '
                   'Make sure the forming bond goes first\n').split()]
    for i in range(0, len(atoms)):
        print('atom ', str(i + 1), ' ', str(atoms[i]))
    judge = input('Do you think the indices are reasonable?(y/n): ')
    if judge != 'y':
        print('I have to quit, sorry!')
        exit(1)
    if (len(atoms) % 2) != 0:
        raise TypeError('Odd number of atomic indices have been received–this is ODD!')
    elif (mode == 1 and len(atoms) != 6) or (mode == 2 and len(atoms) != 4):
        print('Invalid number of atoms. Exiting program')
        exit(1)
    print('Run ...')
    return atoms


def log_results_mode1(bond1_times, bond2_times, bond3_times, gapX, gapY, concertedX, concertedY, totalX, totalY):
    out = open('./trajTS/time_log', 'w+')
    output = 'Results\n' \
             'Median time for bond 1: {0}\n' \
             'Median time for bond 2: {1}\n' \
             'Median time for gapX: {2}\n' \
             'Median time for bond 3: {3}\n' \
             'Median time for gapY: {4}\n' \
             'Percentage dynamically concerted X:: {5:.1f}%\n' \
             'Percentage dynamically concerted Y: {6:.1f}'.format(
        np.median(list(map(abs, bond1_times))),
        np.median(list(map(abs, bond2_times))),
        np.median(list(map(abs, gapX))),
        np.median(list(map(abs, bond3_times))),
        np.median(list(map(abs, gapY))),
        concertedX / totalX * 100,
        concertedY / totalY * 100)
    out.write(output)
    out.close()
    return output


def log_results_mode2(bond1_times, bond2_times):
    out = open('./trajTS/time_log', 'w+')

    output = 'Results\n' \
             'Median time for bond 1: {0}\n' \
             'Median time for bond 2: {1}'.format(
        np.median(list(map(abs, bond1_times))),
        np.median(list(map(abs, bond2_times))))
    out.write(output)
    out.close()
    return output


# main func
def main():
    if os.path.exists('traj.conf'):
        reset, mode, atom = read_conf_file()
        if reset == 'y': mkdir()
    else:
        judge = input('Do you want to start analyzing from the very beginning? (y/n) '
                      'Type y to remove all analysis folders and n to keep the current folder '
                      '(e.g. reorder, etc.) for analysis: ')
        if judge == 'y': mkdir()
        mode = get_mode()
        atom = get_atom_index(mode)

    for filename in glob.glob('./ntraj/*.xyz'):
        Trajectories(filename, atom, mode).TS_finder()
        Trajectories(filename, atom, mode).rearrangement()

    time_file = open('./trajTS/trajTimes.txt', 'w')
    if mode == 1:
        bond1_times, bond2_times, bond3_times = [], [], []
        gapX, gapY = [], []
        totalX, totalY, concertedX, concertedY = 0, 0, 0, 0
        for filename in glob.glob('./reorder/*.xyz'):
            result = Trajectories(filename, atom, mode).formation_time()
            if type(result) is list:
                time_file.write(', '.join([os.path.basename(filename)] + list(map(str, result))) + '\n')

                if result[0] == 'X':
                    bond1_times.append(result[1])
                    bond2_times.append(result[2])
                    gapX.append(result[2] - result[1])
                    totalX += 1
                    if result[2] - result[1] < 60:
                        concertedX += 1
                elif result[0] == 'Y':
                    bond1_times.append(result[1])
                    bond3_times.append(result[2])
                    gapY.append(result[2] - result[1])
                    totalY += 1
                    if result[2] - result[1] < 60:
                        concertedY += 1
        print('Trajectory analysis complete!')
        output = log_results_mode1(bond1_times, bond2_times, bond3_times, gapX, gapY, concertedX, concertedY, totalX,
                                   totalY)
        print(output)

        # plot_gaps(gapA, gapB)
    elif mode == 2:
        bond1_times, bond2_times = [], []
        for filename in glob.glob('./reorder/*.xyz'):
            result = Trajectories(filename, atom, mode).formation_time()
            if type(result) is list:
                time_file.write(', '.join([os.path.basename(filename)] + list(map(str, result))) + '\n')
                if result[0] == 'X':
                    bond1_times.append(result[1])
                elif result[0] == 'Y':
                    bond2_times.append(result[1])

        print('Trajectory analysis complete!')
        output = log_results_mode2(bond1_times, bond2_times)
        print(output)
        # plot_bonds(bond1_times, bond2_times)
    time_file.close()


if __name__ == '__main__':
    main()
