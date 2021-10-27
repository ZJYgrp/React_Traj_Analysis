import os
import glob
from trajectories import Trajectories


def mkdir():
    os.system('rm -r ./trajTS')
    os.system('rm -r ./reorder')
    os.system('rm -r ./TDD')
    os.system('rm -r ./TDD_r2p*')
    os.system('rm -r ./TDD_r2r')
    os.system('rm -r ./TDD_p2p*')
    os.system('rm -r ./TDD_inter')
    os.system('mkdir ./trajTS')
    os.system('mkdir ./reorder')
    os.system('mkdir ./TDD')
    os.system('mkdir ./TDD_r2pX')
    os.system('mkdir ./TDD_r2pY')
    os.system('mkdir ./TDD_r2r')
    os.system('mkdir ./TDD_p2pX')
    os.system('mkdir ./TDD_p2pY')
    os.system('mkdir ./TDD_inter')
    return 1


def read_conf_file():
    conf_file = open('traj.conf')
    lines = [a.strip() for a in conf_file.readlines()]
    conf_file.close()
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
    if (mode == 1 and len(atoms) != 6) or (mode == 2 and len(atoms) != 4):
        print('Invalid number of atoms. Exiting program')
        exit(1)
    return reset, mode, atoms


def get_mode():
    print('Available modes:\n'
          'Mode 1:\t1 bond always forms, then 1 of 2 other bonds forms to create 2 products\n'
          'Mode 2:\t1 of 2 possible bonds form to create 2 products\n')
    mode = input('Enter mode: ')
    while not mode.isdigit() or mode < 1 or mode > 2:
        mode = input('Invalid mode. Try again: ')
    return mode


def get_atomindex(mode):
    line = input("Enter space-separated atom indices: ")
    atoms = []
    for num in line.split():
        if not num.isdigit():
            print('Invalid atom index')
            exit(1)
        else:
            atoms.append(int(num))
    if (len(atoms) % 2) != 0:
        raise TypeError('Odd number of atomic indices have been received–this is ODD!')
    if (mode == 1 and len(atoms) != 6) or (mode == 2 and len(atoms) != 4):
        print('Invalid number of atoms. Exiting program')
        exit(1)
    return atoms


def log_results(X, Y, revert, inter, total):
    out = open('./trajTS/traj_log', 'w+')
    output = 'Results\n' \
             'Total number of trajectories: {0}\n' \
             'Total forming product: {1}\n' \
             'X: {2} Y: {3} Reactant: {4} Intermediate: {5}\n' \
             'Percent product X: {6}%\n' \
             'Percent product Y: {7}%'.format(
                 str(total),
                 str(X + Y),
                 str(X), str(Y), str(revert), str(inter),
                 str((X * 100 / (X + Y))),
                 str(Y * 100 / (X + Y)))
    out.write(output)
    out.close()
    return output


# main func
def main():
    # Remember to add a choice function regarding the removal of current folders
    if os.path.exists('traj.conf'):
        reset, mode, atom = read_conf_file()
        if reset == 'y':
            mkdir()
    else:
        judge = input('Do you want to start analyzing from the very beginning? (y/n)\n'
                      'Type y to remove all analysis folders and n to keep the current folder'
                      ' (e.g. reorder, etc.) for analysis: ')
        if judge == 'y':
            mkdir()
        mode = get_mode()
        atom = get_atomindex(mode)
    print("Run...")
    # The attribute of Trajectories Class involves
    #    self.Get_distance(1)
    #    self.TS_finder()
    #    self.Rearrangement()
    for filename in glob.glob('./ntraj/*.'):
        Trajectories(filename, atom, mode).TS_finder()
        Trajectories(filename, atom, mode).rearrangement()
    total, X, Y, revert, inter = 0, 0, 0, 0, 0
    for filename in glob.glob('./reorder/*.'):
        result = Trajectories(filename, atom, mode).classification()
        total += 1
        if result == 'X':
            X += 1
        elif result == 'Y':
            Y += 1
        elif result == 'R':
            revert += 1
        else:
            inter += 1
    print('Trajectory analysis complete!')
    if X + Y == 0:
        print('Neither product X nor Y was formed')
    else:
        output = log_results(X, Y, revert, inter, total)
        print(output)


if __name__ == '__main__':
    main()