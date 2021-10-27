from Trajectories import Trajectories
import os
import glob

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
def log_results(A, B, revert, inter, total):
  out = open('./trajTS/traj_log', 'w+')
  out.write('Results\nTotal number of trajectories: '+str(total)+'\nTotal forming product: '+str(A+B)+'\nA: '+str(A)+' B: '+str(B)+' Reactant: '+str(revert)+' Intermediate: '+str(inter)+'\nPercent product A: '+str((A*100/(A+B)))+'%\nPercent product B: '+str(B*100/(A+B))+'%')


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
    for filename in glob.glob('./ntraj/*.):
        Trajectories(filename, atom, mode).TS_finder()
        Trajectories(filename, atom, mode).Rearrangement()
#       self.Classification() Note: Classification has to be called after Rearrangement or when reorder has been performed.
    total, A, B, revert, inter = 0, 0, 0, 0, 0
    for filename in glob.glob('./reorder/*.):
        result = Trajectories(filename, atom, mode).Classification()
        total += 1
        if result == 'X':
            A += 1
        elif result == 'Y':
            B += 1
        elif result == 'R':
            revert += 1
        else:
            inter += 1
    print('Trajectory analysis complete!')
    if (A+B == 0):
        print('Neither product A nor B was formed')
    else:
        log_results(A, B, revert, inter, total)
        print('Results\nTotal number of trajectories: '+str(total)+'\nTotal forming product: '+str(A+B)+'\nA: '+str(A)+' B: '+str(B)+' Reactant: '+str(revert)+' Intermediate: '+str(inter)+'\nPercent product A: '+str((A*100/(A+B)))+'%\nPercent product B: '+str(B*100/(A+B))+'%')
#    TS = hist('./trajTS/trajTs.txt')
#    print('Done histogram :', TS)
#    print('Complete !')


if __name__ == '__main__':
    main()
