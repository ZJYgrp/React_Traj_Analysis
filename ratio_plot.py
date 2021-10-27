import XYZtraj_Trajs as xyz
from XYZtraj_Trajs import Trajectories
import os
import glob

def log_results(A, B, revert, inter, total):
  out = open('./trajTS/traj_log', 'w+')
  out.write('Results\nTotal number of trajectories: '+str(total)+'\nTotal forming product: '+str(A+B)+'\nA: '+str(A)+' B: '+str(B)+' Reactant: '+str(revert)+' Intermediate: '+str(inter)+'\nPercent product A: '+str((A*100/(A+B)))+'%\nPercent product B: '+str(B*100/(A+B))+'%')


# main func
def main():
# Remember to add a choice function regarding the removal of current folders
    if os.path.exists('traj.conf'):
        reset, mode, atom = xyz.read_conf_file()
        if reset == 'y': xyz.mkdir()
    else:
        judge = input('Do you want to start analyzing from the very beginning? (y/n) Type y to remove all analysis folders and n to keep the current folder (e.g. reorder, etc.) for analysis: ')
        if judge == 'y': xyz.mkdir()
        mode = xyz.Get_mode()
        atom = xyz.Get_atomindex(mode)
# The attribute of Trajectories Class involves
#        self.Get_distance(1)
#        self.TS_finder()
#        self.Rearrangement()
    for filename in glob.glob('./ntraj/*.xyz'):
        Trajectories(filename, atom, mode).TS_finder()
        Trajectories(filename, atom, mode).Rearrangement()
#       self.Classification() Note: Classification has to be called after Rearrangement or when reorder has been performed.
    total, A, B, revert, inter = 0, 0, 0, 0, 0
    for filename in glob.glob('./reorder/*.xyz'):
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
