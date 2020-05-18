# This script should be revised to an object-oriented program, which defines TS extraction class, trajectories classfication class, etc.
import numpy as np
import sys
import os
import glob
import Trajectories_txt as TX
# The attribute of Trajectories Class involves

# main func
def main():
    TX.Trajectories_txt('./trajTS/','').TS_plot(custom_bin_num=0,custom_range=[])
#    TX.Trajectories_txt('','./TDD_r2p/',N_col=3).Traj_plot()
#    for filename in glob.glob('./TDD_r2p/*.txt'):
#        TX.Trajectories_txt(filename).Classification()
    print('Done trajectory plottnig!')

if __name__ == "__main__":
    main()
