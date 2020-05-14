import numpy as np
import sys
import os
import glob

class Trajectories_xyz:
    def __init__(self, file, atoms):
        # trajectory format for ProgDyn output
        self.name = os.path.basename(file)
        print ('Working on '+self.name)
        if '.xyz' not in file:
            raise TypeError('A ProgDyn .xyz file must be provided')
        if os.stat(file).st_size == 0:
            print('The file '+self.name+' is empty.')
            return
        #Creating new folders for the following analysis
        #Open new file handles and get parameters
        self.atoms = atoms
        self.lines = open(file).readlines()
        self.n_lines = len(self.lines)
        if self.n_lines == 0:
            print('The file ' + self.name + ' is empty.')
            return
        self.n_atoms = int(self.lines[0].split()[0])
        self.n_idx = int((self.n_lines / (self.n_atoms + 2)))

    # This function is used to get distance from coordinates
    def Get_distance(self,n):
        if not hasattr(self, 'atoms'):
            print('The xyz file ' + self.name + ' has not been successfully initiated.')
            return
        elif self.atoms == 0:
            print('The xyz file ' + self.name + ' has zero atoms.')
            return
        X = np.zeros(len(self.atoms))
        Y = np.zeros(len(self.atoms))
        Z = np.zeros(len(self.atoms))
        Bonds = np.zeros(int(len(self.atoms)/2))
        for i in range(0, len(self.atoms)):
            X[i] = float(self.lines[(self.n_atoms + 2) * n + self.atoms[i] + 1].split()[1])
            Y[i] = float(self.lines[(self.n_atoms + 2) * n + self.atoms[i] + 1].split()[2])
            Z[i] = float(self.lines[(self.n_atoms + 2) * n + self.atoms[i] + 1].split()[3])
        for j in range(0, len(Bonds)):
            Bonds[j] = round(((X[j*2+1] - X[j*2]) ** 2 + (Y[j*2+1] - Y[j*2]) ** 2 + (Z[j*2+1] - Z[j*2]) ** 2) ** .5, 3)
        return Bonds

    # TS finder is used to collect the sampled TS geometries from trajectories. The sampled TS is usually the starting point of each trajectory
    def TS_finder(self):
        if not hasattr(self, 'n_idx'):
            print('The xyz file ' + self.name + ' has not been successfully initiated.')
            return
        elif self.n_idx == 0:
            print('The xyz file ' + self.name + ' has zero snapshots.')
            return
        fileout_TS_xyz = open('./trajTS/trajTs.xyz', 'a')
        fileout_TS = open('./trajTS/trajTs.txt', 'a')
        for i in range(0, self.n_idx):
            if len(self.lines[1].split()) < 7:
                print('The xyz file ' + self.name + ' does not have the snapshot numeration on the 7th word of the title line.')
                break
            elif int(self.lines[1 + i * (self.n_atoms + 2)].split()[6]) == 1:
                bond_TS= self.Get_distance(i)
                fileout_TS.write(self.name + ', ')
                for j in range(0, len(bond_TS)):
                    fileout_TS.write(str(bond_TS[j]) + ', ')
                fileout_TS.write("\n")
                fileout_TS.close()
                for i in range(0, self.n_atoms + 2):
                    fileout_TS_xyz.write(self.lines[i])
                fileout_TS_xyz.close()
                break
            else:
                print('The xyz file ' + self.name + ' does not have the TS geometry!')

    def Rearrangement(self):
        if not hasattr(self, 'lines'):
            print('The xyz file ' + self.name + ' has not been successfully initiated.')
            return
        elif len(self.lines) == 0:
            print('The xyz file ' + self.name + ' has zero lines.')
            return
        fileout_reorder = open('./reorder/' + self.name, 'w')
        if len(self.lines[1].split()) < 7:
            print('The xyz file ' + self.name + ' does not have the snapshot numeration on the 7th word of the title line.')
            return
        elif int(self.lines[1].split()[6]) != 1:
            print('I cannot find the first TS and reorder is not feasible; break!')
            return
        else:
            for i in range(1, self.n_idx):
                if self.lines[1].split()[6] == self.lines[1 + i * (self.n_atoms + 2)].split()[6]: break
            n1 = i
            n2 = self.n_idx - i
            if n1 == n2 == 1:
                print('The file ' + self.name + ' only has two TS points.')
                return
            else:
                bond_TS = self.Get_distance(0)
                bond_D1 = self.Get_distance(n1 - 1)
                bond_D2 = self.Get_distance(self.n_idx - 1)
                print('Bond 1 changes from D1:', bond_D1[0], ' to TS:', bond_TS[0], ' then to D2:',
              bond_D2[0])
                print('Assuming bond 1 forms from R to P')
                if (bond_D2[0] > bond_D1[0]):
                    for i in range(0, n2):
                        for j in range(0, self.n_atoms + 2):
                            fileout_reorder.write(self.lines[(self.n_idx - 1 - i) * (self.n_atoms + 2) + j])
                    for i in range(0, n1):
                        for j in range(0, self.n_atoms + 2):
                            fileout_reorder.write(self.lines[i * (self.n_atoms + 2) + j])
                if (bond_D1[0] > bond_D2[0]):
                    for i in range(0, n1):
                        for j in range(0, self.n_atoms + 2):
                            fileout_reorder.write(self.lines[(n1 - 1 - i) * (self.n_atoms + 2) + j])
                    for i in range(0, n2):
                        for j in range(0, self.n_atoms + 2):
                            fileout_reorder.write(self.lines[(i + n1) * (self.n_atoms + 2) + j])
                fileout_reorder.close()
## classification function take reordered trajectories to prcess, generating distance/angle/dihedral time series that inform where the trajectories come from and end up.
    def Classification(self):
        if not hasattr(self, 'name'):
            print('The xyz file ' + self.name + ' has not been successfully initiated.')
            return
        elif self.n_idx == 0:
            print('The xyz file ' + self.name + ' has zero snapshots.')
            return
        fileout_traj = open('./TDD/' + self.name + '.txt', 'w')
        for i in range(0, self.n_idx):
            if int(self.lines[1 + i * (self.n_atoms + 2)].split()[6]) == 1: break
        n1=i
        bond_R = self.Get_distance(0)
        bond_TS = self.Get_distance(n1)
        bond_P = self.Get_distance(self.n_idx-1)
# now writing every snapshots to TDD
        for i in range(0,self.n_idx):
            runpoint = int(self.lines[1 + i * (self.n_atoms + 2)].split()[6])
            bond = self.Get_distance(i)
            if i<n1:
                fileout_traj.write(str(-runpoint+1)+ ' , ')
                for j in range(0, len(bond)):
                    fileout_traj.write(str(bond[j]) + ', ')
                fileout_traj.write("\n")
            elif i>n1:
                fileout_traj.write(str(runpoint-1) + ' , ')
                for j in range(0, len(bond)):
                    fileout_traj.write(str(bond[j]) + ', ')
                fileout_traj.write("\n")
        fileout_traj.close()
#Now start classifying trajectories
        if (bond_R[0] > bond_TS[0] > bond_P[0]):
            os.system('cp ./TDD/' + self.name + '.txt '+'./TDD_r2p/' + self.name + '.txt')
            print('go to r2p')
        elif (bond_R[0] >= bond_TS[0]) and (bond_P[0] >= bond_TS[0]):
            os.system('cp ./TDD/' + self.name + '.txt '+'./TDD_r2r/' + self.name + '.txt')
            print('go to r2r')
        elif (bond_R[0] <= bond_TS[0]) and (bond_P[0] <= bond_TS[0]):
            os.system('cp ./TDD/' + self.name + '.txt '+'./TDD_p2p/' + self.name + '.txt')
            print('go to p2p')
        else:
            os.system('cp ./TDD/' + self.name + '.txt '+'./TDD_inter/' + self.name + '.txt')
            print('go to intermediate')