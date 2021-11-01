# TS Extraction and Trajectory classification class
import os
import matplotlib as mpl
import numpy as np

mpl.use('Agg')


class Trajectories:
    def __init__(self, file, atoms, mode):
        # trajectory format for ProgDyn output
        self.name = os.path.basename(file)
        # print ('Working on '+self.name)
        if '.xyz' not in file:
            raise TypeError('A ProgDyn .xyz file must be provided')
        if os.stat(file).st_size == 0:
            print('The file ' + self.name + ' is empty.')
            return
        # Creating new folders for the following analysis
        # Open new file handles and get parameters
        self.atoms = atoms
        self.mode = mode
        self.lines = open(file).readlines()
        self.n_lines = len(self.lines)
        if self.n_lines == 0:
            print('The file ' + self.name + ' is empty.')
            return
        self.n_atoms = int(self.lines[0].split()[0])
        self.n_idx = int((self.n_lines / (self.n_atoms + 2)))

    # This function is used to get distance from coordinates
    def get_distance(self, n):
        if not hasattr(self, 'atoms'):
            print('The xyz file ' + self.name + ' has not been successfully initiated.')
            return -1
        if self.atoms == 0:
            print('The xyz file ' + self.name + ' has zero atoms.')
            return -1
        X = np.zeros(len(self.atoms))
        Y = np.zeros(len(self.atoms))
        Z = np.zeros(len(self.atoms))
        Bonds = np.zeros(int(len(self.atoms) / 2))
        for i in range(0, len(self.atoms)):
            X[i] = float(self.lines[(self.n_atoms + 2) * n + self.atoms[i] + 1].split()[1])
            Y[i] = float(self.lines[(self.n_atoms + 2) * n + self.atoms[i] + 1].split()[2])
            Z[i] = float(self.lines[(self.n_atoms + 2) * n + self.atoms[i] + 1].split()[3])
        for j in range(0, len(Bonds)):
            x_partial = (X[j * 2 + 1] - X[j * 2]) ** 2
            y_partial = (Y[j * 2 + 1] - Y[j * 2]) ** 2
            z_partial = (Z[j * 2 + 1] - Z[j * 2]) ** 2
            dist = (x_partial+y_partial+z_partial) ** .5
            Bonds[j] = round(dist, 3)

        return Bonds

    # TS finder is used to collect the sampled TS geometries from trajectories.
    # The sampled TS is usually the starting point of each trajectory
    def TS_finder(self):
        if not hasattr(self, 'n_idx'):
            print('The xyz file ' + self.name + ' has not been successfully initiated.')
            return
        if self.n_idx == 0:
            print('The xyz file ' + self.name + ' has zero snapshots.')
            return
        fileout_TS_xyz = open('./trajTS/trajTs.xyz', 'a')
        fileout_TS = open('./trajTS/trajTs.txt', 'a')
        for i in range(0, self.n_idx):
            if len(self.lines[1].split()) < 7:
                print('The xyz file ' + self.name
                      + ' does not have the snapshot enumeration '
                        'on the 7th word of the title line.')
                break
            elif int(self.lines[1 + i * (self.n_atoms + 2)].split()[6]) == 1:
                bond_TS = self.get_distance(i)
                fileout_TS.write(self.name + ', ')
                fileout_TS.write(', '.join([str(bond_TS[j]) for j in range(len(bond_TS))]))
                # for j in range(0, len(bond_TS)):
                #    fileout_TS.write(str(bond_TS[j]) + ', ')
                fileout_TS.write('\n')
                fileout_TS.close()
                for j in range(0, self.n_atoms + 2):
                    fileout_TS_xyz.write(self.lines[j])
                fileout_TS_xyz.close()
                break
            else:
                print('The xyz file ' + self.name + ' does not have the TS geometry!')

    def rearrangement(self):
        if not hasattr(self, 'lines'):
            print('The xyz file ' + self.name + ' has not been successfully initiated.')
            return
        if self.lines == []:
            print('The xyz file ' + self.name + ' has zero lines.')
            return
        fileout_reorder = open('./reorder/' + self.name, 'w')
        if len(self.lines[1].split()) < 7:
            print('The xyz file ' + self.name
                  + ' does not have the snapshot numeration on the 7th word of the title line.')
            return
        if int(self.lines[1].split()[6]) != 1:
            print('I cannot find the first TS and reorder is not feasible; break!')
            return
        for i in range(1, self.n_idx):
            if self.lines[1].split()[6] == self.lines[1 + i * (self.n_atoms + 2)].split()[6]:
                break
        n1 = i
        n2 = self.n_idx - i
        if n1 == n2 == 1:
            print('The file ' + self.name + ' only has two TS points.')
            return

        bond_TS = self.get_distance(0)
        bond_D1 = self.get_distance(n1 - 1)
        bond_D2 = self.get_distance(self.n_idx - 1)
        # print('Bond 1 changes from D1:', bond_D1[0],
        # ' to TS:', bond_TS[0], ' then to D2:',bond_D2[0])
        # print('Assuming bond 1 forms from R to P')
        if bond_D2[0] > bond_D1[0]:
            for i in range(0, n2):
                for j in range(0, self.n_atoms + 2):
                    fileout_reorder.write(self.lines[(self.n_idx - 1 - i) * (self.n_atoms + 2) + j])
            for i in range(0, n1):
                for j in range(0, self.n_atoms + 2):
                    fileout_reorder.write(self.lines[i * (self.n_atoms + 2) + j])
        if bond_D1[0] > bond_D2[0]:
            for i in range(0, n1):
                for j in range(0, self.n_atoms + 2):
                    fileout_reorder.write(self.lines[(n1 - 1 - i) * (self.n_atoms + 2) + j])
            for i in range(0, n2):
                for j in range(0, self.n_atoms + 2):
                    fileout_reorder.write(self.lines[(i + n1) * (self.n_atoms + 2) + j])
            fileout_reorder.close()

    # classification function take reordered trajectories to process, generating distance/angle/dihedral time series
    # that inform where the trajectories come from and end up.
    def classification(self):
        if not hasattr(self, 'name'):
            print('The xyz file ' + self.name + ' has not been successfully initiated.')
            return
        if self.n_idx == 0:
            print('The xyz file ' + self.name + ' has zero snapshots.')
            return
        fileout_traj = open('./TDD/' + self.name + '.txt', 'w')
        for i in range(0, self.n_idx):
            if int(self.lines[1 + i * (self.n_atoms + 2)].split()[6]) == 1:
                break
        n1 = i
        bond_R = self.get_distance(0)
        bond_TS = self.get_distance(n1)
        bond_P = self.get_distance(self.n_idx - 1)
        # now writing every snapshots to TDD
        for i in range(0, self.n_idx):
            runpoint = int(self.lines[1 + i * (self.n_atoms + 2)].split()[6])
            bond = self.get_distance(i)
            if i < n1:
                fileout_traj.write(str(-runpoint + 1) + ', ')
                fileout_traj.write(', '.join([str(bond[j]) for j in range(len(bond))]))
                # for j in range(0, len(bond)):
                #    fileout_traj.write(str(bond[j]) + ', ')
                fileout_traj.write('\n')
            elif i > n1:
                fileout_traj.write(str(runpoint - 1) + ', ')
                fileout_traj.write(', '.join([str(bond[j]) for j in range(len(bond))]))
                # for j in range(0, len(bond)):
                #    fileout_traj.write(str(bond[j]) + ', ')
                fileout_traj.write('\n')
        fileout_traj.close()

        # Now start classifying trajectories

        if self.mode == 1:
            if bond_R[0] > bond_TS[0] > bond_P[0]:
                if bond_P[1] < bond_P[2]:
                    os.system('cp ./TDD/' + self.name + '.txt ' + './TDD_r2pX/' + self.name + '.txt')
                    # print('go to r2pX')
                    return 'XR'

                os.system('cp ./TDD/' + self.name + '.txt ' + './TDD_r2pY/' + self.name + '.txt')
                # print('go to r2pY')
                return 'YR'

            if bond_R[0] >= bond_TS[0] and bond_P[0] >= bond_TS[0]:
                os.system('cp ./TDD/' + self.name + '.txt ' + './TDD_r2r/' + self.name + '.txt')
                # print('go to r2r')
                return 'R'
            if bond_R[0] <= bond_TS[0] and bond_P[0] <= bond_TS[0]:
                if bond_P[1] < bond_P[2]:
                    os.system('cp ./TDD/' + self.name + '.txt ' + './TDD_p2pX/' + self.name + '.txt')
                    # print('go to p2pX')
                    return 'XP'

                os.system('cp ./TDD/' + self.name + '.txt ' + './TDD_p2pY/' + self.name + '.txt')
                # print('go to p2pY')
                return 'YP'

            os.system('cp ./TDD/' + self.name + '.txt ' + './TDD_inter/' + self.name + '.txt')
            print('go to intermediate')
            return 'I'
        if self.mode == 2:
            if bond_R[0] > bond_TS[0] > bond_P[0] or bond_R[1] > bond_TS[1] > bond_P[1]:
                if bond_P[0] < bond_P[1]:
                    os.system('cp ./TDD/' + self.name + '.txt ' + './TDD_r2pX/' + self.name + '.txt')
                    print('go to r2pX')
                    return 'XR'

                os.system('cp ./TDD/' + self.name + '.txt ' + './TDD_r2pY/' + self.name + '.txt')
                print('go to r2pY')
                return 'YR'

            if bond_R[0] >= bond_TS[0] and bond_P[0] >= bond_TS[0] and bond_R[1] >= bond_TS[1] and bond_P[1] >= \
                    bond_TS[1]:
                os.system('cp ./TDD/' + self.name + '.txt ' + './TDD_r2r/' + self.name + '.txt')
                print('go to r2r')
                return 'R'
            if bond_R[0] <= bond_TS[0] and bond_P[0] <= bond_TS[0] or bond_R[1] <= bond_TS[1] and bond_P[1] <= \
                    bond_TS[1]:
                if bond_P[0] < bond_P[1]:
                    os.system('cp ./TDD/' + self.name + '.txt ' + './TDD_p2pX/' + self.name + '.txt')
                    print('go to p2pX')
                    return 'XP'

                os.system('cp ./TDD/' + self.name + '.txt ' + './TDD_p2pY/' + self.name + '.txt')
                print('go to p2pY')
                return 'YP'

            os.system('cp ./TDD/' + self.name + '.txt ' + './TDD_inter/' + self.name + '.txt')
            print('go to intermediate')
            return 'I'

    def formation_time(self):
        if not hasattr(self, 'name'):
            print('The xyz file ' + self.name + ' has not been successfully initiated.')
            return -1
        if self.n_idx == 0:
            print('The xyz file ' + self.name + ' has zero snapshots.')
            return -1
        fileout_traj = open('./TDD/' + self.name + '.txt', 'w')
        for i in range(0, self.n_idx):
            if int(self.lines[1 + i * (self.n_atoms + 2)].split()[6]) == 1:
                break
        n1 = i

        bond1, bond2, bond3 = -1, -1, -1
        product = ''
        # now writing every snapshots to TDD
        for i in range(n1, self.n_idx):
            runpoint = int(self.lines[1 + i * (self.n_atoms + 2)].split()[6])
            bond = self.get_distance(i)
            if self.mode == 1:
                if bond1 == -1 and bond[0] <= 1.6:
                    bond1 = runpoint - 1
                if product == '' and bond[1] <= 1.6:
                    bond2 = runpoint - 1
                    product = 'X'
                if product == '' and bond[2] <= 1.6:
                    bond3 = runpoint - 1
                    product = 'Y'
            elif self.mode == 2:
                if product == '' and bond[0] <= 1.6:
                    bond1 = runpoint - 1
                    product = 'X'
                if product == '' and bond[1] <= 1.6:
                    bond2 = runpoint - 1
                    product = 'Y'
            if i < n1:
                fileout_traj.write(str(-runpoint + 1) + ', ')
                fileout_traj.write(', '.join([str(bond[j]) for j in range(len(bond))]))
                # for j in range(0, len(bond)):
                #    fileout_traj.write(str(bond[j]) + ', ')
                fileout_traj.write('\n')
            elif i > n1:
                fileout_traj.write(str(runpoint - 1) + ', ')
                fileout_traj.write(', '.join([str(bond[j]) for j in range(len(bond))]))
                # for j in range(0, len(bond)):
                #    fileout_traj.write(str(bond[j]) + ', ')
                fileout_traj.write('\n')
        fileout_traj.close()
        if self.mode == 1:
            if bond1 < -1 or bond2 < -1 or bond3 < -1:
                return -1
            if product == 'X':
                return [product, bond1, bond2]
            if product == 'Y':
                return [product, bond1, bond3]
            return -1
        if self.mode == 2:
            if bond1 < -1 or bond2 < -1:
                return -1
            if product == 'X':
                return [product, bond1]
            if product == 'Y':
                return [product, bond2]
            return -1
