import pandas as pd
import numpy as np
from matplotlib.pyplot import MultipleLocator
import matplotlib.pyplot as plt
import glob
class Trajectories_txt:
    def __init__(self,filename):
        time = []
        bond_1 = []
        bond_2 = []
        bond_3 = []
        with open(filename,'r') as fr:
            line_list = fr.readlines()
            for line in line_list:
                words = line.split(',')
                time.append(words[0])
                bond_1.append(float(words[1]))
                bond_2.append(float(words[2]))
                bond_3.append(float(words[3]))
        data = {"time":time,"bond_1":bond_1,"bond_2":bond_2,"bond_3":bond_3}
        f = pd.DataFrame(data,columns=["time","bond_1","bond_2","bond_3"])
        data = np.array(f)
        #array：time
        self.time_data = data[:, 0]
        #arry:bond_1 & bond2 &bond3
        self.bond_data1 = data[:, 1]
        self.bond_data2 = data[:, 2]
        self.bond_data3 = data[:, 3]
    def traj_plot(self,file):
        plt.plot(self.time_data,self.bond_data1,label = file )
        #plt.plot(self.time_data,self.bond_data2,label = file)
        #plt.plot(self.time_data,self.bond_data3,label = file)
        plt.title('Time-Bond',fontsize=18)
        plt.tick_params(axis='both',which='major',labelsize=12)
        plt.xlabel('Time(fs)',fontsize=14)
        plt.ylabel('Bond_length(Å)',fontsize=14)
        plt.legend()
        x_major_locator=MultipleLocator(100)
        y_major_locator=MultipleLocator(1)
        ax=plt.gca()
        ax.xaxis.set_major_locator(x_major_locator)
        ax.yaxis.set_major_locator(y_major_locator)
        
def main():  
    path = './TDD/'  
    for filename in glob.glob(path +'\\' '*.xyz.txt'):
        file_a = filename.split('\\')
        file = file_a[1]
        T = Trajectories_txt(filename)
        T.traj_plot(file)
    plt.show()
if __name__ =="__main__":
    main()
