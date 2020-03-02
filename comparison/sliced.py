import numpy as np
import h5py
import illustris_python as il
import matplotlib.pyplot as plt


#TNG slice
# coor = il.snapshot.loadSubset('/Raid1/Illustris/TNG-100',99, 1, fields='Coordinates')
#runtime test

for add in range(10000, 70000, 10000):
    for i in range(512):
        with h5py.File('/Raid1/Illustris/TNG-100/Snapshot/snap_99/snap_099.%d.hdf5'%i,'r') as f:
            coor = np.array(f['PartType0']['Coordinates'])
            select = (coor[:,0] > (add - 5)) & (coor[:,0] <= add) 
            sliced = (coor[select])
            ss = np.zeros(len(sliced[:,0]))
            for i in range(len(ss)):
                ss[i] = 0.1

            plt.scatter(sliced[:,1], sliced[:,2], c='c', marker = '.',s = ss)
            print("part %d finished"%i)



    plt.savefig('/Raid0/zhouzb/fig_TNG/sliced/TNG_%d.png'%add, dpi = 1000)
    plt.close()

