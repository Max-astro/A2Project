import numpy as np
import h5py
import sys
sys.path.append('F:\Linux')
import illustris_python as il
from sklearn.neighbors import KDTree
from matplotlib import ticker
import matplotlib.pyplot as plt


cm_TNG = il.groupcat.loadHalos('f:/Linux/data/TNG/TNG-100',99,'GroupPos')
cm_il1 = il.groupcat.loadHalos('f:/Linux/data/illustris-1',135,'GroupPos')

velo_TNG = il.groupcat.loadHalos('f:/Linux/data/TNG/TNG-100',99,'GroupVel')
velo_il1 = il.groupcat.loadHalos('f:/Linux/data/illustris-1',135,'GroupVel')

mass_TNG = il.groupcat.loadHalos('f:/Linux/data/TNG/TNG-100',99,'GroupMassType')[:,1]
mass_il1 = il.groupcat.loadHalos('f:/Linux/data/illustris-1',135,'GroupMassType')[:,1]

tree = KDTree(cm_TNG, leaf_size = 2)
#dis, index = tree.query(cm_il1, k=3)

#if halo in il1 have the similar coordinates and velosity, put haloID into samehalo[]
#samehalo[ [il1_haloID], [TNG_haloID] ]
sameHalo = [[],[]]
count_coor = 0
count = 0
for haloID in range(len(cm_il1)):
    dis, index = tree.query(cm_il1[haloID].reshape(1,-1), k=3)
    for i in range(3):
        if dis[0,i] < 75:
            count_coor += 1
            vel = velo_il1[haloID] - velo_TNG[index[0,i]]
            if (vel < velo_il1[haloID] / 10).all():
                sameHalo[0].append(haloID)
                sameHalo[1].append(index[0,i])
                count += 1
                print('find one', count, count_coor)

np.save('/Raid0/zhouzb/samehalo.npy',sameHalo)
print('All done')

#
mass_TNG = np.log10(mass_TNG*10**10)
mass_il1 = np.log10(mass_il1*10**10)
mass_TNG[np.isinf(mass_TNG)] = 0
mass_il1[np.isinf(mass_il1)] = 0

#Unique id list
utng , ind = np.unique(sameHalo[1], return_index=True)
uil1 = sameHalo[0][ind]
utng = sameHalo[1][ind]

#Plot mass function
fig = plt.figure()
ax = fig.add_subplot(111)
n,bins,others = ax.hist(mass_m200_TNG[utng], 40, rwidth=0.9)
ax.set_xlim(7,13)
ax.set_ylim(0,150000)
ax.set_xlabel('Halo Mass( log(M_sun) )')
ax.set_ylabel('Halo number')
plt.title('Illustris-TNG')
formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True) 
formatter.set_powerlimits((-1,1)) 
ax.yaxis.set_major_formatter(formatter)  
plt.savefig('f:/Linux/data/result/TNG_match_totalmass.png',dpi=300)