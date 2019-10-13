import numpy as np
import h5py
import os
import illustris_python as il
import matplotlib.pyplot as plt

def SelectDisk(snap_num):
    '''
    Input Snapshot number, like snap_num = 99 (z=0)
    Select disk galaxies, return haloID of them.
    '''
    #select halo stellar particles > 40000
    with h5py.File('/Raid0/zhouzb/TNGdata/offsets_0%d.hdf5'%snap_num,'r') as offset:
        haloSBT = (np.array(offset['Subhalo']['SnapByType']))[:,4]
        
    halolen = haloSBT[1:]
    halolen = np.append(halolen, halolen.max()) - haloSBT
    

    with h5py.File('/Raid0/zhouzb/TNGdata/stellar_circs.hdf5','r') as cir:
        haloID = np.array(cir['Snapshot_%d'%snap_num]['SubfindID'])
        cir07frac = np.array(cir['Snapshot_%d'%snap_num]['CircAbove07Frac'])
        MassTensor = np.array(cir['Snapshot_%d'%snap_num]['MassTensorEigenVals'])
    #circularity parameter ϵ > 0.2
    cir_mask = cir07frac > 0.2
    
    #flatness of the galaxy is defined as the ratio M1=(M2M3)**0.5 , disk galaxy's flatness < 0.7
    MassTensor = MassTensor[cir_mask]
    haloID = haloID[cir_mask]

    flat = MassTensor[:,0]/(MassTensor[:,1]*MassTensor[:,2])**0.5
    flat_mask = flat < 0.7

    haloID = haloID[flat_mask]

    mas_mask = halolen[haloID] >= 20000
    haloID = haloID[mas_mask]
    return haloID

snap_num = 99
ids = SelectDisk(99)
StellarMass = il.groupcat.loadSubhalos('f:/Linux/data/TNG/Groupcatalog/', 99, 'SubhaloMassInRadType')[:,4]
Gas = il.groupcat.loadSubhalos('f:/Linux/data/TNG/Groupcatalog/', 99, 'SubhaloMassInRadType')[:,0]
#Disk halo's mass
Diskmass = StellarMass[ids]
Diskmass = np.log10(Diskmass*10**10)
DiskGas = Gas[ids]
#Create figer
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_xlabel('Stellar Mass')
ax1.set_ylabel('Halo number N')
ax2 = ax1.twinx()
ax2.set_ylabel('Bar Fraction')

#load disk halo's A2 file
A2list = []
for haloID in ids:
    tmp = np.load('f:/Linux/data/TNG/A2_snap99/%d.npy'%haloID)
    A2 = np.array(tmp[0])
    A2 = A2[int(len(A2) / 100):]
    A2list.append(A2.max())

#plot histogram
n,bins,others = ax1.hist(Diskmass, 20, rwidth=0.9)
# fig.xlim(10.0,12)

A2list = np.array(A2list)

for i in range(20):
    low = bins[i]
    high = bins[i+1]

    mask = (Diskmass >= low) & (Diskmass < high)
    A2tmp = A2list[mask]
    if len(A2tmp) == 0:
        continue 
    Barfraction = len(A2tmp[A2tmp > 0.15]) / len(A2tmp)
    ax2.plot(low, Barfraction, 'o', c = 'r')






