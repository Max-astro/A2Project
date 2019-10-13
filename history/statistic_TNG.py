import numpy as np
import h5py
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
    #Total halo number, it also be an index of haloID
    ids = np.arange(len(haloSBT))

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

    mas_mask = halolen[haloID] >= 40000
    haloID = haloID[mas_mask]
    return haloID

snap_num = 99

ids = SelectDisk(99)
StellarHalfmassRads = (il.groupcat.loadSubhalos('/Raid1/Illustris/TNG-100',snap_num,'SubhaloHalfmassRadType'))[:,4]

new_bar02 = 0
new_bar04 = 0

new_bar = [ [], [] ]
new_Strong = [ [], [] ]
new_disk = []
did = 0
for haloID in ids:
    did += 1
    if did % 10 ==0:
        print('%d halo finished'%did)
    tmp = np.load('/Raid0/zhouzb/TNG_a2/snap_%d/%d.npy'%(snap_num ,haloID),'r')
    A2 = np.array(tmp[0])
    rlist = np.array(tmp[1])
    half_r = StellarHalfmassRads[haloID]


    #Ignore the most center stellar particles (First 0.1%)
    A2 = A2[int(len(A2) / 100):]
    rlist = rlist[int(len(rlist) / 100):]
    new_disk.append(A2.max())
    if A2.max() > 0.2:
        new_bar[0].append(haloID)
        new_bar[1].append(A2.max())
        new_bar02 += 1
        if A2.max() > 0.4:
            new_Strong[0].append(haloID)
            new_Strong[1].append(A2.max())
            new_bar04 += 1


print('Disk halo number: %d'%len(ids))
print('A2 > 0.2 number: '+str(new_bar02))
print('A2 > 0.4 number: '+str(new_bar04))

