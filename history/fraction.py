import h5py
import numpy as np
import matplotlib.pyplot as plt

def finddisk_135():
    '''Return a list of all Disk galaxies in this snapshot'''
    with h5py.File('/Raid0/zhouzb/stellar_circs.hdf5','r') as cir:
        return np.array(cir['Snapshot_135']['SubfindID'])

def groups_135(data):
    '''
    parameter: data= 'len','cm','half_r','rad_r'. 
    'len' return SubhaloID corresponded stellar numbers, shape:(halonumber+1,). 
    'cm' return Subhalo's center of mass's 3D-coordinates, shape:(halonumber, 3).
    'half_r' return Sum of masses of all particles/cells (split by type) within the stellar half mass radius, shape(halonumber, 6) .
    'rad_r' return Sum of masses of all particles/cells (split by type) within twice the stellar half mass radius, shape(halonumber, 6)
    'max_r' retrun Sum of masses of all particles/cells (split by type) within the radius of Vmax.
    '''
    if data == 'len' :
        with h5py.File('/Raid1/Illustris/Illustris-1/Groupcat/groups_135.0.hdf5','r') as f:
            Subhalolen = np.array(f['Offsets']['Subhalo_SnapByType'][:,4])
        for n in range(1,8):
            with h5py.File('/Raid1/Illustris/Illustris-1/Groupcat/groups_135.%d.hdf5'%n,'r') as f:
                Subhalolen = np.concatenate((Subhalolen, np.array(f['Offsets']['Subhalo_SnapByType'][:,4])))
        return Subhalolen

    elif data == 'cm':
        with h5py.File('/Raid1/Illustris/Illustris-1/Groupcat/groups_135.0.hdf5','r') as f:
            SubhaloCM_list = np.array(f['Subhalo']['SubhaloCM'])
        for n in range(1,8):
            with h5py.File('/Raid1/Illustris/Illustris-1/Groupcat/groups_135.%d.hdf5'%n,'r') as f:
                SubhaloCM_list = np.concatenate((SubhaloCM_list, np.array(f['Subhalo']['SubhaloCM'])))
        return SubhaloCM_list
    
    elif data == 'half_r':
        with h5py.File('/Raid1/Illustris/Illustris-1/Groupcat/groups_135.0.hdf5','r') as f:
            half_r = np.array(f['Subhalo']['SubhaloHalfmassRadType'][:,4])
        for n in range(1,8):
            with h5py.File('/Raid1/Illustris/Illustris-1/Groupcat/groups_135.%d.hdf5'%n,'r') as f:
                half_r = np.concatenate((half_r, np.array(f['Subhalo']['SubhaloHalfmassRadType'][:,4])))
        return half_r



halo_len = groups_135('len')
half_r = groups_135('half_r')


haloID = finddisk_135()
#a2>0.2 halos number
bar02 = 0
allhalo = 0
for i in haloID:
    if allhalo%100 == 0:
        print(allhalo)
    try:
        temp = np.load('/Raid0/zhouzb/A2_all/%d.npy'%i,'r')
    except FileNotFoundError:
        continue
    else:
        allhalo += 1
        rlist = np.array(temp[1])
        rlist.sort()
        a2 = np.array(temp[0])
        sel = (rlist >= (half_r[i] / 10)) & (rlist <= 2 * half_r[i])
        
        plt.plot(rlist[sel], a2[sel])
        plt.ylim(0,0.7)
        plt.plot([half_r[i], half_r[i]],[0, 0.6])
        plt.savefig('/Raid0/zhouzb/halo_40000/%d.png'%i)
        plt.close()
        