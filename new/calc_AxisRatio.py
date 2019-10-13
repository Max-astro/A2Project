import numpy as np
import h5py
import json
import sys
sys.path.append('F:\Linux')
import illustris_python as il
from datetime import datetime


def LoadMergHist(simu, subhaloID):
    '''
    return subhalo's main progenitor and merger history with snapshot
    '''
    if simu == 'TNG':
        ldir = 'f:/Linux/localRUN/tng_DiskMerTree/%d.json' % subhaloID
    else:
        ldir = 'f:/Linux/localRUN/il1_DiskMerTree/%d.json' % subhaloID
    
    with open(ldir) as f:
        data = json.load(f)
    
    Main = np.array(data['Main'])
    return dict(zip(Main[:, 0], Main[:, 1])), np.array(data['Mergers'])

def MassTensorEigenVals(coor, mas, half_r):
    '''
    Return eigenvalues of the mass tensor, sorted by M1 < M2 < M3
    '''
    r = coor - coor[0]
    r[r > 37500] -= 75000
    r[r < -37500] += 75000
    dis = np.linalg.norm(r, axis=1)
    inside = dis < (half_r * 2)
    r = r[inside]
    mas = mas[inside]

    M_x = ((mas * (r[:, 0]/0.6774)**2).sum())**0.5 / mas.sum()**0.5
    M_y = ((mas * (r[:, 1]/0.6774)**2).sum())**0.5 / mas.sum()**0.5
    M_z = ((mas * (r[:, 2]/0.6774)**2).sum())**0.5 / mas.sum()**0.5
    M = np.array([M_x, M_y, M_z])
    M.sort()
    return M



diskID = np.load('f:/Linux/localRUN/diskID_il1.npy')

'''
Illustris-1 Snapshot-Redshift:
snap_127 z=0.1
snap_120 z=0.2
snap_113 z=0.3
snap_108 z=0.4
snap_103 z=0.5
snap_95 z=0.7
snap_85 z=1.0
snap_75 z=1.5
snap_68 z=2.0
'''
'''
SnapList = [135, 127, 120, 113, 108, 103, 95, 85, 75, 68]
RedShift = [0, 0.1 , 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0]
'''

'''
TNG Snapshot-Redshift:
snap_91 z=0.1
snap_84 z=0.2
snap_78 z=0.3
snap_72 z=0.4
snap_67 z=0.5
snap_59 z=0.7
snap_50 z=1.0
snap_40 z=1.5
snap_33 z=2.0
'''

SnapList = [135, 127, 120, 113, 108, 103, 95, 85, 75, 68]
RedShift = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0]

# move = np.array([0, 1, -1, 2, -2])

MTE = {}
for haloID in diskID:
    a = datetime.now()
    MTE[haloID] = []
    proj = LoadMergHist('il1', haloID)[0]

    for snap in SnapList[1:]:
        try:
            hid = proj[snap]
        except KeyError:
            MTE[haloID].append(-1)
        else:
            print('hid ', hid, 'snap ', snap)
            if hid == -1:
                MTE[haloID].append(-1)
            else:    
                with h5py.File('f:/Linux/il1_cutoff_dm/snap_%d/cutout_%d.hdf5' % (snap, hid)) as f:
                    coor = np.array(f['PartType1']['Coordinates'])
                    DMhalf_r = il.func.loadSubhalos('il1', snap, 'SubhaloHalfmassRad')[hid]
                MTE[haloID].append(MassTensorEigenVals(coor, np.ones(len(coor)), DMhalf_r))
    b = datetime.now()
    print('halo_%d finished. Time: '%haloID,(b-a).seconds)
print('All done')
np.save('f:/Linux/localRUN/MTE_il1.npy', MTE)