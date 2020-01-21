
import numpy as np
import h5py
import json
sys.path.append('F:\Linux')
import illustris_python as il
import matplotlib.pyplot as plt

def Flatness(MassTensor):
    #c / a = (M3)**0.5 / (M1)**0.5
    return np.sqrt(MassTensor[0]) / np.sqrt(MassTensor[2])

def BtoA(MassTensor):
    #b / a = (M2)**0.5 / (M1)**0.5
    return np.sqrt(MassTensor[0]) / np.sqrt(MassTensor[1])

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
barredID = np.load('f:/Linux/localRUN/barredID_il1.npy')

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


SnapList = [127, 120, 113, 108, 103, 95, 85, 75, 68]
RedShift = [0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0]

# move = np.array([0, 1, -1, 2, -2])


for haloID in barredID[::3]:
    proj = LoadMergHist('il1', haloID)[0]
    ca = []
    ba = []
    for snap in SnapList:
        try:
            hid = proj[snap]
        except KeyError:
            apdata = 0
        else:
            print('hid ', hid, 'snap ', snap)
            if hid == -1:
                apdata = 0
            else:    
                with h5py.File('g:/il1_cutoff/snap_%d/cutout_%d.hdf5' % (snap, hid)) as f:
                    coor = np.array(f['PartType1']['Coordinates'])
                    DMhalf_r = il.func.loadSubhalos('il1', snap, 'SubhaloHalfmassRadType')[hid, 1]
                MTE = MassTensorEigenVals(coor, np.ones(len(coor)), DMhalf_r)

                ba.append(BtoA(MTE))
                ca.append(Flatness(MTE))
                apdata = 1

        if apdata == 0:
            ba.append(1)
            ca.append(1)
    CA = []
    x = []
    for i in range(9):
        if ca[i] != 1:
            CA.append(ca[i])
            x.append(RedShift[i])

    plt.ylim(0.5,1)
    plt.plot(x, CA, color='black', linewidth=0.2)

    plt.plot(RedShift, ba)


import numpy as np
import h5py
import illustris_python as il
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

data = {}
barredID = np.load('/Raid0/zhouzb/TNG/barredID_4WP_TNG.npy')
for haloID in barredID:
    coor = il.func.loadgalaxy('TNG', 99, haloID, 1, 'Coordinates')
    r = il.func.loadSubhalos('TNG', 99, 'SubhaloHalfmassRadType')[haloID, 1]
    data[haloID] = MassTensorEigenVals(coor, np.ones(len(coor)), r)

np.save('/Raid0/zhouzb/barMTE_TNGsnap99.npy', data)