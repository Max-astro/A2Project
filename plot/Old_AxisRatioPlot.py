import numpy as np
import h5py
import json
import sys
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


# path_99 = 'f:/Linux/data/TNG/cutoff/disk_99'
diskID = np.load('f:/Linux/localRUN/diskID_4WP.npy')
barredID = np.load('f:/Linux/localRUN/barredID_4WP_TNG.npy')
MTE = np.load('f:/Linux/localRUN/MTE_TNGdisk_4WP.npy').item()
A2list = np.load('f:/Linux/localRUN/TNG_A2withRedshift.npy').item()
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

SnapList = [99, 91, 84, 78, 72, 67, 59, 50, 40, 33]
RedShift = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0]

# move = np.array([0, 1, -1, 2, -2])

nobar = []
for i in diskID:
    if i not in barredID:
        nobar.append(i)

fig = plt.figure()
ax = fig.add_subplot(111)
norm = plt.Normalize(0, 0.5)
ax.set_xlim(-0.05, 2.05)
ax.set_ylim(0.6, 1.0)
ax.set_xlabel('Redshift z')
ax.set_ylabel('c / a')

for haloID in barredID[::7]:
    mte = MTE[haloID]
    A2 = A2list[haloID]
    ca = []
    ba = []
    for i in range(10):
        if type(mte[i]) is not int:
            ca.append(Flatness(mte[i]))
            # ba.append(BtoA(mte[i]))

    if len(ca) == 10:
        seg = xyline(RedShift, ca)
        lc = LineCollection(seg, cmap='rainbow_r', norm=norm)
        lc.set_array(np.array(A2))
        lc.set_linewidth(1)
        line = ax.add_collection(lc)
fig.colorbar(line, ax=ax)
plt.savefig('f:/Linux/local_result/AxisRatio/ARP_tng.png',dpi=400)
