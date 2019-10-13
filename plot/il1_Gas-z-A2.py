import numpy as np
import h5py
import json
import sys
sys.path.append('F:\Linux')
import illustris_python as il
from datetime import datetime

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import BoundaryNorm, ListedColormap


def xyline(x, y):
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments

# def LoadMergHist(simu, subhaloID):
#     '''
#     return subhalo's main progenitor and merger history with snapshot
#     '''
#     if simu == 'TNG':
#         ldir = '/Raid0/zhouzb/merg_data/tng_DiskMerTree/%d.json' % subhaloID
#     else:
#         ldir = '/Raid0/zhouzb/merg_data/il1_DiskMerTree/%d.json' % subhaloID
    
#     with open(ldir) as f:
#         data = json.load(f)
    
#     Main = np.array(data['Main'])
#     return dict(zip(Main[:, 0], Main[:, 1])), np.array(data['Mergers'])

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


'''
snap_135 z=0
snap_127 z=0.1
snap_120 z=0.2
snap_113 z=0.31
snap_108 z=0.4
snap_103 z=0.5
snap_85 z=1.0
snap_75 z=1.53
snap_68 z=2.0
'''

snapshot = [135, 127, 120, 113, 103, 108, 95, 85, 75, 68]
Redshift = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0]

GasFraction = {}
for snap in snapshot:
    mas = il.func.loadSubhalos('il1', snap, 'SubhaloMassInHalfRadType')
    Gf = mas[:, 0] / (mas[:, 4] + mas[:, 0])
    Gf[np.isnan(Gf)] = 0
    GasFraction[snap] = Gf
    
barID = np.load('f:/Linux/localRUN/barredID_il1.npy')
A2list = np.load('f:/Linux/localRUN/il1_A2withRedshift.npy').item()

fig = plt.figure()
ax = fig.add_subplot(111)
norm = plt.Normalize(0, 0.5)
ax.set_xlim(-0.05, 2.05)
ax.set_ylim(0, 0.9)
ax.set_xlabel('Redshift z')
ax.set_ylabel('Gas Fraction')
ax.set_title('il1 Gas Fraction with Redshift')

for subhaloID in barID:
    # if GasFraction[33][subhaloID] < 0.25:
    #     continue
    prog = LoadMergHist('il1', subhaloID)[0]
    GFlist = []
    A2 = A2list[subhaloID]
    for snapnum in snapshot:
        try:
            haloID = prog[snapnum]
        except:
            GFlist.append(-1)
            continue
        GFlist.append(GasFraction[snapnum][haloID])

    xaxis = []
    yaxis = []
    color = []
    for i in range(10):
        if GFlist[i] != -1:
            xaxis.append(Redshift[i])
            yaxis.append(GFlist[i])
            color.append(float(A2[i]))
    
    seg = xyline(xaxis, yaxis)
    lc = LineCollection(seg, cmap='rainbow_r', norm=norm)
    lc.set_array(np.array(color))
    lc.set_linewidth(1)
    line = ax.add_collection(lc)
fig.colorbar(line, ax=ax)
plt.savefig('f:/Linux/local_result/Gas-Z-A2_il1/GF_all.png',dpi=400)






