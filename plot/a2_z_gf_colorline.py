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

# dolist = [0, 2, 5, 6, 7, 8, 9]
# rs = np.array([0, 0.2, 0.5, 0.7, 1.0, 1.5, 2.0])

Redshift = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0]
il1_snap = [135, 127, 120, 113, 103, 108, 95, 85, 75, 68]
tng_snap = [99, 84, 67, 59, 50, 40, 33]
[99, 91, 84, 78, 72, 67, 63, 59, 50, 40, 33]

#-----------------------------------------------------------------------------
tng_GF = {}
for snap in tng_snap:
    mas = il.func.loadSubhalos('TNG', snap, 'SubhaloMassInHalfRadType')
    Gf = mas[:, 0] / (mas[:, 4] + mas[:, 0])
    Gf[np.isnan(Gf)] = 0
    tng_GF[snap] = Gf
il1_GF = {}
for snap in il1_snap:
    mas = il.func.loadSubhalos('il1', snap, 'SubhaloMassInHalfRadType')
    Gf = mas[:, 0] / (mas[:, 4] + mas[:, 0])
    Gf[np.isnan(Gf)] = 0
    il1_GF[snap] = Gf
#-----------------------------------------------------------------------------
    
il1_A2list = np.load('f:/Linux/localRUN/il1_A2dict(135-68_21part).npy', allow_pickle=1).item()
tng_A2list = np.load('f:/Linux/localRUN/tng_A2dict(99-33_21part).npy', allow_pickle=1).item()

il1_barID = np.load('f:/Linux/localRUN/barredID_il1.npy')
il1_diskID = np.load('f:/Linux/localRUN/diskID_il1.npy')

tng_barID = np.load('f:/Linux/localRUN/barredID_4WP_TNG.npy')
tng_diskID = np.load('f:/Linux/localRUN/diskID_4WP.npy')

#---------------TNG--------------------------------------------------------------


fig = plt.figure()
ax = fig.add_subplot(111)
norm = plt.Normalize(0, 0.5)
ax.set_xlim(-0.05, 2.05)
ax.set_ylim(0, 0.9)
ax.set_xlabel('Redshift z')
ax.set_ylabel('A2')
# ax.set_title('il1 Gas Fraction with Redshift')

for subhaloID in tng_barID:
    prog = LoadMergHist('TNG', subhaloID)[0]
    GFlist = []
    A2 = np.array(tng_A2list[subhaloID])[dolist]
    for snapnum in tng_snap:
        try:
            haloID = prog[snapnum]
        except:
            GFlist.append(-1)
            continue
        GFlist.append(tng_GF[snapnum][haloID])

    xaxis = []
    yaxis = []
    color = []
    for i in range(len(GFlist)):
        if GFlist[i] != -1:
            xaxis.append(rs[i])
            color.append(GFlist[i])
            yaxis.append(float(A2[i]))
    
    seg = xyline(xaxis, yaxis)
    lc = LineCollection(seg, cmap='rainbow_r', norm=norm, linewidth=0.1)
    lc.set_array(np.array(color))
    lc.set_linewidth(0.1)
    line = ax.add_collection(lc)
fig.colorbar(line, ax=ax)
plt.savefig('f:/Linux/local_result/Gas-Z-A2/A2-Z-GAS/tng_barred.pdf')


#---------------Illustris--------------------------------------------------------------

fig = plt.figure()
ax = fig.add_subplot(111)
norm = plt.Normalize(0, 0.5)
ax.set_xlim(-0.05, 2.05)
ax.set_ylim(0, 0.9)
ax.set_xlabel('Redshift z')
ax.set_ylabel('A2')
# ax.set_title('il1 Gas Fraction with Redshift')

for subhaloID in il1_barID:
    prog = LoadMergHist('il1', subhaloID)[0]
    GFlist = []
    A2 = np.array(il1_A2list[subhaloID])[dolist]
    for snapnum in il1_snap:
        try:
            haloID = prog[snapnum]
        except:
            GFlist.append(-1)
            continue
        GFlist.append(il1_GF[snapnum][haloID])

    xaxis = []
    yaxis = []
    color = []
    for i in range(len(GFlist)):
        if GFlist[i] != -1:
            xaxis.append(rs[i])
            color.append(GFlist[i])
            yaxis.append(float(A2[i]))
    
    seg = xyline(xaxis, yaxis)
    lc = LineCollection(seg, cmap='rainbow_r', norm=norm)
    lc.set_array(np.array(color))
    lc.set_linewidth(0.1)
    line = ax.add_collection(lc)
fig.colorbar(line, ax=ax)
plt.savefig('f:/Linux/local_result/Gas-Z-A2/A2-Z-GAS/il1_barred.pdf')


