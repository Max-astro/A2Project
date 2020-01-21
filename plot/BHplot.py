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
    if (simu == 'TNG') or (simu == 'tng'):
        ldir = 'f:/Linux/localRUN/tng_DiskMerTree/%d.json' % subhaloID
    else:
        ldir = 'f:/Linux/localRUN/il1_DiskMerTree/%d.json' % subhaloID
    
    with open(ldir) as f:
        data = json.load(f)
    
    Main = np.array(data['Main'])
    return dict(zip(Main[:, 0], Main[:, 1])), np.array(data['Mergers'])

def HistValAndBin(nums, bins, more=0, mask=0):
    if mask == 1:
        reMask = []

    val = []
    tmp = nums[nums < bins[1]]
    if mask == 1:
        reMask.append(nums < bins[1])
    val.append(len(tmp))

    for i in range(1,len(bins)-1):
        tmp = nums[(nums > bins[i]) & (nums <= bins[i+1])]
        val.append(len(tmp))
        if mask == 1:
            reMask.append((nums > bins[i]) & (nums <= bins[i+1]))

    if more == 1:
        tmp = nums[nums > bins[-1]]
        val.append(len(tmp))
        if mask == 1:
            reMask.append(nums > bins[-1])

    if mask == 0:
        return np.array(val)
    else:
        return np.array(val), np.array(reMask)

def ErrorBarMedian(data):
    #return 25%, 50%, 75%
    if len(data) == 0:
        return 0, 0, 0
    elif len(data) < 3:
        return 0, np.median(data), 0
    else:
        data.sort()
        return data[int(len(data) / 4)], np.median(data), data[int(len(data) * 0.75)]

def Ydata(simu, ids, rawdata, SnapList, Redshift):
    dataWithZ = {}
    for x in Redshift:
        dataWithZ[x] = []
    #all each halo's information into 'dataWithZ'
    for subID in ids:
        data = []
        prog = LoadMergHist(simu, subID)[0]
        plot = 1
        for snap in SnapList:
            try:
                haloID = prog[snap]
                data.append(rawdata[snap][haloID])
            except:
                plot = 0       
            
        if plot:
            for i in range(len(data)):
                dataWithZ[Redshift[i]].append(data[i])

    #calculate Error bar
    plotdata = [[], [], []]
    for i in range(len(data)):
        d0, d1, d2 = ErrorBarMedian(dataWithZ[Redshift[i]])
        plotdata[0].append(d0)
        plotdata[1].append(d1)
        plotdata[2].append(d2)
    plotdata = np.array(plotdata)
    Err = np.vstack((plotdata[1,:] - plotdata[0,:], plotdata[2,:] - plotdata[1,:]))
    return plotdata[1,:], Err
    
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

il1_snapshot = [135, 127, 120, 113, 103, 108, 95, 85, 75, 68]
tng_snapshot = [99, 91, 84, 78, 72, 67, 59, 50, 40, 33]
Redshift = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0]



tng2il1 = np.load('F:/Linux/localRUN/Match/tng2il1_allsub.npy', allow_pickle=1).item()
A2list = np.load('f:/Linux/localRUN/il1_A2withRedshift.npy', allow_pickle=1).item()

bar2bar = np.load('F:/npy/bar2bar.npy',allow_pickle=1).item()
bar2disk = np.load('f:/npy/bar2no.npy',allow_pickle=1).item()

il1_barID = np.load('f:/Linux/localRUN/barredID_il1.npy')
il1_diskID = np.load('f:/Linux/localRUN/diskID_il1.npy')

tng_barID = np.load('f:/Linux/localRUN/barredID_4WP_TNG.npy')
tng_diskID = np.load('f:/Linux/localRUN/diskID_4WP.npy')

tng_BHMass = {}
tng_BHdot = {}
for snap in tng_snapshot:
    mas = il.func.loadSubhalos('TNG', snap, 'SubhaloBHMass')
    dot = il.func.loadSubhalos('TNG', snap, 'SubhaloBHMdot')
    tng_BHMass[snap] = mas
    tng_BHdot[snap] = dot

il1_BHMass = {}
il1_BHdot = {}
for snap in il1_snapshot:
    mas = il.func.loadSubhalos('il1', snap, 'SubhaloBHMass')
    dot = il.func.loadSubhalos('il1', snap, 'SubhaloBHMdot')
    il1_BHMass[snap] = mas
    il1_BHdot[snap] = dot



def BHmassWithZ():
    tng_ids = tng_diskID
    il1_ids = il1_diskID

    tng_mass_Y, tng_Err = Ydata('TNG', tng_ids, tng_BHMass, tng_snapshot, Redshift)
    il1_mass_Y, il1_Err = Ydata('il1', il1_ids, il1_BHMass, il1_snapshot, Redshift)


    #plot info
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Z')
    ax.set_ylabel(r'Subhalo BHmass $10^{10} M_\odot/h$')
    ax.set_yscale("log")
    # ax.set_xlim(-0.4, 0.4)
    # ax.set_ylim(-0.4, 0.4)
    ax.set_title("Disk galaxies BlackHole mass")

    #lines
    ax.errorbar(Redshift, tng_mass_Y, yerr=tng_Err, elinewidth=2, capthick=2, capsize=3, color='r', fmt='o', label='TNG')
    ax.errorbar(Redshift, il1_mass_Y, yerr=il1_Err, elinewidth=2, capthick=2, capsize=3, color='c', fmt='o', label='Illustris-1')
    ax.legend()
    plt.savefig('f:/Linux/local_result/BH/BHmass.png',dpi = 300)



def BHdotWithZ():
    tng_ids = tng_diskID
    il1_ids = il1_diskID

    tng_mass_Y, tng_Err = Ydata('TNG', tng_ids, tng_BHdot, tng_snapshot, Redshift)
    il1_mass_Y, il1_Err = Ydata('il1', il1_ids, il1_BHdot, il1_snapshot, Redshift)


    #plot info
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Z')
    ax.set_ylabel(r'Subhalo BHdot ($(10^{10} M_\odot/h) / (0.978Gyr/h)$)')
    ax.set_yscale("log")
    # ax.set_xlim(-0.4, 0.4)
    # ax.set_ylim(-0.4, 0.4)
    ax.set_title(r"Disk galaxies accretion rates $\dot{M}$")

    #lines
    ax.errorbar(Redshift, tng_mass_Y, yerr=tng_Err, elinewidth=2, capthick=2, capsize=3, color='r', fmt='o', label='TNG')
    ax.errorbar(Redshift, il1_mass_Y, yerr=il1_Err, elinewidth=2, capthick=2, capsize=3, color='c', fmt='o', label='Illustris-1')
    ax.legend()
    plt.savefig('f:/Linux/local_result/BH/BHdot.png',dpi=300)
