import numpy as np
import h5py
import json
import sys
sys.path.append('F:/Linux')
import illustris_python as il


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
    return plotdata[1, :], Err


def zbar(haloID, A2list):
    #return bar origin redshift
    Redshift = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
    A2 = np.array(A2list[haloID])
    z=0
    for i in range(1,len(Redshift)): 
        if A2[i] < 0.15:
            break
        z += 1
    while z != 0:
        if abs((A2[z] - A2[z - 1]) / A2[z]) <= 0.4:
            break
        z -= 1

    return Redshift[z]