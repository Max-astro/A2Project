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
        tmp = nums[(nums >= bins[i]) & (nums < bins[i+1])]
        val.append(len(tmp))
        if mask == 1:
            reMask.append((nums >= bins[i]) & (nums < bins[i+1]))

    if more == 1:
        tmp = nums[nums >= bins[-2]]
        val.pop()
        val.append(len(tmp))
        if mask == 1:
            reMask.append(nums > bins[-1])

    if mask == 0:
        return np.array(val)
    else:
        return np.array(val), np.array(reMask)

def countInBins(nums, bins):
    nn = []
    for i in range(len(bins)-1):
        nn.append(nums[(nums > bins[i]) & (nums < bins[i+1])])
    nn.append(nums[nums > bins[-1]])
    count = []
    for dd in nn:
        count.append(len(dd))
    return np.array(count)



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


def Y_rawdata(data, snapnum):
    plotdata = [[], [], []]
    for i in range(snapnum):
        d0, d1, d2 = ErrorBarMedian(data[i, :])
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

def tForming(simu, haloID):
    il1_snap = [135, 127, 120, 113, 108, 103, 99, 95, 92, 89, 85, 82, 80, 78, 76, 75, 73, 71, 70, 69, 68, 64, 60]
    tng_snap = [99, 91, 84, 78, 72, 67, 63, 59, 56, 53, 50, 47, 45, 43, 41, 40, 38, 36, 35, 34, 33, 29, 25]
    if simu == 'TNG' or simu == 'tng':
        snapList = np.array(tng_snap)
    else:
        snapList = np.array(il1_snap)
    Redshift = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0])
    m0 = il.func.loadSubhalos(simu, snapList[0], 'SubhaloMassType')[haloID, 4]

    prog = LoadMergHist(simu, haloID)[0]
    for snap in snapList:
        try:
            subID = prog[snap]
        except KeyError:
            continue
        if il.func.loadSubhalos(simu, snap, 'SubhaloMassType')[subID, 4] < m0 * 0.5:
            return Redshift[snapList == snap]
    return 3.0

def getMultiData(simu, snapList, fields, haloType='sub'):
    res = {}
    for i in fields:
        res[i] = {}
    for snap in snapList:
        if haloType == 'fof':
            tmp = il.func.loadhalos(simu, snap, fields)
        else:
            tmp = il.func.loadSubhalos(simu, snap, fields)
        for i in fields:
            res[i][snap] = tmp[i]
    return res

def getData(simu, snapList, fields, haloType='sub'):
    raw = {}
    for snap in snapList:
        if haloType == 'fof':
            tmp = il.func.loadhalos(simu, snap, fields)
        else:
            tmp = il.func.loadSubhalos(simu, snap, fields)
        raw[snap] = tmp[i]
    return raw

def log_Msun(data):
    if type(data) != type(np.array(0)):
        data = np.array(data)
    data = np.log10(data * 10 ** 10)
    data[np.isinf(data)] = 0
    data[np.isnan(data)] = 0
    return data

def log10_safe(data):
    if type(data) != type(np.array(0)):
        data = np.array(data)
    data = np.log10(data)
    data[np.isinf(data)] = 0
    data[np.isnan(data)] = 0
    return data