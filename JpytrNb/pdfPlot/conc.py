import numpy as np
import h5py
import json
import sys
sys.path.append('F:/Linux')
sys.path.append("C:/Users/qq651/OneDrive/Codes/A2project/")
import illustris_python as il
import matplotlib.pyplot as plt
# from plotTools.plot import *


def LoadConc(simu, snap):
    if simu == 'TNG' or simu == 'tng':
        path = 'f:/Linux/localRUN/concRes/tng/tng_%d.conc' % snap
    else:
        path = 'f:/Linux/localRUN/concRes/il1/il1_%d.conc' % snap

    return np.fromfile(path, dtype='d')

def Y_rawdata(data, snapnum):
    plotdata = [[], [], []]
    for snap in snapnum:
        d0, d1, d2 = ErrorBarMedian(data[snap])
        plotdata[0].append(d0)
        plotdata[1].append(d1)
        plotdata[2].append(d2)
    plotdata = np.array(plotdata)
    Err = np.vstack((plotdata[1,:] - plotdata[0,:], plotdata[2,:] - plotdata[1,:]))
    return plotdata[1, :], Err

def findGid(simu, snap, ids):
    res = []
    GID = il.func.loadSubhalos(simu, snap, 'SubhaloGrNr')
    for haloID in ids:
        try:
            prog = LoadMergHist(simu, haloID)[0]
            subID = prog[snap]
        except:
            continue

        res.append(GID[subID])
    return np.unique(res)

def ConcData(simu, snap, diskID, ids, data):
    diskGID = findGid(simu, snap, diskID)
    gid = findGid(simu, snap, ids)
    
    ind = []
    for haloID in gid:
        ind.append(np.where(diskGID == haloID)[0][0])
    
    return data[ind]




rs = np.array([0, 0.2, 0.5, 0.7, 1.0, 1.5, 2.0])
il1_snap = [135, 120, 108, 95, 85, 75, 68]
tng_snap = [99, 84, 67, 59, 50, 40, 33]

il1_barID = np.load('f:/Linux/localRUN/barredID_il1.npy', allow_pickle=1)
il1_diskID = np.load('f:/Linux/localRUN/diskID_il1.npy', allow_pickle=1)

tng_barID = np.load('f:/Linux/localRUN/barredID_4WP_TNG.npy', allow_pickle=1)
tng_diskID = np.load('f:/Linux/localRUN/diskID_4WP.npy', allow_pickle=1)
tng_unbarID = []
for haloID in tng_diskID:
    if haloID not in tng_barID:
        tng_unbarID.append(haloID)

il1_unbarID = []
for haloID in il1_diskID:
    if haloID not in il1_barID:
        il1_unbarID.append(haloID)


# GID = il.func.loadSubhalos(simu, snap, 'SubhaloGrNr')

tng_c = {}
for snap in tng_snap:
    tng_c[snap] = LoadConc('tng', snap)

il1_c = {}
for snap in il1_snap:
    il1_c[snap] = LoadConc('il1', snap)



tng_bar_c = {}
tng_unbar_c = {}

for snap in tng_snap:
    tng_bar_c[snap]= (ConcData('TNG', snap, tng_diskID, tng_barID, tng_c[snap]))
    tng_unbar_c[snap] = (ConcData('TNG', snap, tng_diskID, tng_unbarID, tng_c[snap]))
    

il1_bar_c = {}
il1_unbar_c = {}

for snap in il1_snap:
    il1_bar_c[snap]= (ConcData('il1', snap, il1_diskID, il1_barID, il1_c[snap]))
    il1_unbar_c[snap] = (ConcData('il1', snap, il1_diskID, il1_unbarID, il1_c[snap]))
    

def separatePlot():
    plt.figure(figsize=(8,8))
    bar_Y, bar_err = Y_rawdata(tng_bar_c, tng_snap)
    un_Y, un_err = Y_rawdata(tng_unbar_c, tng_snap)

    il1_bar_Y, il1_bar_err = Y_rawdata(il1_bar_c, il1_snap)
    il1_un_Y, il1_un_err = Y_rawdata(il1_unbar_c, il1_snap)


    plt.errorbar(rs, bar_Y, bar_err, elinewidth=2, label='TNG barred', capthick=2, capsize=3, color='r', fmt='o',ms=5, ls='-')
    plt.errorbar(rs, un_Y, un_err, elinewidth=2, label='TNG unbarred', capthick=2, capsize=3, color='orange', fmt='o', ms=5, ls='-')

    plt.errorbar(rs, il1_bar_Y, il1_bar_err, elinewidth=2, label='TNG barred', capthick=2, capsize=3, color='blue', fmt='o',ms=5, ls='-')
    plt.errorbar(rs, il1_un_Y, il1_un_err, elinewidth=2, label='TNG unbarred', capthick=2, capsize=3, color='c', fmt='o', ms=5, ls='-')

    plt.xlim(-0.05, 2.1)
    plt.xlabel('Z', fontsize=22)
    plt.ylabel('c', fontsize=22)
    plt.legend(fontsize=15)
    plt.tick_params(labelsize=16)
    plt.savefig('f:/Linux/local_result/conc/conc_sepa.pdf')


def totalConc():
    tng_c, tng_err = Y_rawdata(tng_c, tng_snap)
    il1_c, il1_err = Y_rawdata(il1_c, il1_snap)
    plt.errorbar(rs, tng_c, tng_err, label='TNG disk galaxies concentration c')
    plt.errorbar(rs, il1_c, il1_err, label='Illustris-1 disk galaxies concentration c')











