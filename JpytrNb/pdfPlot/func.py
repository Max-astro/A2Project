import numpy as np
import h5py
import json
import sys
sys.path.append("C:/Users/qq651/OneDrive/Codes/")
sys.path.append("C:/Users/qq651/OneDrive/Codes/A2project")
import illustris_python as il
import matplotlib.pyplot as plt
from plotTools.plot import *
import illustrisAPI as iapi


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

def e_const(snap, simuData):
    M=1000 #g
    L=100 #cm
    T=1 #s
    V = 100  #cms-1
    h = simuData['h']
    a = simuData['Redshifts'][snap, 2]
    # t = simuData['Redshifts'][snap, 3] * 1e3 * 1e9 * 31556926 * T / h
    return (1e10*1.989e30*M/h)*(a*3.086e+19*L/h)**2/((0.978*1e9*31556926*T/h)**2)

def selec(left, right, data):
    #return a mask
    return (data >= left) & (data < right)

#get data function
def getBH_cutoff(simu, IDs, barreds, field):
    if (simu == 'TNG') or (simu =='tng'):
        cons = tng_cons
        t = tng_t
        snapList = tng_snap
        path = 'f:/Linux/TNG_cutoff/bhs/'
    else:
        cons = il1_cons
        t = il1_t
        snapList = il1_snap
        path = 'f:/Linux/il1_bh_cutoff/'

    #init result data
    unbar_rawdata = []
    bar_rawdata = []
    for i in snapList:
        unbar_rawdata.append([])
        bar_rawdata.append([])
    #get data by id
    for subID in IDs:
        isdata = 1
        prog = LoadMergHist(simu, subID)[0]
        tmp = []
        last = 0
        t_last = 0
        for snap in snapList[::-1]:
            try:
                haloID = prog[snap]
                f = h5py.File(path + 'snap_%d/cutout_%d.hdf5'%(snap, haloID), 'r')
                engy = np.array(f['PartType5'][field]).sum() * cons[snap] - last
                last = np.array(f['PartType5'][field]).sum() * cons[snap]
                delta_t = t[snap] - t_last
                t_last = t[snap]
            except:
                # print(sys.exc_info()[0])
                isdata = 0
                break
            if delta_t < 1e-9:
                tmp.append(0)
            else:
                tmp.append(engy / delta_t)
        #only use halos which data in every snapshot
        if isdata:
            if subID in barreds:
                for i in range(len(tmp)):
                    bar_rawdata[i].append(tmp[len(tmp)-1-i])
            else:
                for i in range(len(tmp)):
                    unbar_rawdata[i].append(tmp[len(tmp)-1-i])
    return bar_rawdata, unbar_rawdata

def Y_rawdataProcess(data, n_snap):
    plotdata = [[], [], []]
    for i in range(n_snap):
        d0, d1, d2 = ErrorBarMedian(data[i, :])
        plotdata[0].append(d0)
        plotdata[1].append(d1)
        plotdata[2].append(d2)
    plotdata = np.array(plotdata)
    Err = np.vstack((plotdata[1,:] - plotdata[0,:], plotdata[2,:] - plotdata[1,:]))
    return plotdata[1, :], Err

def getGrData(simu, snapList, fields, haloType='sub'):
    raw = {}
    for snap in snapList:
        if haloType == 'fof':
            tmp = il.func.loadhalos(simu, snap, fields)
        else:
            tmp = il.func.loadSubhalos(simu, snap, fields)
        raw[snap] = tmp
    return raw

def a2threshold(haloID, thres, a2list):
    if a2list[haloID][0] >= thres:
        return True
    return False