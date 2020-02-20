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


def logMasSun(data):
    if type(data) != type(np.array(0)):
        data = np.array(data)
    data = np.log10(data * 10 ** 10)
    data[np.isinf(data)] = 0
    return data

def logmass(data):
    if type(data) != type(np.array(0)):
        data = np.array(data)
    data = np.log10(data)
    data[np.isinf(data)] = 0
    return data

def e_const(snap, simuData):
    h = simuData['h']
    a = simuData['Redshifts'][snap, 2]
    # t = simuData['Redshifts'][snap, 3] * 1e3 * 1e9 * 31556926 * T / h
    return (1e10*1.989e30*M/h)*(a*3.086e+19*L/h)**2/((0.978*1e9*31556926*T/h)**2)

def selec(left, right, data):
    #return a mask
    return (data >= left) & (data < right)

#get data function
def getCF_Data(simu, IDs, barreds):
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
                engy = np.array(f['PartType5']['BH_CumEgyInjection_QM']).sum() * cons[snap] - last
                last = np.array(f['PartType5']['BH_CumEgyInjection_QM']).sum() * cons[snap]
                delta_t =t[snap] - t_last
                t_last = t[snap]
            except:
                # print(sys.exc_info()[0])
                isdata = 0
                break
            tmp.append(engy / delta_t)
        if isdata:
            if subID in barreds:
                for i in range(len(tmp)):
                    bar_rawdata[i].append(tmp[len(tmp)-1-i])
            else:
                for i in range(len(tmp)):
                    unbar_rawdata[i].append(tmp[len(tmp)-1-i])
    return bar_rawdata, unbar_rawdata


M=1000 #g
L=100 #cm
T=1 #s
V = 100  #cms-1

rs = np.array([0, 0.2, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5, 3.0])
il1_snap = [135, 120, 108, 95, 85, 75, 68, 64, 60]
tng_snap = [99, 84, 67, 59, 50, 40, 33, 29, 25]

il1 = np.load('f:/Linux/localRUN/il1SimuData.npy', allow_pickle=True).item()
tng = np.load('f:/Linux/localRUN/tngSimuData.npy', allow_pickle=True).item()
il1_barID = np.load('f:/Linux/localRUN/barredID_il1.npy')
il1_diskID = np.load('f:/Linux/localRUN/diskID_il1.npy')
tng_barID = np.load('f:/Linux/localRUN/barredID_4WP_TNG.npy')
tng_diskID = np.load('f:/Linux/localRUN/diskID_4WP.npy')

tng_unbar = []
for i in tng_diskID:
    if i not in tng_barID:
        tng_unbar.append(i)

il1_unbar = []
for i in il1_diskID:
    if i not in il1_barID:
        il1_unbar.append(i)

tng_cons = {}
tng_t = {}
for snap in tng_snap:
    tng_cons[snap] = e_const(snap, tng)
    tng_t[snap] = tng['Redshifts'][snap, 3] * 1e3 * 1e9 * 31556926 * T / tng['h']

il1_cons = {}
il1_t = {}
for snap in il1_snap:
    il1_cons[snap] = e_const(snap, il1)
    il1_t[snap] = il1['Redshifts'][snap, 3] * 1e3 * 1e9 * 31556926 * T / il1['h']


#Mass bins
tng_sMass = il.func.loadSubhalos('TNG', 99, 'SubhaloMassType')[:, 4]
tng_sMass = logMasSun(tng_sMass)
il1_sMass = il.func.loadSubhalos('il1', 135, 'SubhaloMassType')[:, 4]
il1_sMass = logMasSun(il1_sMass)

b1 = tng_diskID[selec(10.5, 10.75, sMass[tng_diskID])]
b2 = tng_diskID[selec(10.75, 11, sMass[tng_diskID])]
b3 = tng_diskID[selec(11, 14, sMass[tng_diskID])]

bi1 = il1_diskID[selec(10.5, 10.75, il1_sMass[il1_diskID])]
bi2 = il1_diskID[selec(10.75, 11, il1_sMass[il1_diskID])]
bi3 = il1_diskID[selec(11, 14, il1_sMass[il1_diskID])]