import numpy as np
import h5py
import json
import sys
sys.path.append('F:/Linux')
sys.path.append("C:/Users/qq651/OneDrive/Codes/A2project/plot/")
import illustris_python as il
import matplotlib.pyplot as plt
from Tools import *

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

il1_A2list = np.load('f:/Linux/localRUN/il1_A2dict(135-68_21part).npy', allow_pickle=1).item()
tng_A2list = np.load('f:/Linux/localRUN/tng_A2dict(99-33_21part).npy', allow_pickle=1).item()

il1_snap = [135, 127, 120, 113, 103, 108, 99, 95, 92, 89, 85, 82, 80, 78, 76, 73, 71, 70, 69, 68]
tng_snap = [99, 91, 84, 78, 72, 67, 63, 59, 56, 53, 50, 47, 45, 43, 41, 40, 38, 36, 35, 34, 33]
Redshift = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]


il1_tbar = []
tng_tbar = []

il1_id = []
tng_id = []

for i in il1_A2list.keys():
    il1_tbar.append(zbar(i, il1_A2list))
    prog = LoadMergHist('il1', i)[0]
    try:
        il1_id.append(prog[68])
    except KeyError:
        continue


for i in tng_A2list.keys():
    prog = LoadMergHist('TNG', i)[0]
    try:
        tng_id.append(prog[33])
        tng_tbar.append(zbar(i, tng_A2list))
    except KeyError:
        continue


mas = il.func.loadSubhalos('TNG', 33, 'SubhaloMassInHalfRadType')
Gf = mas[:, 0] / (mas[:, 4] + mas[:, 0])
Gf[np.isnan(Gf)] = 0
tng_gf = Gf


mas = il.func.loadSubhalos('il1', 68, 'SubhaloMassInHalfRadType')
Gf = mas[:, 0] / (mas[:, 4] + mas[:, 0])
Gf[np.isnan(Gf)] = 0
il1_gf = Gf

il1_tbar = np.array(il1_tbar)
tng_tbar = np.array(tng_tbar)

def zbarplot():
    plt.scatter(tng_gf[tng_id], tng_tbar, s=4, label='TNG-100')
    plt.scatter(il1_gf[il1_id], il1_tbar, s=6, color='r', label='Illustris-1')



    plt.xlabel('Gas fraction at z=2')
    plt.ylabel('zbar')
    plt.xlim(-0.02, 1)
    plt.legend()
    plt.savefig('F:/Linux/local_result/zbar.eps')


def zbarErrPlot():
    tng_dots = np.vstack((tng_gf[tng_id], tng_tbar))
    il1_dots = np.vstack((il1_gf[il1_id], il1_tbar))

    bins = np.linspace(0, 1, 10)
    tng_plotdata = [[], [], []]
    for i in range(9):
        mask = (tng_dots[0, :] >= bins[i]) & (tng_dots[0, :] < bins[i + 1])
        tmp = tng_dots[1, :][mask]
        d0, d1, d2 = ErrorBarMedian(tmp)
        tng_plotdata[0].append(d0)
        tng_plotdata[1].append(d1)
        tng_plotdata[2].append(d2)
    tng_plotdata = np.array(tng_plotdata)
    tng_Err = np.vstack((tng_plotdata[1, :] - tng_plotdata[0, :], tng_plotdata[2, :] - tng_plotdata[1, :]))

    il1_plotdata = [[], [], []]
    for i in range(9):
        mask = (il1_dots[0, :] >= bins[i]) & (il1_dots[0, :] < bins[i + 1])
        tmp = il1_dots[1, :][mask]
        d0, d1, d2 = ErrorBarMedian(tmp)
        il1_plotdata[0].append(d0)
        il1_plotdata[1].append(d1)
        il1_plotdata[2].append(d2)
    il1_plotdata = np.array(il1_plotdata)
    il1_Err = np.vstack((il1_plotdata[1, :] - il1_plotdata[0, :], il1_plotdata[2, :] - il1_plotdata[1, :]))

    plt.errorbar(bins[:-1], tng_plotdata[1,:], yerr = tng_Err, elinewidth=2, capthick=2, capsize=3, color='c', fmt='o', label='TNG-100')
    plt.errorbar(bins[:-1], il1_plotdata[1, :], yerr=il1_Err, elinewidth=2, capthick=2, capsize=3, color='r', fmt='o', label='Illustris-1')

    plt.xlim(-0.02, 1)
    plt.title('Statistics of zbar in each bin')
    plt.xlabel('Gas fraction at z=2')
    plt.ylabel('zbar')
    plt.legend()
    plt.savefig('F:/Linux/local_result/zbar_err.eps')