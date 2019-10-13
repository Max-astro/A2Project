import numpy as np
import h5py
import json
import sys
sys.path.append('F:\Linux')
import illustris_python as il
from datetime import datetime
sys.path.append(r"C:/Users/qq651/OneDrive/Codes/A2project/plot/")
from Tools import *

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import BoundaryNorm, ListedColormap

il1_A2list = np.load('f:/Linux/localRUN/il1_A2withRedshift.npy', allow_pickle=1).item()
tng_A2list = np.load('f:/Linux/localRUN/tng_A2withRedshift.npy', allow_pickle=1).item()
def zbar(haloID, A2list):
    #return bar origin redshift
    Redshift = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0]
    A2 = np.array(A2list[haloID])
    z=0
    for i in range(1,10): 
        if A2[i] < 0.15:
            break
        z += 1
    while z != 0:
        if abs((A2[z] - A2[z - 1]) / A2[z]) <= 0.4:
            break
        z -= 1
    if z != 9:
        z += 1
    return Redshift[z]

il1_tbar = []
tng_tbar = []
for i in il1_A2list.keys():
    il1_tbar.append(zbar(i, il1_A2list))
for i in tng_A2list.keys():
    tng_tbar.append(zbar(i, tng_A2list))


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

il1_snapshot = [135, 127, 120, 113, 103, 108, 95, 85, 75, 68, 64, 60]
tng_snapshot = [99, 91, 84, 78, 72, 67, 59, 50, 40, 33, 29, 25]
Redshift = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0])

tng2il1 = np.load('F:/Linux/localRUN/Match/tng2il1_allsub.npy',allow_pickle=1).item()
A2list = np.load('f:/Linux/localRUN/il1_A2withRedshift.npy',allow_pickle=1).item()

bar2bar = np.load('F:/npy/bar2bar.npy',allow_pickle=1).item()
bar2disk = np.load('f:/npy/bar2no.npy',allow_pickle=1).item()

il1_barID = np.load('f:/Linux/localRUN/barredID_il1.npy')
il1_diskID = np.load('f:/Linux/localRUN/diskID_il1.npy')

tng_barID = np.load('f:/Linux/localRUN/barredID_4WP_TNG.npy')
tng_diskID = np.load('f:/Linux/localRUN/diskID_4WP.npy')

tng_SFR = {}
tng_GF = {}
for snap in tng_snapshot:
    sfr = il.func.loadSubhalos('TNG', snap, 'SubhaloSFR')
    mas = il.func.loadSubhalos('TNG', snap, 'SubhaloMassInHalfRadType')
    Gf = mas[:, 0] / (mas[:, 4] + mas[:, 0])
    Gf[np.isnan(Gf)] = 0
    tng_GF[snap] = Gf
    tng_SFR[snap] = sfr

il1_SFR = {}
il1_GF = {}
for snap in il1_snapshot:
    sfr = il.func.loadSubhalos('il1', snap, 'SubhaloSFR')
    mas = il.func.loadSubhalos('il1', snap, 'SubhaloMassInHalfRadType')
    Gf = mas[:, 0] / (mas[:, 4] + mas[:, 0])
    Gf[np.isnan(Gf)] = 0
    il1_GF[snap] = Gf
    il1_SFR[snap] = sfr


def SFR_GAS_WithZ():
    tng_ids = tng_diskID
    il1_ids = il1_diskID

    tng_SFR_Y, tng_SFR_Err = Ydata('TNG', tng_ids, tng_SFR, tng_snapshot, Redshift)
    il1_SFR_Y, il1_SFR_Err = Ydata('il1', il1_ids, il1_SFR, il1_snapshot, Redshift)

    # tng_GF_Y, tng_GF_Err = Ydata('TNG', tng_ids, tng_GF, tng_snapshot, Redshift)
    # il1_GF_Y, il1_GF_Err = Ydata('il1', il1_ids, il1_GF, il1_snapshot, Redshift)

    #plot info
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Z')
    ax.set_ylabel(r'Star formation rate $M_\odot/yr$')
    # ax.set_yscale("log")
    # ax.set_xlim(-0.4, 0.4)
    # ax.set_ylim(-0.4, 0.4)
    ax.set_title(r"Star formation rate evolve with redshift")

    #lines
    p1 = ax.errorbar(Redshift, tng_SFR_Y, yerr=tng_SFR_Err, elinewidth=2, capthick=2, capsize=3, color='r', fmt='o', label='TNG SFR')
    # ax2 = ax.twinx()
    ax2.set_ylabel(r'Gas Fraction $10^{10} M_\odot$')
    
    # p2 = ax2.errorbar(Redshift, tng_GF_Y, yerr=tng_GF_Err, elinewidth=2, capthick=2, capsize=3, color='c', fmt='o', label='TNG Gas Fraction')
    # pic = p1 + p2
    
    # plt.savefig('f:/Linux/local_result/tng_SFR.png',dpi=300)

    ax.errorbar(Redshift, il1_SFR_Y, yerr=il1_SFR_Err, elinewidth=2, capthick=2, capsize=3, color='c', fmt='o', label='Illustris-1 SFR')
    ax2.errorbar(Redshift, tng_GF_Y, yerr=il1_GF_Err, elinewidth=2, capthick=2, capsize=3, color='c', fmt='o', label='Illustris-1 Gas Fraction')
    ax.legend()
    # plt.savefig('f:/Linux/local_result/BH/BHdot.png',dpi=300)


def il1_tng_SFR():
    tng_ids = tng_barID
    il1_ids = il1_barID

    tng_ubID = []
    for i in tng_diskID:
        if i not in tng_barID:
            tng_ubID.append(i)
    il1_ubID = []
    for i in il1_diskID:
        if i not in il1_barID:
            il1_ubID.append(i)
            

    tng_SFR_Y, tng_SFR_Err = Ydata('TNG', tng_ids, tng_SFR, tng_snapshot, Redshift)
    il1_SFR_Y, il1_SFR_Err = Ydata('il1', il1_ids, il1_SFR, il1_snapshot, Redshift)

    tng_ub_SFR_GF, tng_ub_SFR_err = Ydata('TNG', tng_ubID, tng_SFR, tng_snapshot, Redshift)
    il1_ub_SFR_GF, il1_ub_SFR_err = Ydata('il1', il1_ubID, il1_SFR, il1_snapshot, Redshift)

    #plot info
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Z')
    ax.set_ylabel(r'Star formation rate $M_\odot/yr$')
    # ax.set_yscale("log")
    # ax.set_xlim(-0.4, 0.4)
    # ax.set_ylim(-0.4, 0.4)
    ax.set_title(r"Star formation rate evolve with redshift")

    #lines

    ax.errorbar(Redshift-0.02, tng_SFR_Y, yerr=tng_SFR_Err, elinewidth=2, capthick=2, capsize=3, color='r', fmt='o', label='TNG barred')
    ax.errorbar(Redshift-0.02, tng_ub_SFR_GF, yerr=tng_ub_SFR_err, elinewidth=2, capthick=2, capsize=3, color='violet', fmt='^', label='TNG no bar')

    ax.errorbar(Redshift+0.02, il1_SFR_Y, yerr=il1_SFR_Err, elinewidth=2, capthick=2, capsize=3, color='c', fmt='o', label='Illustris-1 barred')
    ax.errorbar(Redshift+0.02, il1_ub_SFR_GF, yerr=il1_ub_SFR_err, elinewidth=2, capthick=2, capsize=3, color='blue', fmt='^', label='Illustris-1 nobar')

    ax.legend()
    plt.savefig('f:/Linux/local_result/il1-TNG-SFR.png', dpi=300)

def il1_tng_GF():
    # fig : 'f:/Linux/local_result/il1-TNG-GF.png'
    tng_ids = tng_barID
    il1_ids = il1_barID

    tng_ubID = []
    for i in tng_diskID:
        if i not in tng_barID:
            tng_ubID.append(i)
    il1_ubID = []
    for i in il1_diskID:
        if i not in il1_barID:
            il1_ubID.append(i)
            

    tng_GF_Y, tng_GF_Err = Ydata('TNG', tng_ids, tng_GF, tng_snapshot, Redshift)
    il1_GF_Y, il1_GF_Err = Ydata('il1', il1_ids, il1_GF, il1_snapshot, Redshift)

    tng_ub_GF, tng_ub_GF_err = Ydata('TNG', tng_ubID, tng_GF, tng_snapshot, Redshift)
    il1_ub_GF, il1_ub_GF_err = Ydata('il1', il1_ubID, il1_GF, il1_snapshot, Redshift)

    #plot info
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Z')
    ax.set_ylabel(r'Gas Fraction $10^{10} M_\odot$')
    # ax.set_yscale("log")
    # ax.set_xlim(-0.4, 0.4)
    # ax.set_ylim(-0.4, 0.4)
    ax.set_title(r"Gas fraction evolve with redshift")

    #lines

    ax.errorbar(Redshift-0.02, tng_GF_Y, yerr=tng_GF_Err, elinewidth=2, capthick=2, capsize=3, color='r', fmt='o', label='TNG barred')
    ax.errorbar(Redshift-0.02, tng_ub_GF, yerr=tng_ub_GF_err, elinewidth=2, capthick=2, capsize=3, color='violet', fmt='^', label='TNG no bar')

    ax.errorbar(Redshift+0.02, il1_GF_Y, yerr=il1_GF_Err, elinewidth=2, capthick=2, capsize=3, color='c', fmt='o', label='Illustris-1 barred')
    ax.errorbar(Redshift+0.02, il1_ub_GF, yerr=il1_ub_GF_err, elinewidth=2, capthick=2, capsize=3, color='blue', fmt='^', label='Illustris-1 nobar')

    ax.legend()
    plt.savefig('f:/Linux/local_result/il1-TNG-GF.png', dpi=300)

def trim_axs(axs, N):
    """little helper to massage the axs list to have correct length..."""
    axs = axs.flat
    for ax in axs[N:]:
        ax.remove()
    return axs[:N]

il1_unbar = []
for i in il1_diskID:
    if i not in il1_barID:
        il1_unbar.append(i)

tng_unbar = []
for i in tng_diskID:
    if i not in tng_barID:
        tng_unbar.append(i)

def Global_il1_TNG_GF():
    #Global_il1-TNG_GAS
    il1_Y, il1_Err = Ydata('il1', il1_barID, il1_SFR, il1_snapshot, Redshift)
    il1_Y_2, il1_Err_2 = Ydata('il1', il1_unbar, il1_SFR, il1_snapshot, Redshift)
            
    tng_Y, tng_Err = Ydata('TNG', tng_barID, tng_SFR, tng_snapshot, Redshift)
    tng_Y_2, tng_Err_2 = Ydata('TNG', tng_unbar, tng_SFR, tng_snapshot, Redshift)


    #plot info
    fig = plt.figure(figsize=(8,5))
    ax = fig.add_subplot(111)
    ax.set_xlabel('Z')
    ax.set_ylabel(r'Star formation rate $M_\odot/yr$')
    # ax.set_yscale("log")
    # ax.set_xlim(-0.4, 0.4)
    # ax.set_ylim(0.1, 0.00000001)
    ax.set_title(r"il1&TNG galaxies star formation rate")


    #lines
    ax.errorbar(Redshift-0.015, tng_Y, yerr=tng_Err, elinewidth=2, capthick=2, capsize=3, color='blue', fmt='o', ms=5, label='TNG barred')
    ax.errorbar(Redshift-0.015, tng_Y_2, yerr=tng_Err_2, elinewidth=2, capthick=2, capsize=3, color='grey', fmt='^', ms=5,label='TNG unbar')
    ax.errorbar(Redshift+0.015, il1_Y, yerr=il1_Err, elinewidth=2, capthick=2, capsize=3, color='blue', fmt='o',ms=5, label='il1 barred')
    ax.errorbar(Redshift+0.015, il1_Y_2, yerr=il1_Err_2, elinewidth=2, capthick=2, capsize=3, color='grey', fmt='^',ms=5, label='il1 unbar')
    ax.legend()