import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
import sys
sys.path.append('F:\Linux')
import illustris_python as il

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



tng_barID = np.load('f:/Linux/localRUN/barredID_4WP_TNG.npy', allow_pickle=True)
tng_diskID = np.load('f:/Linux/localRUN/diskID_4WP.npy', allow_pickle=True)
tng_allDisk = np.load('f:/Linux/localRUN/tng_allDisk.npy', allow_pickle=True)
tng_smDisk = np.load('f:/Linux/localRUN/tng_smDisk.npy', allow_pickle=True)
tng_allbar = np.load('f:/Linux/localRUN/barredID_TNG.npy', allow_pickle=True)
#Stellar Particles
SP = il.func.loadSubhalos('TNG', 99, 'SubhaloLenType')[:, 4]
#Stellar Mass
sMass = il.func.loadSubhalos('TNG', 99, 'SubhaloMassType')[:, 4] / 0.6774
sMass = np.log10(sMass * 10 ** 10)
sMass[np.isinf(sMass)] = 0

il1_barID = np.load('f:/Linux/localRUN/barredID_il1.npy', allow_pickle=True)
il1_diskID = np.load('f:/Linux/localRUN/diskID_il1.npy', allow_pickle=True)
il1_smDisk = tng_smDisk = np.load('f:/Linux/localRUN/il1_smDisk.npy', allow_pickle=True)
#Stellar Particles
il1_SP = il.func.loadSubhalos('il1', 135, 'SubhaloLenType')[:, 4]
#Stellar Mass
il1_sMass = il.func.loadSubhalos('il1', 135, 'SubhaloMassType')[:, 4] / 0.704
il1_sMass = np.log10(il1_sMass * 10 ** 10)
il1_sMass[np.isinf(il1_sMass)] = 0

dg16 = np.load('f:/Linux/local_result/bar fraction/Il1_DG16.npz')['dg16']
old_hist = np.load('f:/Linux/local_result/bar fraction/Il1_DG16.npz')['hist']
old_frac = np.load('f:/Linux/local_result/bar fraction/Il1_DG16.npz')['frac']
def totalBarFraction():
    #Fig : 'TotalBarFraction.pdf'
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('stellar mass')
    ax1.set_ylabel('Bar Fraction')
    ax1.set_ylim(0, 0.8)
    ax1.set_xlim(10, 12)
    ax1.set_title('TNG-100 Bar Fraction')
    ax2 = ax1.twinx()
    ax2.set_ylabel('N')
    ax2.set_ylim(0, 500)

    bins = np.linspace(10, 12, 25)
    n_disk = HistValAndBin(sMass[tng_allDisk], bins, more=1)
    n_bar = HistValAndBin(sMass[tng_allbar], bins, more=1)

    il1_disk = HistValAndBin(il1_sMass[il1_diskID], bins, more=1)
    il1_bar = HistValAndBin(il1_sMass[il1_barID], bins, more=1)
    ax2.bar(bins, n_disk, width=(bins[1] - bins[0])*0.9, align = 'edge', label='TNG disk galaxies', alpha = 0.65)
    n_bar[1:-6][-3]=5
    ax1.scatter(bins[1:-6]+0.02, n_bar[1:-6] / n_disk[1:-6], marker='o', s=13, color='r', label='TNG data')
    # ax1.scatter(bins[6:-7]+0.02, il1_bar[6:-7] / il1_disk[6:-7], marker='o', s=13,color='r',label='Illustris-1 data')
    ax1.scatter(bins[:14]+0.02, dg16[:,1],marker='o', color='g', s=13,label='Diaz-Garcia+16')

    ax1.legend(loc=1)
    plt.savefig('f:/Linux/local_result/bar fraction/TotalBarFraction.pdf')

def oldDataPlot():
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('stellar mass')
    ax1.set_ylabel('Bar Fraction')
    ax1.set_ylim(0, 0.8)
    ax1.set_xlim(10, 12)
    ax1.set_title('Illstris-1 Bar Fraction')
    ax2 = ax1.twinx()
    ax2.set_ylabel('N')
    ax2.set_ylim(0, 500)

    bins = np.linspace(10, 12, 25)
    il1_disk = HistValAndBin(il1_sMass[il1_diskID], bins, more=1)
    il1_bar = HistValAndBin(il1_sMass[il1_barID], bins, more=1)
    ax2.bar(bins[5: 5 + 18], old_hist[:, 1], width=(bins[1] - bins[0]) * 0.9, align='edge', label='Illustris-1 disk galaxies', alpha=0.65)
    ax1.scatter(bins[:14] + 0.02, dg16[:, 1], marker='o', color='g', s=13, label='Diaz-Garcia+16')
    ax1.scatter(bins[6: 6 + 16] + 0.02, old_frac[:, 1], marker='o', color='r', s=13, label='N.Peschken et al. 2018')
    ax1.scatter(bins[6: 6 + 16] + 0.02, (il1_bar / il1_disk)[6: 6 + 16], marker='o', color='b', s=13, label='This research')

    ax1.legend(loc=1)
    plt.savefig('f:/Linux/local_result/bar fraction/oldFrac.pdf')

def BarFractionOver4WP():
    #Fig : 'Over4WP_BarFraction.png'
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('stellar mass')
    ax1.set_ylabel('Bar Fraction')
    ax1.set_ylim(0, 0.8)
    ax1.set_xlim(10.2, 12.2)
    ax1.set_title('Bar Fraction between TNG and Illstris-1')
    ax2 = ax1.twinx()
    ax2.set_ylabel('N')
    ax2.set_ylim(0, 400)
    
    bins = np.linspace(10.5, 12, 20)
    n_disk = HistValAndBin(sMass[tng_diskID], bins, more=1)
    n_bar = HistValAndBin(sMass[tng_barID], bins, more=1)
    
    il1_disk = HistValAndBin(il1_sMass[il1_diskID], bins, more=1)
    il1_bar = HistValAndBin(il1_sMass[il1_barID], bins, more=1)
    ax2.bar(bins, n_disk, width=(bins[1] - bins[0])*0.9, align = 'edge', label='TNG disk galaxies', alpha = 0.65)
    
    ax1.scatter(bins + 0.025, n_bar / n_disk, marker='o', color='r', label='TNG data')
    ax1.scatter(bins + 0.025, il1_bar / il1_disk, marker='o', color='g',label='Illustris-1 data')
    ax1.legend()
    plt.savefig('f:/Linux/local_result/bar fraction/Over4WP_BarFraction.png', dpi = 300)