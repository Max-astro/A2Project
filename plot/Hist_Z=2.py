import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
import sys
import json
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

def ErrorBarMedian(data):
    #return 25%, 50%, 75%
    if len(data) == 0:
        return 0, 0, 0
    elif len(data) < 3:
        return np.median(data), np.median(data), np.median(data)
    else:
        data.sort()
        return data[int(len(data) / 4)], np.median(data), data[int(len(data) * 0.75)]


#TNG data
snap99_barID = np.load('f:/Linux/localRUN/barredID_4WP_TNG.npy')
snap99_diskID = np.load('f:/Linux/localRUN/diskID_4WP.npy')
#Gas Fraction Gf
mas = il.func.loadSubhalos('TNG', 33, 'SubhaloMassInRadType')
Gf = mas[:, 0] / (mas[:, 4] + mas[:, 0])
Gf[np.isnan(Gf)] = 0
#Stellar Particles
SP = il.func.loadSubhalos('TNG', 33, 'SubhaloLenType')[:, 4]
#Stellar Mass
sMass = il.func.loadSubhalos('TNG', 33, 'SubhaloMassType')[:, 4] / 0.6774
sMass = np.log10(sMass * 10 ** 10)
sMass[np.isinf(sMass)] = 0

#Illsutris-1 data
snap135_barID = np.load('f:/Linux/localRUN/barredID_il1.npy')
snap135_diskID = np.load('f:/Linux/localRUN/diskID_il1.npy')
il1_mas = il.func.loadSubhalos('il1', 68, 'SubhaloMassInRadType')
#Gas Fraction
il1_gf = il1_mas[:, 0] / (il1_mas[:, 4] + il1_mas[:, 0])
il1_gf[np.isnan(il1_gf)] = 0
#Stellar Particles
il1_SP = il.func.loadSubhalos('il1', 68, 'SubhaloLenType')[:, 4]
#Stellar Mass
il1_sMass = il.func.loadSubhalos('il1', 68, 'SubhaloMassType')[:, 4] / 0.704
il1_sMass = np.log10(il1_sMass * 10 ** 10)
il1_sMass[np.isinf(il1_sMass)] = 0

#ID in z=2
barID = []
diskID = []
for haloID in snap99_diskID:
    prog = LoadMergHist('TNG', haloID)[0]
    try:
        progID = prog[33]
    except:
        continue
    if progID != -1:
        diskID.append(progID)
        if haloID in snap99_barID:
            barID.append(progID)

il1_barID = []
il1_diskID = []
for haloID in snap135_diskID:
    prog = LoadMergHist('il1', haloID)[0]
    try:
        progID = prog[68]
    except:
        continue
    if progID != -1:
        il1_diskID.append(progID)
        if haloID in snap135_barID:
            il1_barID.append(progID)

def Plot_TNG_sMassAndGasFraction():
    #
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Stellar Mass')
    ax1.set_ylabel('N')
    ax1.set_title('TNG Gas Fraction with Stellar Mass bin in Z=2')
    ax1.set_xlim(8,11.5)
    ax1.set_ylim(0, 400)
    ax2 = ax1.twinx()
    ax2.set_ylabel('Gas Fraction')
    ax2.set_ylim(0, 1.1)

    bins = np.linspace(8, 11.3, 20)
    n1, m1 = HistValAndBin(sMass[diskID], bins, more=1, mask=1)
    n2, m2 = HistValAndBin(sMass[barID], bins, more=1, mask=1)
    ax1.bar(bins, n1, width=(bins[1] - bins[0]) * 0.9, align='edge', label='disk nobar')
    ax1.bar(bins, n2, width=(bins[1] - bins[0]) * 0.9, align='edge', label='barred galaxies')
    ax1.legend()

    #GasFraction ErrorBar
    Dgas = Gf[diskID]
    data = [[], [], []]
    mean = []
    for i in m1:
        mean.append(np.mean(Dgas[i]))
        d0, d1, d2 = ErrorBarMedian(Dgas[i])
        data[0].append(d0)
        data[1].append(d1)
        data[2].append(d2)
    data = np.array(data)
    yerr = np.vstack((data[1,:] - data[0,:], data[2,:] - data[1,:]))
    x = bins + (bins[1] - bins[0]) / 2
    #ax2.plot(x, mean, color='c', marker='o', label='Mean Gas Fraction')
    ax2.errorbar(x, data[1,:], yerr=yerr, elinewidth=2, capthick=2, capsize=3, color='r', fmt='o', label='Median Gas Fraction Errorbar')
    ax2.legend(loc=2)
    plt.savefig('f:/Linux/local_result/Z=2/TNG_sMassAndGas_Z=2.png', dpi = 300)


def Plot_il1_sMassAndGasFraction():
    #
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Stellar Mass')
    ax1.set_ylabel('N')
    ax1.set_title('Illustris Gas Fraction with Stellar Mass bin')
    ax1.set_xlim(8,11.5)
    ax1.set_ylim(0, 230)
    ax2 = ax1.twinx()
    ax2.set_ylabel('Gas Fraction')
    ax2.set_ylim(0, 1.1)

    bins = np.linspace(8, 11.3, 20)
    n1, m1 = HistValAndBin(il1_sMass[il1_diskID], bins, more=1, mask=1)
    n2, m2 = HistValAndBin(il1_sMass[il1_barID], bins, more=1, mask=1)
    ax1.bar(bins, n1, width=(bins[1] - bins[0]) * 0.9, align='edge', label='disk nobar')
    ax1.bar(bins, n2, width=(bins[1] - bins[0]) * 0.9, align='edge', label='barred galaxies')
    ax1.legend()

    #GasFraction ErrorBar
    Dgas = il1_gf[il1_diskID]
    data = [[], [], []]
    mean = []
    for i in m1:
        mean.append(np.mean(Dgas[i]))
        d0, d1, d2 = ErrorBarMedian(Dgas[i])
        data[0].append(d0)
        data[1].append(d1)
        data[2].append(d2)
    data = np.array(data)
    yerr = np.vstack((data[1,:] - data[0,:], data[2,:] - data[1,:]))
    x = bins + (bins[1] - bins[0]) / 2
    ax2.errorbar(x, data[1,:], yerr=yerr, elinewidth=2, capthick=2, capsize=3, color='r', fmt='o', label='Median Gas Fraction Errorbar')
    ax2.plot(x, mean, color='c', marker='o', label='Mean Gas Fraction')
    ax2.legend(loc=2)
    plt.savefig('f:/Linux/local_result/Z=2/Illustris_sMassAndGas_Mean.png', dpi = 300)


def plot_GasF_TNGAndil1():
    #Fig : 'TNG_Illustris-1_GF_Z=0.png'
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('GasFraction')
    ax1.set_ylabel('N')
    ax1.set_title('TNG & Illustris-1 Gas Fraction & Bar Fraction at Z=2')
    ax1.set_ylim(0, 350)
    ax2 = ax1.twinx()
    ax2.set_ylabel('Total Bar Fraction')
    ax2.set_ylim(0, 0.9)

    bins = np.linspace(0, 1, 15)
    disk_tng_n = HistValAndBin(Gf[diskID], bins)
    bar_tng_n = HistValAndBin(Gf[barID], bins)
    disk_il1_n = HistValAndBin(il1_gf[il1_diskID], bins)
    bar_il1_n = HistValAndBin(il1_gf[il1_barID], bins)
    ax1.bar(bins[:-1], disk_tng_n, width=(bins[1] - bins[0]) * 0.35,align = 'edge', label='TNG disk galaxies')
    ax1.bar(bins[:-1], bar_tng_n, width=(bins[1] - bins[0]) * 0.35, align='edge', label='TNG barred galaxies', color='c')
    ax1.bar(bins[:-1] + 0.02, disk_il1_n, width=(bins[1] - bins[0]) * 0.35, align = 'edge', label='Illustris-1 disk galaxies')
    ax1.bar(bins[:-1] + 0.02, bar_il1_n, width=(bins[1] - bins[0]) * 0.35, align='edge', label='Illustris-1 barred galaxies', color='r')

    frac = bar_tng_n / disk_tng_n
    frac[-3:] = 0
    ax2.plot(bins[:-1] + 0.021, frac, marker='o', label='TNG bar fraction', color='b')
    ax2.plot(bins[:-1] + 0.021, bar_il1_n / disk_il1_n, marker='o', label='Illustris-1 bar fraction', color='k')

    ax1.legend()
    ax2.legend(loc=2)
    plt.savefig('f:/Linux/local_result/TNG_Illustris-1_GF_Z=2.png', dpi=300)

def plot_TNGAndil1_GasFraction():
    #Fig : 'TNG_il1_sMassAndGas_Z=2.png'
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Stellar Mass')
    ax1.set_ylabel('N')
    ax1.set_title('TNG Gas Fraction with Stellar Mass bin in Z=2')
    ax1.set_xlim(8,11.5)
    ax1.set_ylim(0, 250)
    ax2 = ax1.twinx()
    ax2.set_ylabel('Gas Fraction')
    ax2.set_ylim(0, 1.2)

    bins = np.linspace(8, 11.3, 20)
    n1, m1 = HistValAndBin(sMass[diskID], bins, more=1, mask=1)
    n2, m2 = HistValAndBin(il1_sMass[il1_diskID], bins, more=1, mask=1)
    ax1.bar(bins, n1, width=(bins[1] - bins[0]) * 0.35, align='edge', label='TNG disk in Z=0')
    ax1.bar(bins+0.064, n2, width=(bins[1] - bins[0]) * 0.35, align='edge', label='ill-1 disk in Z=0')

    #TNG GasFraction ErrorBar
    Dgas = Gf[diskID]
    data = [[], [], []]
    mean = []
    for i in m1:
        mean.append(np.mean(Dgas[i]))
        d0, d1, d2 = ErrorBarMedian(Dgas[i])
        data[0].append(d0)
        data[1].append(d1)
        data[2].append(d2)
    data = np.array(data)
    yerr = np.vstack((data[1,:] - data[0,:], data[2,:] - data[1,:]))
    x = bins + 0.026
    ax2.errorbar(x[1:], data[1,1:], yerr=yerr[:,1:], elinewidth=1.8, capthick=1.8, capsize=2, color='r', fmt='.', label='TNG Gas Fraction Errorbar')

    #il1 GasFraction ErrorBar
    Dgas = il1_gf[il1_diskID]
    data = [[], [], []]
    mean = []
    for i in m2:
        mean.append(np.mean(Dgas[i]))
        d0, d1, d2 = ErrorBarMedian(Dgas[i])
        data[0].append(d0)
        data[1].append(d1)
        data[2].append(d2)
    data = np.array(data)
    yerr = np.vstack((data[1,:] - data[0,:], data[2,:] - data[1,:]))
    x = bins + 0.093
    ax2.errorbar(x[1:], data[1,1:], yerr=yerr[:,1:], elinewidth=1.8, capthick=1.8, capsize=2, color='c', fmt='.', label='il-1 Gas Fraction Errorbar')


    ax1.legend()
    ax2.legend(loc=2)
    plt.savefig('f:/Linux/local_result/Z=2/TNG_il1_sMassAndGas_Z=2.png', dpi = 300)