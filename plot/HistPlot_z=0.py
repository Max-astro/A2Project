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
        return 0, np.median(data), 0
    else:
        data.sort()
        return data[int(len(data) / 4)], np.median(data), data[int(len(data) * 0.75)]


#TNG data
barID = np.load('f:/Linux/localRUN/barredID_4WP_TNG.npy')
diskID = np.load('f:/Linux/localRUN/diskID_4WP.npy')
#Gas Fraction Gf
mas = il.func.loadSubhalos('TNG', 99, 'SubhaloMassInRadType')
Gf = mas[:, 0] / (mas[:, 4] + mas[:, 0])
Gf[np.isnan(Gf)] = 0
#Stellar Particles
SP = il.func.loadSubhalos('TNG', 99, 'SubhaloLenType')[:, 4]
#Stellar Mass
sMass = il.func.loadSubhalos('TNG', 99, 'SubhaloMassType')[:, 4] / 0.6774
sMass = np.log10(sMass * 10 ** 10)
sMass[np.isinf(sMass)] = 0

#Illsutris-1 data
il1_barID = np.load('f:/Linux/localRUN/barredID_il1.npy')
il1_diskID = np.load('f:/Linux/localRUN/diskID_il1.npy')
il1_mas = il.func.loadSubhalos('il1', 135, 'SubhaloMassInRadType')
#Gas Fraction
il1_gf = il1_mas[:, 0] / (il1_mas[:, 4] + il1_mas[:, 0])
il1_gf[np.isnan(il1_gf)] = 0
#Stellar Particles
il1_SP = il.func.loadSubhalos('il1', 135, 'SubhaloLenType')[:, 4]
#Stellar Mass
il1_sMass = il.func.loadSubhalos('il1', 135, 'SubhaloMassType')[:, 4] / 0.704
il1_sMass = np.log10(il1_sMass * 10 ** 10)
il1_sMass[np.isinf(il1_sMass)] = 0

def plot_fig_1():
    #Fig : 'TNG-4WP_GasFraction.png'
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('GasFraction')
    ax1.set_ylabel('N')
    ax1.set_title('TNG Gas Fraction & Bar Fraction')
    ax2 = ax1.twinx()
    ax2.set_ylabel('Bar Fraction')
    ax2.set_ylim(0, 0.8)
    

    bins = np.linspace(0, 0.45, 9)
    n1 = HistValAndBin(Gf[diskID], bins)
    n2 = HistValAndBin(Gf[barID], bins)
    ax1.bar(bins[:-1], n1, width=(bins[1] - bins[0])*0.9,align = 'edge', label='TNG disk galaxies')
    ax1.bar(bins[:-1], n2, width=(bins[1] - bins[0])*0.9,align = 'edge', label='TNG barred galaxies')
    ax2.plot(bins[:-1] + 0.025, n2 / n1, marker='o', color='r')
    ax1.legend()
    plt.savefig('f:/Linux/local_result/tng-4WP_GasFraction.png', dpi = 300)

def plot_fig_2():
    #Fig : 'TNG_Illustris-1_GF_Z=0.png'
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('GasFraction')
    ax1.set_ylabel('N')
    ax1.set_title('TNG & Illustris-1 Gas Fraction & Bar Fraction')
    ax1.set_ylim(0, 1000)
    ax2 = ax1.twinx()
    ax2.set_ylabel('Total Bar Fraction')
    ax2.set_ylim(0, 1.1)

    bins = np.linspace(0, 0.6, 12)
    disk_tng_n = HistValAndBin(Gf[diskID], bins)
    bar_tng_n = HistValAndBin(Gf[barID], bins)
    disk_il1_n = HistValAndBin(il1_gf[il1_diskID], bins)
    bar_il1_n = HistValAndBin(il1_gf[il1_barID], bins)
    ax1.bar(bins[:-1], disk_tng_n, width=(bins[1] - bins[0]) * 0.35,align = 'edge', label='TNG disk galaxies')
    ax1.bar(bins[:-1], bar_tng_n, width=(bins[1] - bins[0]) * 0.35, align='edge', label='TNG barred galaxies')
    ax1.bar(bins[:-1] + 0.02, disk_il1_n, width=(bins[1] - bins[0]) * 0.35, align = 'edge', label='Illustris-1 disk galaxies', color='c')
    ax1.bar(bins[:-1] + 0.02, bar_il1_n, width=(bins[1] - bins[0]) * 0.35, align='edge', label='Illustris-1 barred galaxies', color='r')

    frac = bar_tng_n / disk_tng_n
    frac[-3:] = 0
    ax2.plot(bins[:-1] + 0.021, frac, marker='o', label='TNG', color='b')
    ax2.plot(bins[:-1] + 0.021, bar_il1_n / disk_il1_n, marker='o', label='Illsutis-1', color='k')

    ax1.legend()
    ax2.legend(loc=2)
    plt.savefig('f:/Linux/local_result/TNG_Illustris-1_GF_Z=0.png', dpi=300)
    
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
    plt.savefig('f:/Linux/local_result/TNG_Illustris-1_GF_Z=0.png', dpi=300)

'''
#Galaxies with more than 100k stellar particles
diskID_10WP = []
barID_10WP = []
for haloID in diskID:
    if SP[haloID] >= 100000:
        diskID_10WP.append(haloID)
        if haloID in barID:
            barID_10WP.append(haloID)
il1_diskID_10WP = []
il1_barID_10WP = []
for haloID in il1_diskID:
    if il1_SP[haloID] >= 100000:
        il1_diskID_10WP.append(haloID)
        if haloID in il1_barID:
            il1_barID_10WP.append(haloID)
diskID_10WP = np.array(diskID_10WP)
barID_10WP = np.array(barID_10WP)
il1_diskID_10WP = np.array(il1_diskID_10WP)
il1_barID_10WP = np.array(il1_barID_10WP)
'''

def plot_sMass():
    #Campare TNG and il1 galaxies distribution in bins
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Stellar Mass')
    ax1.set_ylabel('N')
    ax1.set_title('Stellar Mass distribution between TNG and Illustirs-1')
    ax1.set_xlim(10.4, 12)
    ax1.set_ylim(0, 550)
    ax2 = ax1.twinx()
    ax2.set_ylabel('Times')
    ax2.set_ylim(0, 2.5)

    y = np.ones(100)
    x = np.linspace(10, 15, 100)
    ax2.plot(x, y, c='b', linestyle='-.')

    bins = np.linspace(10.2, 11.6, 20)
    n1 = HistValAndBin(sMass[diskID], bins, more=1)
    n2 = HistValAndBin(il1_sMass[il1_diskID], bins, more=1)
    ax1.bar(bins, n1+n2, width=(bins[1] - bins[0])*0.9,align = 'edge', label='TNG disk galaxies')
    ax1.bar(bins, n2, width=(bins[1] - bins[0])*0.9,align = 'edge', label='il1 disk galaxies')
    ax2.plot(bins + (bins[1]-bins[0])/2, n1 / n2, marker='o', color='r', label='Number TNG / Il-1')
    ax1.legend()
    ax2.legend(loc=2)
    plt.savefig('f:/Linux/local_result/TNG_Illustris-1_MassDistribution.png', dpi = 300)


def Plot_TNG_sMass_BarFraction():
    #
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Stellar Mass')
    ax1.set_ylabel('N')
    ax1.set_title('TNG Bar Fraction with Stellar Mass bin')
    ax1.set_xlim(10.4, 12)
    ax1.set_ylim(0, 400)
    ax2 = ax1.twinx()
    ax2.set_ylabel('Fraction')
    ax2.set_ylim(0, 1)

    bins = np.linspace(10.2, 11.7, 20)
    n1 = HistValAndBin(sMass[diskID], bins, more=1)
    n2 = HistValAndBin(sMass[barID], bins, more=1)
    ax1.bar(bins, n1, width=(bins[1] - bins[0]) * 0.9, align='edge', label='disk nobar')
    ax1.bar(bins, n2, width=(bins[1] - bins[0]) * 0.9, align='edge', label='barred galaxies')
    ax1.legend()

    ax2.plot(bins + (bins[1] - bins[0]) / 2, n2 / n1, marker='o', color='r', label='Bar Fraction')
    ax2.legend(loc=2)
    plt.savefig('f:/Linux/local_result/TNG_StellarMass.png', dpi = 300)

def Plot_il1_sMass_BarFraction():
    #
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Stellar Mass')
    ax1.set_ylabel('N')
    ax1.set_title('il1_ Bar Fraction with Stellar Mass bin')
    ax1.set_xlim(10.4, 12)
    ax1.set_ylim(0, 400)
    ax2 = ax1.twinx()
    ax2.set_ylabel('Fraction')
    ax2.set_ylim(0, 1)

    bins = np.linspace(10.2, 11.7, 20)
    n1 = HistValAndBin(il1_sMass[il1_diskID], bins, more=1)
    n2 = HistValAndBin(il1_sMass[il1_barID], bins, more=1)
    ax1.bar(bins, n1, width=(bins[1] - bins[0]) * 0.9, align='edge', label='disk nobar')
    ax1.bar(bins, n2, width=(bins[1] - bins[0]) * 0.9, align='edge', label='barred galaxies')
    ax1.legend()

    ax2.plot(bins + (bins[1] - bins[0]) / 2, n2 / n1, marker='o', color='r', label='Bar Fraction')
    ax2.legend(loc=2)
    plt.savefig('f:/Linux/local_result/il1_StellarMass.V2.png', dpi = 300)

def Plot_TNG_sMassAndGasFraction():
    #
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Stellar Mass')
    ax1.set_ylabel('N')
    ax1.set_title('TNG disk galaxies Gas Fraction in different Stellar Mass bin')
    ax1.set_xlim(10.4, 12)
    ax1.set_ylim(0, 400)
    ax2 = ax1.twinx()
    ax2.set_ylabel('Gas Fraction')
    ax2.set_ylim(0, 0.3)

    bins = np.linspace(10.2, 11.7, 20)
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
    ax2.errorbar(x, data[1,:], yerr=yerr, elinewidth=2, capthick=2, capsize=3, color='r', fmt='o', label='Median Gas Fraction Errorbar')
    ax2.plot(x, mean, color='c', marker='o', label='Mean Gas Fraction')
    ax2.legend(loc=2)
    plt.savefig('f:/Linux/local_result/TNG_sMassAndGas_Mean.png', dpi = 300)


def Plot_il1_sMassAndGasFraction():
    #
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Stellar Mass')
    ax1.set_ylabel('N')
    ax1.set_title('Illustris-1 disk galaxies Gas Fraction in different Stellar Mass bin')
    ax1.set_xlim(10.4, 12)
    ax1.set_ylim(0, 400)
    ax2 = ax1.twinx()
    ax2.set_ylabel('Gas Fraction')
    ax2.set_ylim(0, 0.6)

    bins = np.linspace(10.2, 11.7, 20)
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
    plt.savefig('f:/Linux/local_result/Illustris_sMassAndGas_Mean.png', dpi = 300)

