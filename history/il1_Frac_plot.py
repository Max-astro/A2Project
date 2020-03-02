import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
import sys
sys.path.append('F:\Linux')
import illustris_python as il

def FilepathList(path, suffix='.hdf5'):
    L = []
    for files in os.listdir(path):
        if os.path.splitext(files)[1] == '%s'%suffix:
            L.append(int(os.path.splitext(files)[0]))
    return L    

# particleNum = il.groupcat.loadSubhalos('f:/Linux/data/illustris-1',135,'SubhaloLenType')[:,4]

# barID = np.load('F:/Linux/data/barID_il1.npy')
barID = FilepathList('F:/Linux/data/135_4WP_ALL','.png')
diskID = np.load('F:/Linux/data/diskID_il1.npy')
StellarMass = il.groupcat.loadSubhalos('f:/Linux/data/illustris-1',135,'SubhaloMassType')[:,4]

# minmass = 10.5
#Disk halo's mass
diskmass = StellarMass[diskID]
diskmass = np.log10(diskmass*10**10)
# diskmass = diskmass[diskmass > minmass]
#Barred halo's mass
barmass = StellarMass[barID]
barmass = np.log10(barmass*10**10)
# barmass = barmass[barmass > minmass]

#
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_xlabel('Stellar Mass')
ax1.set_ylabel('Bar Fraction')
ax2 = ax1.twinx()
ax2.set_ylabel('Halo number N')

#plot histogram
n,bins,others = ax2.hist(diskmass, 17, rwidth=0.9, alpha = 0.5)
ax2.set_xlim(10.0,12)

Fraction = []
x_point = []
for i in range(len(bins)-1):
    low = bins[i]
    high = bins[i+1]
    x_point.append((low + high)/2)

    disknum = len(diskmass[(diskmass >= low) & (diskmass < high)])
    barred = len(barmass[(barmass >= low) & (barmass < high)])
    if disknum == 0:
        Fraction.append(0)
        continue
    Barfraction = barred / disknum
    Fraction.append(Barfraction)

ax1.plot(x_point, Fraction, 'o', c = 'r')
plt.savefig('F:/Linux/data/il1_DISK-BarFraction.png',dpi=300)

