import numpy as np
import h5py
import os
import illustris_python as il
import matplotlib.pyplot as plt


# snap_num = 99
diskID = np.load('F:/Linux/data/diskID.npy')
StellarMass = il.groupcat.loadSubhalos('F:/Linux/data/TNG/Groupcatalog', 99, 'SubhaloMassType')[:,4]
#load barred galaxies' ID
bigID = np.load('F:/Linux/data/bigID.npy')
smallID = np.load('F:/Linux/data/smallID.npy')
ids = np.concatenate((smallID, bigID))

#Barred halo's mass
halomass = StellarMass[ids]
halomass = np.log10(halomass*10**10)
diskmass = StellarMass[diskID]
StellarMass = np.log10(diskmass*10**10)

#Create figer
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_xlabel('Stellar Mass')
ax1.set_ylabel('Bar Fraction')
ax2 = ax1.twinx()
ax2.set_ylabel('Halo number N')

#plot histogram
n,bins,others = ax2.hist(halomass, 20, rwidth=0.9)
ax2.set_xlim(9.8,12)

Fraction = []
x_point = []
for i in range(20):
    low = bins[i]
    high = bins[i+1]
    x_point.append((low + high)/2)

    disknum = len(diskmass[(diskmass >= low) & (diskmass < high)])
    barred = len(halomass[(halomass >= low) & (halomass < high)])

    Barfraction = barred / disknum
    Fraction.append(Barfraction)

ax1.plot(x_point, Fraction, 'o', c = 'r')




