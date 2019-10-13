import numpy as np
import h5py
import os
import illustris_python as il
import matplotlib.pyplot as plt

def FilepathList(path, suffix='.hdf5'):
    L = []
    for files in os.listdir(path):
        if os.path.splitext(files)[1] == '%s'%suffix:
            L.append(int((os.path.splitext(files)[0])[5:]))
    return L    
    
path = '/Raid0/zhouzb/diskHalo_Sublink/'
ids = FilepathList(path, suffix = '.npy')

'''
snap_91 z=0.1
snap_84 z=0.2
snap_78 z=0.3
snap_72 z=0.4
snap_67 z=0.5
snap_59 z=0.7
snap_50 z=1.0
snap_40 z=1.5
snap_33 z=2.0
'''
snapshot = [91, 84, 78, 72, 67, 59, 50]

ToList = {'99': 0, '91' : 1, '84' : 2, '78' : 3, '72' : 4, '67' : 5, '59' : 6, '50' : 7}
Redshift = {'99': 0, '91' : 0.1, '84' : 0.2, '78' : 0.3, '72' : 0.4, '67' : 0.5, '59' : 0.7, '50' : 1.0}
#Redshift = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0]
'''
GasFraction / StellarMass[0] = snap_99
GasFraction / StellarMass[1] = 91,
GasFraction / StellarMass[2] = 84,
GasFraction / StellarMass[3] = 78,
GasFraction / StellarMass[4] = 72,
GasFraction / StellarMass[5] = 67,
GasFraction / StellarMass[6] = 59,
GasFraction / StellarMass[7] = 50
'''


StellarMass = il.groupcat.loadSubhalos('/Raid1/Illustris/TNG-100', 99, 'SubhaloMassType')[:,4]
StellarMass = np.pad(StellarMass, (0, 5127294-len(StellarMass)), 'constant')
for snap_num in snapshot:
    group = il.groupcat.loadSubhalos('/Raid1/Illustris/TNG-100', snap_num, 'SubhaloMassType')[:,4]
    StellarMass = np.vstack((StellarMass, np.pad(group, (0, abs(5127294-len(group))), 'constant')))

GasFraction = il.groupcat.loadSubhalos('/Raid1/Illustris/TNG-100', 99, 'SubhaloMassType')[:,0]
GasFraction = np.pad(GasFraction, (0, 5127294-len(GasFraction)), 'constant')
for snap_num in snapshot:
    group = il.groupcat.loadSubhalos('/Raid1/Illustris/TNG-100', snap_num, 'SubhaloMassType')[:,0]
    GasFraction = np.vstack((GasFraction, np.pad(group, (0, abs(5127294-len(group))), 'constant')))

StellarParticleLen = il.groupcat.loadSubhalos('/Raid1/Illustris/TNG-100', 99, 'SubhaloLenType')[:,4]
StellarParticleLen = np.pad(StellarParticleLen, (0, 5127294-len(StellarParticleLen)), 'constant')
for snap_num in snapshot:
    group = il.groupcat.loadSubhalos('/Raid1/Illustris/TNG-100', snap_num, 'SubhaloLenType')[:,4]
    StellarParticleLen = np.vstack((StellarParticleLen, np.pad(group, (0, abs(5127294-len(group))), 'constant')))

#Painting
for diskID in ids[:25]:
    
    x_point = []
    mass_point = []
    gas_point = []
    halolen = []
    snaplist = []

    
    
    plt.xlabel('A2')
    plt.ylabel('Stellar Mass')

    #Load A2 parameter inside the stellar half mass radius*2, and pick A2.max() as x point
    tmp = np.load('/Raid0/zhouzb/TNG_a2/disk/disk_99/%d.npy'%diskID)
    A2 = np.array(tmp[0])
    A2 = A2[int(len(A2) / 100):]
    x_point.append(A2.max())

    tmp_mas = StellarMass[0, diskID] * 10**10
    tmp_gas = GasFraction[0, diskID] * 10**10
    mass_point.append(np.log10(tmp_mas))
    gas_point.append(tmp_gas / (tmp_mas + tmp_gas))
    
    plt.annotate('z=0', xy=(A2.max(), np.log10(tmp_mas)), xycoords='data', xytext=(-20, +20), textcoords='offset points', arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))
   

    SnapDict = (np.load('/Raid0/zhouzb/diskHalo_Sublink/halo_%d.npy'%diskID)).item()
    #Load data with redshift
    for snap_num in snapshot:
        try:
            haloID = SnapDict['%d'%snap_num]
            tmp = np.load('/Raid0/zhouzb/TNG_a2/disk/disk_%d/%d.npy'%(snap_num, haloID))
            snaplist.append(snap_num)
        except:
            print('Halo %d no found in snapshot %d'%(diskID, snap_num))
            continue

        A2 = np.array(tmp[0])
        A2 = A2[int(len(A2) / 100):]
        x_point.append(A2.max())

        tmp_mas = StellarMass[ToList['%d'%snap_num], haloID] * 10**10
        tmp_gas = GasFraction[ToList['%d'%snap_num], haloID] * 10**10
        mass_point.append(np.log10(tmp_mas))
        gas_point.append(tmp_gas / (tmp_mas + tmp_gas))

        #if particle number less than 40000, use carmine and yellow
        halolen.append(StellarParticleLen[ToList['%d'%snap_num], haloID])

        # ax1.annotate('z=%.1f'%Redshift['%d'%snap_num], xy=(A2.max(), tmp_mas), xycoords='data', xytext=(+8, +18), textcoords='offset points', arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))
        # plt.annotate('z=%.1f'%Redshift['%d'%snap_num], xy=(A2.max(), tmp_mas), xycoords='data', xytext=(+8, -18), textcoords='offset points', arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))


    #Normal Color: blue = StellarMass, red = GasFraction
    #Small halo color: cyan = StellarMass, carmine = GasFraction
    #if particle number less than 40000, use carmine and yellow
    # mask = np.array(halolen) < 40000

    #Plot StellarMass in ax1 with black point, and GasFraction in ax2 with blue point
    plt.scatter(x_point, mass_point, c = gas_point, label = 'StellarMass')
   
    # ax1.plot(x_point, mass_point)
    plt.plot(x_point, mass_point)
    plt.colorbar()

    plt.legend()  
    
    #Set xy axis' limit

    plt.ylim(9.5, 11.8)
    plt.xlim(0, 0.6)

    plt.savefig('/Raid0/zhouzb/fig_TNG/disk_Evolve/%d.png'%diskID)
    print('%d saved'%diskID)
    plt.close()

