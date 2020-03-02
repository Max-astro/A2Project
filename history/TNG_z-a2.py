import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
import sys
import requests
sys.path.append('F:\Linux')
import illustris_python as il

#Calculate subhalo's A2 parameter, return a list of A2 inside the half mass radius
def subhaloA2(haloID, snapnum):
    with h5py.File('F:/Linux/data/TNG/cutoff/disk_%d/cutout_%d.hdf5'%(snapnum, haloID),'r') as f:
        coor = np.array(f['PartType4']['Coordinates'])
        mas = np.array(f['PartType4']['Masses'])
        vel = np.array(f['PartType4']['Velocities'])
        sf_time = np.array(f['PartType4']['GFM_StellarFormationTime'])

    half_r = (il.groupcat.loadSubhalos('/Raid1/Illustris/Illustris-1',snapnum,'SubhaloHalfmassRadType'))[haloID,4]
    halo_position = il.groupcat.loadSubhalos('/Raid1/Illustris/Illustris-1',snapnum,'SubhaloPos')[haloID]

    is_star = (sf_time>=0.0) # don't use wind particles


    r = coor - halo_position
    dis = ((r**2).sum(1))**0.5
    inside = dis < (half_r*2)

    vel = vel[(inside & is_star)]
    mas = mas[(inside & is_star)]
    coor = r[(inside & is_star)]

    #Calculate angular momentum
    for i in range(len(vel)):
        vel[i] = vel[i] * mas[i]
    J = np.cross(coor, vel)
    Jz = J.sum(axis = 0)

    #Set halo's angular momentum Jz as z axis
    coor = rothalo(coor, Jz)

    #r_list is a list about the distance between stellar particles and halo_position
    #r_sort a index about paritcles' distance
    r_list = ((coor[:,:2]**2).sum(1))**0.5
    r_sort = r_list.argsort()

    #a0 = Sigma(M[i])
    a0 = np.zeros(len(r_sort))
    ptr = 0
    for i in r_sort:
        if ptr == 0:
            a0[ptr] = mas[i]
        else:
            a0[ptr] = mas[i] + a0[ptr-1]
        ptr += 1

    #Position angle: Theta = arctan(y/x)
    #Creat a list of a_m(R)
    a2_R = np.zeros(len(r_sort))
    b2_R = np.zeros(len(r_sort))
    list_i = 0
    for star in r_sort:
        Theta = np.arctan(coor[star,0] / coor[star,1])
        #a_i = M[numb] * cos( m * Theta[numb] )  , m=2
        a_i = mas[star] * np.cos(2*Theta)
        #b_i = M[numb] * sin( m * Theta[numb] ) , m=2
        b_i = mas[star] * np.sin(2*Theta) 

        #a2_R[i] = a_i.sum(:i)
        if list_i > 0:
            a2_R[list_i] = a_i + a2_R[list_i - 1]
            b2_R[list_i] = b_i + b2_R[list_i - 1]
        else:
            a2_R[list_i] = a_i
            b2_R[list_i] = b_i

        list_i += 1
        
    #It's a list about A2 parameter, sorted by distance between halo and center
    A2_R = (a2_R**2 + b2_R**2)**0.5 / a0
    return A2_R


def Progenitors_dictionary(haloID):
    with h5py.File('F:/Linux/data/TNG/cutoff/sub_mpb/halo_%d.npy'%diskID,'r') as f:
        snap_num = np.array(f['SnapNum'])
        subfind_ID = np.array(f['SubfindID'])

        Progenitors_dict = {}
        for i in range(len(snap_num)):
            Progenitors_dict['%d'%snap_num[i]] = subfind_ID[i]
    return Progenitors_dict


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


StellarParticleLen = il.groupcat.loadSubhalos('F:/Linux/data/TNG/TNG-100', 99, 'SubhaloLenType')[:,4]
StellarParticleLen = np.pad(StellarParticleLen, (0, 5127294-len(StellarParticleLen)), 'constant')
for snap_num in snapshot:
    group = il.groupcat.loadSubhalos('F:/Linux/data/TNG/TNG-100', snap_num, 'SubhaloLenType')[:,4]
    StellarParticleLen = np.vstack((StellarParticleLen, np.pad(group, (0, abs(5127294-len(group))), 'constant')))

ids = np.load('F:/Linux/data/099fig/barredID.npy')

#Painting
for diskID in ids:
    A2list = []
    halolen = []
    snaplist = []

    plt.figure()
    plt.xlabel('Red Shift')
    plt.ylabel('A2')

    #Load A2 parameter inside the stellar half mass radius*2, and pick A2.max() as x point
    tmp = subhaloA2(diskID,)
    A2 = np.array(tmp[0])
    A2 = A2[int(len(A2) / 100):]
    if A2.max() < 0.15:
        continue

    A2list.append(A2.max())
    snaplist.append(0)


    #Load halo's haloID with redshift
    SnapDict = Progenitors_dictionary(diskID)
    #Load data with redshift
    for snap_num in snapshot:
        try:
            haloID = SnapDict['%d'%snap_num]
            tmp = np.load('F:/Linux/data/TNG/cutoff/disk/disk_%d/%d.npy'%(snap_num, haloID))
            snaplist.append(Redshift['%d'%snap_num])
        except:
            print('Halo %d no found in snapshot %d'%(diskID, snap_num))
            continue

        A2 = np.array(tmp[0])
        A2 = A2[int(len(A2) / 100):]
        A2list.append(A2.max())

        #if particle number less than 40000, use carmine and yellow
        halolen.append(StellarParticleLen[ToList['%d'%snap_num], haloID])



    #Normal Color: blue = StellarMass, red = GasFraction
    #Small halo color: cyan = StellarMass, carmine = GasFraction
    #if particle number less than 40000, use carmine and yellow
    mask = np.array(halolen) < 40000

    #Plot StellarMass in ax1 with black point, and GasFraction in ax2 with blue point
    # red = ax1.plot(x_point, mass_point, 'ob', color = 'r', label = 'StellarMass')
    plt.scatter(snaplist, A2list, marker='o', color='b')
    plt.plot(snaplist, A2list)
    if mask.any() != False:
        A2list = np.array(A2list[1:])
        snaplist = np.array(snaplist[1:])


        # carmine = ax1.plot(x_point[mask], mass_point[mask], 'ob', color = 'c', label = 'StellarMass, P<40000')
        yellow = plt.plot(snaplist[mask], A2list[mask], 'ob', color = 'y', label = 'Stellar particles <40000')
        plt.legend()
    plt.xlim(-0.1, 1.1)
    plt.ylim(0, 0.6)
    plt.savefig('/Raid0/zhouzb/fig_TNG/z-a2/%d.png'%diskID)
    print('%d saved'%diskID)

