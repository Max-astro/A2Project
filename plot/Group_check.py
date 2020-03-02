import numpy as np
import h5py
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import sys
import requests
import illustris_python as il
from datetime import datetime

def Progenitors_dictionary(haloID):
    with h5py.File('/Raid0/zhouzb/sublink/sublink_mpb_%d.hdf5'%haloID,'r') as f:
        snap_num = np.array(f['SnapNum'])
        subfind_ID = np.array(f['SubfindID'])

        Progenitors_dict = {}
        for i in range(len(snap_num)):
            Progenitors_dict['%d'%snap_num[i]] = subfind_ID[i]
    return Progenitors_dict

def subhaloA2(haloID, snapnum, simu='TNG'):
    data = il.func.loadgalaxy(simu, snapnum, haloID, partType=4, fields=['Coordinates','Masses','Velocities','GFM_StellarFormationTime'])
    coor = np.array(data['Coordinates'])
    mas = np.array(data['Masses'])
    vel = np.array(data['Velocities'])
    sf_time = np.array(data['GFM_StellarFormationTime'])

    half_r = (il.func.loadSubhalos(simu,snapnum,'SubhaloHalfmassRadType'))[haloID,4]
    halo_position = coor[0] + 0.00000001

    is_star = (sf_time>=0.0) # don't use wind particles

    r = coor - halo_position
    r[r > 37500] -= 75000
    r[r < -37500] += 75000
    dis = ((r**2).sum(1))**0.5
    inside = dis < (half_r*2)

    vel = vel[(inside & is_star)]
    mas = mas[(inside & is_star)]
    coor = r[(inside & is_star)]

    #Calculate angular momentum
    V = np.sum(vel * mas[:, np.newaxis], 0) / mas.sum()
    vel = vel - V
    Jz = np.sum(np.cross(coor, vel * mas[:, np.newaxis]), axis=0)

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
    A2_R = (a2_R ** 2 + b2_R ** 2)** 0.5 / a0
    return A2_R

def LoadMergHist(simu, subhaloID):
    '''
    return subhalo's main progenitor and merger history with snapshot
    '''
    if simu == 'TNG':
        ldir = '/Raid0/zhouzb/merg_data/tng_DiskMerTree/%d.json' % subhaloID
    else:
        ldir = '/Raid0/zhouzb/merg_data/il1_DiskMerTree/%d.json' % subhaloID
    
    with open(ldir) as f:
        data = json.load(f)
    
    Main = np.array(data['Main'])
    return dict(zip(Main[:, 0], Main[:, 1])), np.array(data['Mergers'])

snapshot = [91, 84, 78, 72, 67, 59, 50, 40, 33]
ToList = {'99': 0, '91' : 1, '84' : 2, '78' : 3, '72' : 4, '67' : 5, '59' : 6, '50' : 7, '40' : 8, '33' : 9}
Redshift = {'99': 0, '91' : 0.1, '84' : 0.2, '78' : 0.3, '72' : 0.4, '67' : 0.5, '59' : 0.7, '50' : 1.0, '40' : 1.5, '33' : 2.0}


tng_barredID = np.load('/Raid0/zhouzb/TNG/barredID.npy')

barFromedSnap = []
secular = 0
for haloID in tng_barredID:
    time_start=time.time()
    Main, Mergers = LoadMergHist('TNG', haloID)
    A2list = []
    for snap in SnapList[1:]:
        try:
            if Main[snap] == -1:
                A2list.append(1)
            else:
                A2 = subhaloA2(Main[snap], snap)
        except KeyError:
            A2list.append(0)
        else:
            A2 = A2[int(len(A2) / 100):]
            A2list.append(A2.max())

    haloMass = il.func.loadSubhalos('TNG', 99, 'SubhaloMass')[haloID]
    progID = []
    for i in SnapList:
        try:
            progID.append(Main[i])
        except KeyError:
            progID.append(-1)
    oriSnap = isSecular(progID, SnapList, A2list, Mergers)
    barFromedSnap.append([haloID, oriSnap[0], oriSnap[1]])
    if oriSnap == -1:
        secular += 1
    time_end = time.time()
    print('time: %.3fs'%(time_end - time_start))


#ids = np.load('/Raid0/zhouzb/TNG/barredID.npy')
for diskID in ids:
    
    #Load halo's haloID with redshift
    SnapDict = Progenitors_dictionary(diskID)
    os.system('mkdir /Raid0/zhouzb/fig_TNG/group_check/halo_%d'%diskID)
    #Load data with redshift
    half_r = (il.groupcat.loadSubhalos('/Raid1/Illustris/TNG-100',snapnum,'SubhaloHalfmassRadType'))[diskID,4]
    for snapnum in snapshot:
        try:
            haloID = SnapDict['%d'%snapnum]
        except:
            print('Halo %d progeniter in snapshot %d no found.'%(diskID, snapnum))
            continue


        #load subhalo information
        coor = il.snapshot.loadSubhalo('/Raid1/Illustris/TNG-100', snapnum, haloID, partType=4, fields='Coordinates')

        #load groups particles inside r*40
        GrNr = il.groupcat.loadSubhalos('/Raid1/Illustris/TNG-100', snapnum,'SubhaloGrNr')[haloID]
        group_particles = il.snapshot.loadHalo('/Raid1/Illustris/TNG-100', snapnum, GrNr, partType=4, fields='Coordinates')
        p0 = coor[0]
        coor -= p0
        coor[coor > 37500] -= 75000
        coor[coor < -37500] += 75000
        group_particles -= p0
        group_particles[group_particles > 37500] -= 75000
        group_particles[group_particles < -37500] += 75000
        
        boxsize =80
        halosize = 18

        inside = (group_particles[:, 0] < half_r * boxsize) & (group_particles[:, 0] > -half_r * boxsize) & (group_particles[:, 1] < half_r * boxsize) & (group_particles[:, 1] > -half_r * boxsize) & (group_particles[:, 2] < half_r * boxsize) & (group_particles[:, 2] > -half_r * boxsize)
        dis = np.linalg.norm(group_particles, axis=1)
        exclude = dis < half_r*halosize
        #exclude = (group_particles[:, 0] < half_r * 15) & (group_particles[:, 0] > -half_r * 15) & (group_particles[:, 1] < half_r * 15) & (group_particles[:, 1] > -half_r * 15) & (group_particles[:, 2] < half_r * 15) & (group_particles[:, 2] > -half_r * 15)
        
        group_particles = group_particles[inside & (~exclude)]

        coor = coor[::3]
        group_particles = group_particles[::15]

        a = datetime.now()
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(coor[:,0], coor[:,1], coor[:,2],c='r',marker='.')
        ax.scatter(group_particles[:,0], group_particles[:,1], group_particles[:,2],c='b',marker='.', alpha = 0.05, s=1)
        plt.title('halo%d in RedShift_%.1f.png'%(diskID, Redshift['%d'%snapnum]))
        
        for ang in range(60, 181,60):
            ax.view_init(azim = ang)
            fig.savefig('/Raid0/zhouzb/fig_TNG/group_check/halo_%d/z=%.1f_Theta=%d.png'%(diskID, Redshift['%d'%snapnum], ang), dpi=300)
        
        plt.close('all')
        b = datetime.now()
        print('halo in z=%.1f finished. Time: '%snapnum,(b-a).seconds)