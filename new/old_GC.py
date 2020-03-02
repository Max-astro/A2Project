import numpy as np
import h5py
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import json
import os
import sys
import requests
import illustris_python as il
from datetime import datetime

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


snapshot = [91, 84, 78, 72, 67, 59, 50]

ToList = {'99': 0, '91' : 1, '84' : 2, '78' : 3, '72' : 4, '67' : 5, '59' : 6, '50' : 7, '40' : 8, '33' : 9}
Redshift = {'99': 0, '91' : 0.1, '84' : 0.2, '78' : 0.3, '72' : 0.4, '67' : 0.5, '59' : 0.7, '50' : 1.0, '40' : 1.5, '33' : 2.0}

newfig = []
ids = np.load('/Raid0/zhouzb/TNG/barredID_4WP_TNG.npy')
done=0
for diskID in ids[692:693]:

    #Load halo's haloID with redshift
    Main, Mergers = LoadMergHist('TNG', diskID)
    s = os.system('mkdir /Raid0/zhouzb/fig_TNG/group_check/halo_%d'%diskID)
    # if s != 0:
    #     continue

    #Load data with redshift
    for snapnum in snapshot:
        newfig.append(snapnum)
        try:
            haloID = Main[snapnum]
        except:
            print('Halo %d progeniter in snapshot %d no found.'%(diskID, snapnum))
            continue

        if haloID == -1:
            continue
        #load subhalo information
        coor = il.snapshot.loadSubhalo('/Raid1/Illustris/TNG-100', snapnum, haloID, partType=4, fields='Coordinates')

        half_r = (il.groupcat.loadSubhalos('/Raid1/Illustris/TNG-100',snapnum,'SubhaloHalfmassRadType'))[haloID,4]

        #load groups particles inside r*40
        GrNr = il.groupcat.loadSubhalos('/Raid1/Illustris/TNG-100', snapnum,'SubhaloGrNr')[haloID]
        group_particles = il.snapshot.loadHalo('/Raid1/Illustris/TNG-100', snapnum, GrNr, partType=4, fields='Coordinates')

        p0 = coor[0]

        group_particles -= p0
        group_particles[group_particles > 37500] -= 75000
        group_particles[group_particles < -37500] += 75000

        coor -= p0
        coor[coor > 37500] -= 75000
        coor[coor < -37500] += 75000

        boxsize = 80
        halosize = 15

        inside = (group_particles[:, 0] < half_r * boxsize) & (group_particles[:, 0] > -half_r * boxsize) & (group_particles[:, 1] < half_r * boxsize) & (group_particles[:, 1] > -half_r * boxsize) & (group_particles[:, 2] < half_r * boxsize) & (group_particles[:, 2] > -half_r * boxsize)
        dis = np.linalg.norm(group_particles, axis=1)
        exclude = dis < half_r*halosize
        #exclude = (group_particles[:, 0] < half_r * 15) & (group_particles[:, 0] > -half_r * 15) & (group_particles[:, 1] < half_r * 15) & (group_particles[:, 1] > -half_r * 15) & (group_particles[:, 2] < half_r * 15) & (group_particles[:, 2] > -half_r * 15)

        group_particles = group_particles[inside & (~exclude)]

        coor = coor[::3]
        group_particles = group_particles[::15]

        print(group_particles.shape)

        a = datetime.now()
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(coor[:,0], coor[:,1], coor[:,2],c='r',marker='.', s=3)
        ax.scatter(group_particles[:,0], group_particles[:,1], group_particles[:,2],c='b',marker='.', alpha=0.5, s=0.05)
        plt.title('halo%d in RedShift_%.1f.png'%(diskID, Redshift['%d'%snapnum]))

        for ang in range(0,361,120):
            ax.view_init(azim = ang)
            fig.savefig('/Raid0/zhouzb/fig_TNG/group_check/halo_%d/z=%.1f_Theta=%d.png'%(diskID, Redshift['%d'%snapnum], ang), dpi=300)
        plt.close('all')
        b = datetime.now()
        print('halo in z=%.1f finished. Time: '%snapnum,(b-a).seconds)
    done += 1
    print(done)

np.save('/Raid0/zhouzb/newfig.npy',newfig)
