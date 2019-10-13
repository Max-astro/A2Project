import numpy as np
import h5py
import json
import sys
import csv
import illustris_python as il
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def SelectDisk(simu, Pmin, Pmax):
    if (simu == 'TNG') | (simu ==  'tng'):
        filePath = '/Raid0/zhouzb/TNGdata/supplementary/stellar_circs.hdf5'
        snap_num = 99
    else:
        filePath = '/Raid0/zhouzb/il1_data/stellar_circs.hdf5'
        snap_num = 135
    
    with h5py.File(filePath,'r') as cir:
        haloID = np.array(cir['Snapshot_%d'%snap_num]['SubfindID'])
        cir07frac = np.array(cir['Snapshot_%d'%snap_num]['CircAbove07Frac'])
        MassTensor = np.array(cir['Snapshot_%d'%snap_num]['MassTensorEigenVals'])

    #select halo stellar Pmax > particles > Pmin
    halolen = il.func.loadSubhalos(simu, snap_num, 'SubhaloLenType')[:, 4]
    #circularity parameter ϵ > 0.2
    cir_mask = cir07frac > 0.2
    
    #flatness of the galaxy is defined as the ratio M1=(M2M3)**0.5 , disk galaxy's flatness < 0.7
    MassTensor = MassTensor[cir_mask]
    haloID = haloID[cir_mask]

    flat = MassTensor[:,0]/(MassTensor[:,1]*MassTensor[:,2])**0.5
    flat_mask = flat < 0.7

    haloID = haloID[flat_mask]

    mas_mask = (halolen[haloID] > Pmin) & (halolen[haloID] < Pmax )
    haloID = haloID[mas_mask]
    return haloID

def rothalo(coor, Jz):
    '''
    set halo's angular momentum as z axis
    '''
    rot = RotMatrix(Norm(Jz))

    for i in range(len(coor)):
        coor[i] = np.dot(rot, coor[i])    

    return coor

#Vector Normalized
def Norm(vect):
    '''Input a vector, normalize it in-place'''
    mol=np.linalg.norm(vect)
    vect /= mol
    return vect

#Find the rotation matrix when input vector rotated to z_axis
def RotMatrix(vect):
    '''Input a vector, return a (3,3) rotation matrix '''
    x=vect[0]
    y=vect[1]
    z=vect[2]

    sinA= y/(z**2+y**2)**0.5
    cosA= z/(z**2+y**2)**0.5
    sinB= x/(x**2+(y*sinA+z*cosA)**2)**0.5
    cosB= (y*sinA+z*cosA)/(x**2+(y*sinA+z*cosA)**2)**0.5

    return np.array([[cosB, -sinA*sinB, -cosA*sinB], [0, cosA, -sinA], [sinB, sinA*cosB, cosA*cosB]])

def subhaloA2(simu, snapnum, haloID):
    if haloID == -1:
        return -1
    try:
        data = il.func.loadgalaxy(simu, snapnum, haloID, partType=4, fields=['Coordinates', 'Masses', 'Velocities', 'GFM_StellarFormationTime'])
        coor= data['Coordinates']
        mas = data['Masses']
        vel = data['Velocities']
        sf_time = data['GFM_StellarFormationTime']
    except:
        print(sys.exc_info()[0])
        print('halo %d in snap_%d load faild.'%(haloID,snapnum))
        print(' ')
        return -1


    half_r = (il.func.loadSubhalos(simu, snapnum, 'SubhaloHalfmassRadType'))[haloID,4]
    halo_position = il.func.loadSubhalos(simu, snapnum, 'SubhaloPos')[haloID]

    is_star = (sf_time>=0.0) # don't use wind particles

    r = coor - halo_position
    dis = ((r**2).sum(1))**0.5
    inside = dis < (half_r*2)

    vel = vel[(inside & is_star)]
    mas = mas[(inside & is_star)]
    coor = r[(inside & is_star)]

    coor[coor > 37500] -= 75000
    coor[coor < -37500] += 75000
    #Calculate angular momentum
    V = np.sum(vel * mas[:, np.newaxis], 0) / mas.sum()
    vel = vel - V
    L = np.sum(np.cross(coor, vel * mas[:, np.newaxis]), axis=0)

    #Set halo's angular momentum Jz as z axis
    coor = rothalo(coor, L)

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
        if coor[star, 1] == 0:
            continue
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
    A2_R = A2_R[int(len(A2_R) / 100):]
    return A2_R.max()

def plotBarFig(simu, snapnum, haloID)：
    data = il.func.loadgalaxy(simu, snapnum, haloID, partType=4, fields=['Coordinates', 'Masses', 'Velocities', 'GFM_StellarFormationTime'])
    coor= data['Coordinates']
    mas = data['Masses']
    vel = data['Velocities']
    sf_time = data['GFM_StellarFormationTime']
    half_r = (il.func.loadSubhalos(simu, snapnum, 'SubhaloHalfmassRadType'))[haloID,4]
    is_star = (sf_time>=0.0) # don't use wind particles

    r = coor - coor[0]
    dis = ((r**2).sum(1))**0.5
    inside = dis < (half_r*6)

    vel = vel[(inside & is_star)]
    mas = mas[(inside & is_star)]
    coor = r[(inside & is_star)]
    coor[coor > 37500] -= 75000
    coor[coor < -37500] += 75000
    #Calculate angular momentum
    V = np.sum(vel * mas[:, np.newaxis], 0) / mas.sum()
    vel = vel - V
    L = np.sum(np.cross(coor, vel * mas[:, np.newaxis]), axis=0)
    #Set halo's angular momentum Jz as z axis
    coor = rothalo(coor, L)
    #plot
    plt.rcParams['figure.figsize'] = (36, 36)
    plt.tick_params(labelsize=70)
    plt.hist2d(coor[::2, 0], coor[::2, 1], bins=500, norm=LogNorm(0.035), cmap='jet')
    if simu == 'TNG':
        savedir = '/Raid0/zhouzb/barfig/small/tng/%d.png' % haloID
    else:
        savedir = '/Raid0/zhouzb/barfig/small/il1/%d.png' % haloID
    plt.savefig(savedir)
    plt.close()

il1_smDisk = SelectDisk('il1', 20000, 40000)
tng_smDisk = SelectDisk('TNG', 20000, 40000)

il1fig = 0
for haloID in il1_smDisk:
    if subhaloA2(simu='il1', snapnum=135, haloID=haloID) > 0.15:
        plotBarFig(simu='il1', snapnum=135, haloID=haloID)
        il1fig += 1
        print('il1 small bar : %d / %d' % (il1fig, len(il1_smDisk)))
        
tngfig = 0
for haloID in tng_smDisk:
    if subhaloA2(simu='TNG', snapnum=99, haloID=haloID) > 0.15:
        plotBarFig(simu='TNG', snapnum=99, haloID=haloID)
        tngfig += 1
        print('TNG small bar : %d / %d'%(tngfig, len(il1_smDisk)))
print('All Done')